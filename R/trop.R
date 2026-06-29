# Triply RObust Panel (TROP) estimator.
# Athey, Imbens, Qu & Viviano (2025), arXiv:2508.21536.
#
# Working model:  Y_it(0) = alpha_i + beta_t + L_it + eps_it
# Estimation (eq. 2):
#   minimise  sum_{j,s} theta_s omega_j (1 - W_js) (Y_js - a_j - b_s - L_js)^2
#             + lambda_nn ||L||_*
#   then      tau_it = Y_it - a_i - b_t - L_it
# Weights (eq. 3):
#   theta_s  = exp(-lambda_time * |t - s|)
#   omega_j  = exp(-lambda_unit * dist_unit(j, i))
# Tuning (eq. 4-5): leave-one-out CV on control cells, choosing lambda to
#   minimise the squared "as-if-treated" prediction error on controls.

# ---- weight building blocks -------------------------------------------------

#' Unit distance to a single anchor unit (eq. 3)
#' @keywords internal
#' @noRd
.unit_distance_to <- function(Y, W, i, t = NULL) {
  N <- nrow(Y); Tt <- ncol(Y)
  ci <- (W[i, ] == 0)
  if (!is.null(t)) ci[t] <- FALSE
  yi <- Y[i, ]
  bi <- ci & is.finite(yi)                       # usable columns from unit i
  if (!any(bi)) return(rep(NA_real_, N))
  # vectorised over control units j (same result as the old per-j loop)
  useM <- (W == 0) & is.finite(Y) & matrix(bi, N, Tt, byrow = TRUE)
  diff2 <- (Y - matrix(yi, N, Tt, byrow = TRUE))^2
  diff2[!useM] <- 0
  n <- rowSums(useM); ssq <- rowSums(diff2)
  d <- rep(NA_real_, N); pos <- n > 0
  d[pos] <- sqrt(ssq[pos] / n[pos])
  d
}

#' Unit distances anchored to a set of treated units (pooled weighting)
#' @keywords internal
#' @noRd
.unit_distance_pooled <- function(Y, W, treated_units) {
  ds <- vapply(treated_units, function(i) .unit_distance_to(Y, W, i, NULL),
               numeric(nrow(Y)))
  if (is.null(dim(ds))) ds <- matrix(ds, ncol = 1)
  rowMeans(ds, na.rm = TRUE)
}

#' Build the N x T weight matrix theta_s * omega_j for an anchor.
#' @keywords internal
#' @noRd
.trop_weight_matrix <- function(du, t_anchor, Tt, lam) {
  du[!is.finite(du)] <- if (any(is.finite(du))) max(du[is.finite(du)]) else 0
  omega <- exp(-lam$unit * du)
  # time distance: min distance to the anchor period(s)
  s <- seq_len(Tt)
  if (length(t_anchor) == 1L) {
    dtime <- abs(s - t_anchor)
  } else {
    dtime <- vapply(s, function(ss) min(abs(ss - t_anchor)), numeric(1))
  }
  theta <- exp(-lam$time * dtime)
  outer(omega, theta)
}

# ---- single-cell counterfactual prediction ---------------------------------

#' Predict the control counterfactual for cells, with anchored weights.
#'
#' Solves eq. (2) once for the supplied weights and returns the fitted matrix.
#' @keywords internal
#' @noRd
.trop_solve <- function(Y, W, wmat, lam, ctrl, L_init = NULL, X = NULL) {
  mask <- (W == 0) * 1
  .mcnnm_fit(Y, mask, wmat, lam$nn,
             max_iter = ctrl$max_iter, tol = ctrl$tol,
             L_init = L_init, svd_method = ctrl$svd %||% "truncated", X = X)
}

# ---- leave-one-out cross-validation (eq. 4-5) ------------------------------

#' LOO-CV criterion Q(lambda) on a subsample of control cells.
#'
#' For each sampled control cell (i, t) we treat it as if it were treated
#' (mask it out of the loss), build weights anchored to that cell, solve eq. (2),
#' and accumulate the squared prediction error -- a subsampled realisation of
#' eq. (5).
#' @keywords internal
#' @noRd
.trop_cv_Q <- function(Y, W, lam, ctrl, cv_cells, du_list = NULL, X = NULL) {
  Tt <- ncol(Y)
  one <- function(k) {
    i <- cv_cells[k, 1]; t <- cv_cells[k, 2]
    Wk <- W
    Wk[i, t] <- 1                       # hold this control cell out
    # unit distances do not depend on the penalties: reuse a cached value when
    # one is supplied (see .trop_select_lambda), else compute it here.
    du <- if (is.null(du_list)) .unit_distance_to(Y, Wk, i, t) else du_list[[k]]
    wmat <- .trop_weight_matrix(du, t, Tt, lam)
    fit <- .trop_solve(Y, Wk, wmat, lam, ctrl, X = X)
    (Y[i, t] - fit$M[i, t])^2
  }
  par <- (ctrl$workers %||% 1L) > 1L
  errs <- unlist(.par_lapply(seq_len(nrow(cv_cells)), one, parallel = par),
                 use.names = FALSE)
  mean(errs)
}

#' Coordinate-descent search over (lambda_time, lambda_unit, lambda_nn).
#'
#' Follows the warm-start scheme of footnote 2 of the paper: start from
#' lambda_nn = Inf, lambda_unit = 0, optimise each penalty marginally, then
#' cycle.
#' @keywords internal
#' @noRd
.trop_select_lambda <- function(Y, W, grids, ctrl, cv_cells, verbose = FALSE,
                                X = NULL) {
  lam <- list(time = 0, unit = 0, nn = Inf)
  # Unit distances for each held-out CV cell do not depend on the penalties, so
  # compute them once and reuse across the whole penalty search (big saving:
  # the search evaluates the CV criterion dozens of times).
  du_cache <- lapply(seq_len(nrow(cv_cells)), function(k) {
    i <- cv_cells[k, 1]; t <- cv_cells[k, 2]
    Wk <- W; Wk[i, t] <- 1
    .unit_distance_to(Y, Wk, i, t)
  })
  eval_one <- function(lam) .trop_cv_Q(Y, W, lam, ctrl, cv_cells,
                                       du_list = du_cache, X = X)

  pick <- function(lam, which, grid) {
    qs <- vapply(grid, function(v) {
      l2 <- lam; l2[[which]] <- v
      eval_one(l2)
    }, numeric(1))
    best <- grid[which.min(qs)]
    if (verbose) {
      message(sprintf("  %-4s -> %.4g (Q=%.5g)", which, best, min(qs)))
    }
    best
  }

  # marginal initialisation
  lam$time <- pick(lam, "time", grids$time)
  lam$nn   <- pick(lam, "nn",   grids$nn)
  lam$unit <- pick(lam, "unit", grids$unit)

  # cycles
  for (cyc in seq_len(ctrl$cv_cycles)) {
    lam$time <- pick(lam, "time", grids$time)
    lam$unit <- pick(lam, "unit", grids$unit)
    lam$nn   <- pick(lam, "nn",   grids$nn)
  }
  lam$Q <- eval_one(lam)
  lam
}

#' Default penalty grids derived from the data scale.
#' @keywords internal
#' @noRd
.trop_default_grids <- function(Y, W, n_unit = 6L, n_time = 7L, X = NULL) {
  # nuclear-norm grid from singular values of the demeaned control matrix
  Yc <- Y
  Yc[W == 1 | is.na(Yc)] <- NA
  fill <- mean(Yc, na.rm = TRUE)
  Z <- Yc
  Z[is.na(Z)] <- fill
  Z <- Z - rowMeans(Z)
  Z <- t(t(Z) - colMeans(Z))
  # With covariates the nuclear norm penalises the residual R = L - X phi, so
  # scale the grid off the covariate-residualised matrix (partial the covariates
  # out of the already two-way-demeaned Z).
  Xl <- .as_cov_list(X, nrow(Y), ncol(Y))
  if (length(Xl)) {
    Xdd <- matrix(unlist(lapply(Xl, function(M) as.numeric(.double_demean(M))),
                         use.names = FALSE), ncol = length(Xl))
    cf <- qr.coef(qr(Xdd), as.numeric(Z)); cf[is.na(cf)] <- 0
    Z <- matrix(as.numeric(Z) - as.numeric(Xdd %*% cf), nrow(Y), ncol(Y))
  }
  s1 <- tryCatch(max(svd(Z)$d), error = function(e) 1)
  nn_grid <- c(Inf, s1 * c(0.5, 0.25, 0.1, 0.05, 0.02))

  # unit grid scaled by the median pairwise distance among control units
  tu <- which(rowSums(W) > 0)
  anchor <- if (length(tu)) tu[1] else 1
  du <- .unit_distance_to(Y, W, anchor, NULL)
  med <- stats::median(du[is.finite(du) & du > 0], na.rm = TRUE)
  if (!is.finite(med) || med <= 0) med <- 1
  unit_grid <- c(0, c(0.25, 0.5, 1, 2, 4)[seq_len(n_unit - 1)]) / med

  # time grid (decay per period)
  time_grid <- c(0, 0.02, 0.05, 0.1, 0.25, 0.5, 1)[seq_len(n_time)]

  list(time = time_grid, unit = unit_grid, nn = nn_grid)
}

# ---- main estimator ---------------------------------------------------------

#' Triply RObust Panel (TROP) estimator
#'
#' Fits the TROP estimator of Athey, Imbens, Qu & Viviano (2025) on a long
#' panel. TROP combines a low-rank-plus-two-way-fixed-effects outcome model with
#' exponential-decay unit weights (upweighting controls similar to the treated)
#' and time weights (upweighting periods near the treated periods). Penalties are
#' chosen by leave-one-out cross-validation on the control cells. The estimator
#' nests DID/TWFE, matrix completion and synthetic-control-type weighting as
#' special cases.
#'
#' @param data A long `data.frame` with one row per unit-time.
#' @param outcome,treatment,unit,time Column names (strings). `treatment` must be
#'   a 0/1 indicator of *active* treatment in that cell (works for block,
#'   staggered, and non-absorbing designs).
#' @param covariates Optional character vector of column names for time-varying
#'   covariates. When supplied, the low-rank term is augmented additively as
#'   \eqn{L_{js} = X_{js}'\phi + R_{js}} (Athey, Imbens, Qu & Viviano 2025,
#'   Section 6.2): the covariate-linear part is unpenalised and the nuclear norm
#'   applies to the residual \eqn{R}. Penalty cross-validation is performed on the
#'   covariate-residualised model. Covariates must be fully observed (including in
#'   treated cells, where they enter the counterfactual). The fitted coefficients
#'   are returned in `$phi`.
#' @param lambda Optional named list `list(time=, unit=, nn=)` to fix the
#'   penalties and skip cross-validation.
#' @param anchor How weights are anchored to treated cells: `"per_cell"`
#'   re-solves eq. (2) with cell-specific weights for every treated cell
#'   (faithful to the paper); `"pooled"` solves once with weights anchored to the
#'   treated set (fast); `"auto"` (default) uses `per_cell` when there are at
#'   most `max_cells` treated cells and `pooled` otherwise.
#' @param se Standard-error method: `"bootstrap"` (default; unit-level
#'   stratified block bootstrap), `"jackknife"` (leave-one-treated-unit-out;
#'   needs at least 2 treated units), `"auto"` (jackknife when there are at least
#'   2 treated units, else bootstrap), or `"none"`.
#' @param grids Optional list of penalty grids; see [trop_control()].
#' @param control A list of solver/CV settings from [trop_control()].
#' @param verbose Logical; print CV progress.
#' @return An object of class `trop`: a list with the ATT `estimate`,
#'   `std.error`, `conf.low`/`conf.high`, the selected penalties `lambda`,
#'   any covariate coefficients `phi`, per-cell effects, the estimated
#'   counterfactual matrix, weights, and the reshaped panel.
#' @references Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2025).
#'   Triply Robust Panel Estimators. arXiv:2508.21536.
#' @examples
#' df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, seed = 1)
#' fit <- trop(df, "y", "w", "id", "t", se = "none",
#'             control = trop_control(n_cv_cells = 8L, cv_cycles = 1L))
#' fit
#' @export
trop <- function(data, outcome, treatment, unit, time,
                 covariates = NULL,
                 lambda = NULL,
                 anchor = c("auto", "per_cell", "pooled"),
                 se = c("bootstrap", "auto", "jackknife", "none"),
                 grids = NULL,
                 control = trop_control(),
                 verbose = FALSE) {
  anchor <- match.arg(anchor)
  se <- match.arg(se)
  m <- .panel_to_matrices(data, outcome, treatment, unit, time)
  Y <- m$Y; W <- m$W
  pat <- .assignment_pattern(W)
  X <- .covariate_matrices(data, covariates, unit, time, m$units, m$times)

  n_cells <- pat$n_treated_cells
  if (anchor == "auto") {
    anchor <- if (n_cells <= control$max_cells) "per_cell" else "pooled"
  }

  if (is.null(grids)) grids <- .trop_default_grids(Y, W, X = X)

  eng <- .trop_engine(Y, W, pat, lambda, grids, control, anchor, se,
                      control$conf_level, verbose, X = X)
  if (length(eng$phi)) eng$phi <- stats::setNames(eng$phi, covariates)
  m$X <- X

  structure(
    c(eng, list(
      conf.level = control$conf_level,
      pattern = pat,
      panel = m,
      outcome = outcome,
      covariates = covariates,
      call = match.call()
    )),
    class = "trop"
  )
}

#' Shared TROP engine operating directly on Y / W matrices.
#'
#' Used by [trop()] and by the native DID/MC engines in [panel_compare()] so
#' that penalty selection and inference are identical across methods.
#' @keywords internal
#' @noRd
.trop_engine <- function(Y, W, pat, lambda, grids, control, anchor, se,
                         conf_level, verbose = FALSE, X = NULL) {
 .with_workers(control$workers %||% 1L, {
  # ---- choose penalties via CV (unless supplied) ----
  if (is.null(lambda)) {
    cv_cells <- .sample_control_cells(W, control$n_cv_cells, control$seed)
    lam <- .trop_select_lambda(Y, W, grids, control, cv_cells, verbose, X = X)
  } else {
    lam <- utils::modifyList(list(time = 0, unit = 0, nn = Inf), lambda)
    lam$Q <- NA_real_
  }

  est <- .trop_att(Y, W, lam, control, anchor, pat, X = X)

  if (identical(se, "auto")) {
    se <- if (pat$n_treated_units >= 2) "jackknife" else "bootstrap"
  }
  inf <- .trop_se(Y, W, lam, control, anchor, pat, est, se, conf_level, X = X)

  list(
    estimate = est$att,
    std.error = inf$se,
    conf.low = inf$conf.low,
    conf.high = inf$conf.high,
    se.method = se,
    n_boot = inf$n_boot,
    lambda = lam,
    anchor = anchor,
    tau_cells = est$tau_cells,
    counterfactual = est$Mhat,
    rank = est$rank,
    phi = est$phi
  )
 })
}

#' @keywords internal
#' @noRd
.sample_control_cells <- function(W, n, seed = NULL) {
  ctrl_idx <- which(W == 0, arr.ind = TRUE)
  if (nrow(ctrl_idx) > n) {
    if (!is.null(seed)) {
      old <- .Random.seed_safe()
      on.exit(.Random.seed_restore(old), add = TRUE)
      set.seed(seed)
    }
    ctrl_idx <- ctrl_idx[sample.int(nrow(ctrl_idx), n), , drop = FALSE]
  }
  ctrl_idx
}

#' @keywords internal
#' @noRd
.Random.seed_safe <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else NULL
}
#' @keywords internal
#' @noRd
.Random.seed_restore <- function(old) {
  if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv)
}

#' ATT given fixed penalties.
#' @keywords internal
#' @noRd
.trop_att <- function(Y, W, lam, ctrl, anchor, pat, X = NULL) {
  N <- nrow(Y); Tt <- ncol(Y)
  treated_cells <- which(W == 1, arr.ind = TRUE)
  Mhat <- matrix(NA_real_, N, Tt, dimnames = dimnames(Y))
  rnk <- NA_integer_
  phi <- numeric(0)

  if (anchor == "pooled") {
    du <- .unit_distance_pooled(Y, W, pat$treated_units)
    t_anchor <- sort(unique(treated_cells[, 2]))
    wmat <- .trop_weight_matrix(du, t_anchor, Tt, lam)
    fit <- .trop_solve(Y, W, wmat, lam, ctrl, X = X)
    Mhat <- fit$M
    rnk <- fit$rank
    phi <- fit$phi
  } else {
    ranks <- integer(nrow(treated_cells))
    phis <- vector("list", nrow(treated_cells))
    L_prev <- NULL                      # warm start: reuse the previous cell's L
    for (k in seq_len(nrow(treated_cells))) {
      i <- treated_cells[k, 1]; t <- treated_cells[k, 2]
      du <- .unit_distance_to(Y, W, i, t)
      wmat <- .trop_weight_matrix(du, t, Tt, lam)
      fit <- .trop_solve(Y, W, wmat, lam, ctrl, L_init = L_prev, X = X)
      Mhat[i, t] <- fit$M[i, t]
      ranks[k] <- fit$rank
      phis[[k]] <- fit$phi
      L_prev <- fit$L                   # consecutive cells share dimensions
    }
    rnk <- round(mean(ranks))
    # per-cell anchor re-solves once per treated cell; report the average
    # covariate coefficient across cells as a summary (the cell-specific phi is
    # already baked into each cell's counterfactual above).
    if (length(phis) && length(phis[[1]]))
      phi <- rowMeans(matrix(unlist(phis, use.names = FALSE),
                             nrow = length(phis[[1]])))
  }

  tau <- Y[treated_cells] - Mhat[treated_cells]
  tau_cells <- data.frame(
    unit = rownames(Y)[treated_cells[, 1]],
    time = colnames(Y)[treated_cells[, 2]],
    y = Y[treated_cells],
    y0_hat = Mhat[treated_cells],
    tau = tau,
    stringsAsFactors = FALSE
  )
  list(att = mean(tau), tau_cells = tau_cells, Mhat = Mhat, rank = rnk,
       phi = phi)
}

#' Standard errors for the TROP ATT (jackknife / bootstrap).
#' @keywords internal
#' @noRd
.trop_se <- function(Y, W, lam, ctrl, anchor, pat, est, method, conf_level,
                     X = NULL) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  na_out <- list(se = NA_real_, conf.low = NA_real_, conf.high = NA_real_)
  if (method == "none") return(na_out)
  sub_rows <- function(idx) if (is.null(X)) NULL else
    lapply(.as_cov_list(X, nrow(Y), ncol(Y)),
           function(M) M[idx, , drop = FALSE])

  if (method == "jackknife") {
    tu <- pat$treated_units
    G <- length(tu)
    if (G < 2) return(na_out)
    par <- (ctrl$workers %||% 1L) > 1L
    atts <- unlist(.par_lapply(seq_len(G), function(g) {
      drop <- tu[g]
      keep <- setdiff(seq_len(nrow(Y)), drop)
      Yk <- Y[keep, , drop = FALSE]
      Wk <- W[keep, , drop = FALSE]
      patk <- .assignment_pattern(Wk)
      .trop_att(Yk, Wk, lam, ctrl, anchor, patk, X = sub_rows(keep))$att
    }, parallel = par), use.names = FALSE)
    v <- (G - 1) / G * sum((atts - mean(atts))^2)
    se <- sqrt(v)
  } else if (method == "bootstrap") {
    # Unit-level stratified block bootstrap (Athey, Imbens, Qu & Viviano 2026):
    # resample treated and control units separately, with replacement, keeping
    # each unit's full time series, and recompute the ATT.
    tu <- pat$treated_units
    co <- setdiff(seq_len(nrow(Y)), tu)
    if (length(tu) < 1 || length(co) < 2) return(na_out)
    B <- ctrl$n_boot %||% 200L
    if (!is.null(ctrl$seed)) {
      old <- .Random.seed_safe(); on.exit(.Random.seed_restore(old), add = TRUE)
      set.seed(ctrl$seed)
    }
    # Draw all resample index sets up front (reproducible given the seed); the
    # ATT solves below are deterministic, so they can be evaluated in parallel
    # without affecting the draws.
    idx_list <- lapply(seq_len(B), function(b)
      c(sample(tu, length(tu), replace = TRUE),
        sample(co, length(co), replace = TRUE)))
    par <- (ctrl$workers %||% 1L) > 1L
    one_boot <- function(idx) {
      Yb <- Y[idx, , drop = FALSE]; Wb <- W[idx, , drop = FALSE]
      rownames(Yb) <- rownames(Wb) <- as.character(seq_along(idx))
      patb <- .assignment_pattern(Wb)
      if (patb$n_treated_units < 1 || length(setdiff(seq_len(nrow(Yb)),
          patb$treated_units)) < 1) return(NA_real_)
      tryCatch(.trop_att(Yb, Wb, lam, ctrl, anchor, patb,
                         X = sub_rows(idx))$att,
               error = function(e) NA_real_)
    }
    boot <- unlist(.par_lapply(idx_list, one_boot, parallel = par),
                   use.names = FALSE)
    boot <- boot[is.finite(boot)]
    if (length(boot) < 2) return(na_out)
    se <- stats::sd(boot)
    if (identical(ctrl$boot_ci %||% "percentile", "percentile")) {
      q <- stats::quantile(boot, c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2),
                           names = FALSE, type = 7)
      return(list(se = se, conf.low = q[1], conf.high = q[2],
                  n_boot = length(boot)))
    }
    return(list(se = se, conf.low = est$att - z * se,
                conf.high = est$att + z * se, n_boot = length(boot)))
  } else {
    return(na_out)
  }
  list(se = se,
       conf.low = est$att - z * se,
       conf.high = est$att + z * se)
}

#' Control settings for [trop()]
#'
#' @param max_iter Maximum solver iterations.
#' @param tol Solver convergence tolerance.
#' @param n_cv_cells Number of control cells sampled for the CV criterion.
#' @param cv_cycles Number of coordinate-descent cycles in penalty selection.
#' @param max_cells Threshold for `anchor = "auto"` to switch to pooled weights.
#' @param conf_level Confidence level for intervals.
#' @param n_boot Number of replications for the bootstrap standard error
#'   (`se = "bootstrap"`).
#' @param boot_ci Bootstrap confidence-interval type: `"percentile"` or
#'   `"normal"`.
#' @param svd Singular-value decomposition used by the soft-impute solver:
#'   `"truncated"` (default) computes only the leading singular triplets with
#'   `RSpectra` when it is installed and the matrix is large enough to benefit,
#'   falling back to the full SVD otherwise; `"full"` always uses the exact base
#'   R `svd()`. The two agree to numerical tolerance; `"truncated"` is faster on
#'   large panels, `"full"` is used for exact numerical-agreement checks.
#' @param workers Number of parallel workers for the embarrassingly parallel
#'   loops (cross-validation cells, and the bootstrap / jackknife
#'   replicates). `1` (default) runs serially. Values `> 1` use
#'   `future.apply`/`future` when installed (a transient `multisession` plan is
#'   set up and restored automatically); if those packages are missing it falls
#'   back to serial with no error. Results are reproducible given `seed`.
#' @param seed Optional seed for CV cell sampling (reproducibility).
#' @return A list of control parameters.
#' @export
trop_control <- function(max_iter = 200L, tol = 1e-5,
                         n_cv_cells = 120L, cv_cycles = 2L,
                         max_cells = 60L, conf_level = 0.95,
                         n_boot = 200L, boot_ci = c("percentile", "normal"),
                         svd = c("truncated", "full"),
                         workers = 1L,
                         seed = NULL) {
  boot_ci <- match.arg(boot_ci)
  svd <- match.arg(svd)
  list(max_iter = max_iter, tol = tol, n_cv_cells = n_cv_cells,
       cv_cycles = cv_cycles, max_cells = max_cells,
       conf_level = conf_level, n_boot = n_boot, boot_ci = boot_ci,
       svd = svd, workers = workers,
       seed = seed)
}

#' @export
print.trop <- function(x, ...) {
  cat("Triply RObust Panel (TROP) estimator\n")
  cat(sprintf("  ATT = %.4f", x$estimate))
  if (is.finite(x$std.error)) {
    cat(sprintf("  (SE %.4f, %s)\n", x$std.error, x$se.method))
    cat(sprintf("  %.0f%% CI: [%.4f, %.4f]\n",
                100 * x$conf.level, x$conf.low, x$conf.high))
  } else cat("\n")
  cat(sprintf("  penalties: lambda_time=%.3g, lambda_unit=%.3g, lambda_nn=%s\n",
              x$lambda$time, x$lambda$unit,
              ifelse(is.infinite(x$lambda$nn), "Inf",
                     sprintf("%.3g", x$lambda$nn))))
  cat(sprintf("  estimated rank(L)=%s, anchor=%s, design=%s\n",
              x$rank, x$anchor, x$pattern$type))
  cat(sprintf("  treated cells=%d across %d unit(s)\n",
              x$pattern$n_treated_cells, x$pattern$n_treated_units))
  if (length(x$phi)) {
    nm <- names(x$phi) %||% paste0("X", seq_along(x$phi))
    cat("  covariate coefficients (phi):\n")
    cat(paste0(sprintf("    %-12s % .4g", nm, x$phi), collapse = "\n"), "\n")
  }
  invisible(x)
}

# ---- per-period event study -------------------------------------------------

#' Pooled TROP counterfactual matrix (all cells).
#'
#' The per-cell anchor leaves the counterfactual undefined off the treated
#' cells, but the event study needs fitted values for the pre-treatment periods
#' too. This solves eq. (2) once with weights anchored to the whole treated set,
#' exactly as [autoplot.trop()] does for its trajectory plot, and returns the
#' full N x T fitted matrix.
#' @keywords internal
#' @noRd
.trop_pooled_M <- function(Y, W, lam, ctrl, pat, X = NULL) {
  tu <- pat$treated_units
  du <- .unit_distance_pooled(Y, W, tu)
  t_anchor <- sort(unique(which(W == 1, arr.ind = TRUE)[, 2]))
  wmat <- .trop_weight_matrix(du, t_anchor, ncol(Y), lam)
  .trop_solve(Y, W, wmat, lam, ctrl, X = X)$M
}

#' Average treatment effect by event time (period relative to treatment).
#'
#' Event time is the column-index difference t - t0_i, where t0_i is the first
#' treated period of unit i; this handles irregular spacing and staggered
#' adoption naturally. Returns the mean gap Y - Mhat and the contributing cell
#' count for each event time. With `pre_periods = FALSE` only event times >= 0
#' are kept (otherwise the pre-treatment gaps become placebo / pre-trend
#' points). A pooled `M` may be supplied to avoid re-solving.
#' @keywords internal
#' @noRd
.trop_period_effects <- function(Y, W, lam, ctrl, pat, pre_periods = TRUE,
                                 M = NULL, X = NULL) {
  tu <- pat$treated_units
  if (length(tu) < 1) return(list(effect = numeric(0), n = integer(0)))
  Tt <- ncol(Y)
  if (is.null(M)) M <- .trop_pooled_M(Y, W, lam, ctrl, pat, X = X)
  gap <- Y - M
  t0 <- vapply(tu, function(i) {
    w <- which(W[i, ] == 1); if (length(w)) min(w) else NA_integer_
  }, integer(1))
  ok <- is.finite(t0); tu <- tu[ok]; t0 <- t0[ok]
  if (!length(tu)) return(list(effect = numeric(0), n = integer(0)))
  E <- outer(t0, seq_len(Tt), function(a, b) b - a)
  G <- gap[tu, , drop = FALSE]
  keep <- is.finite(G) & is.finite(E)
  if (!pre_periods) keep <- keep & (E >= 0L)
  if (!any(keep)) return(list(effect = numeric(0), n = integer(0)))
  eff <- tapply(G[keep], E[keep], mean)
  cnt <- tapply(G[keep], E[keep], length)
  list(effect = stats::setNames(as.numeric(eff), names(eff)),
       n = stats::setNames(as.integer(cnt), names(cnt)))
}

#' Per-event-time standard errors (jackknife / bootstrap).
#'
#' Mirrors `.trop_se()` but, instead of collapsing each resample to a single
#' ATT, records the vector of period effects and aggregates it by event time.
#' The expensive per-resample refit is identical; only the aggregation differs.
#' Intervals are pointwise.
#' @keywords internal
#' @noRd
.trop_event_se <- function(Y, W, lam, ctrl, pat, point, method, conf_level,
                           pre_periods, X = NULL) {
  ev <- names(point)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  na_vec <- function() stats::setNames(rep(NA_real_, length(ev)), ev)
  empty <- function() list(se = na_vec(), conf.low = na_vec(),
                           conf.high = na_vec(), n_boot = NULL)
  if (method == "none") return(empty())
  par <- (ctrl$workers %||% 1L) > 1L
  sub_rows <- function(idx) if (is.null(X)) NULL else
    lapply(.as_cov_list(X, nrow(Y), ncol(Y)),
           function(M) M[idx, , drop = FALSE])
  to_row <- function(v) { r <- v[ev]; names(r) <- ev; r }
  col_sd <- function(M) apply(M, 2, function(x) {
    x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else stats::sd(x)
  })

  if (method == "jackknife") {
    tu <- pat$treated_units; G <- length(tu)
    if (G < 2) return(empty())
    reps <- .par_lapply(seq_len(G), function(g) {
      keep <- setdiff(seq_len(nrow(Y)), tu[g])
      Yk <- Y[keep, , drop = FALSE]; Wk <- W[keep, , drop = FALSE]
      patk <- .assignment_pattern(Wk)
      to_row(.trop_period_effects(Yk, Wk, lam, ctrl, patk, pre_periods,
                                  X = sub_rows(keep))$effect)
    }, parallel = par)
    M <- do.call(rbind, reps)
    v <- vapply(seq_len(ncol(M)), function(j) {
      x <- M[, j]; x <- x[is.finite(x)]; g <- length(x)
      if (g < 2) NA_real_ else (g - 1) / g * sum((x - mean(x))^2)
    }, numeric(1))
    se <- sqrt(v)
    return(list(se = stats::setNames(se, ev),
                conf.low = stats::setNames(point - z * se, ev),
                conf.high = stats::setNames(point + z * se, ev),
                n_boot = NULL))
  }

  if (method == "bootstrap") {
    tu <- pat$treated_units; co <- setdiff(seq_len(nrow(Y)), tu)
    if (length(tu) < 1 || length(co) < 2) return(empty())
    B <- ctrl$n_boot %||% 200L
    if (!is.null(ctrl$seed)) {
      old <- .Random.seed_safe(); on.exit(.Random.seed_restore(old), add = TRUE)
      set.seed(ctrl$seed)
    }
    idx_list <- lapply(seq_len(B), function(b)
      c(sample(tu, length(tu), replace = TRUE),
        sample(co, length(co), replace = TRUE)))
    one <- function(idx) {
      Yb <- Y[idx, , drop = FALSE]; Wb <- W[idx, , drop = FALSE]
      rownames(Yb) <- rownames(Wb) <- as.character(seq_along(idx))
      patb <- .assignment_pattern(Wb)
      if (patb$n_treated_units < 1 || length(setdiff(seq_len(nrow(Yb)),
          patb$treated_units)) < 1) return(na_vec())
      tryCatch(
        to_row(.trop_period_effects(Yb, Wb, lam, ctrl, patb, pre_periods,
                                    X = sub_rows(idx))$effect),
        error = function(e) na_vec())
    }
    reps <- .par_lapply(idx_list, one, parallel = par)
    M <- do.call(rbind, reps)
    nfin <- colSums(is.finite(M)); se <- col_sd(M)
    if (identical(ctrl$boot_ci %||% "percentile", "percentile")) {
      qs <- vapply(seq_len(ncol(M)), function(j) {
        x <- M[, j]; x <- x[is.finite(x)]
        if (length(x) < 2) c(NA_real_, NA_real_)
        else stats::quantile(x, c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2),
                             names = FALSE, type = 7)
      }, numeric(2))
      lo <- qs[1, ]; hi <- qs[2, ]
    } else { lo <- point - z * se; hi <- point + z * se }
    return(list(se = stats::setNames(se, ev),
                conf.low = stats::setNames(lo, ev),
                conf.high = stats::setNames(hi, ev),
                n_boot = max(nfin)))
  }

  # any other method (e.g. "none") -> no pointwise standard errors
  return(empty())
}

#' Per-period (event-study) effects from a TROP fit
#'
#' Decomposes a fitted [trop()] ATT into per-period treatment effects indexed by
#' event time (periods relative to each unit's first treated period) and attaches
#' pointwise standard errors and confidence intervals. The same resampling used
#' for the overall SE is reused: each resample's per-cell effects are aggregated
#' by event time rather than into a single mean, so no extra modelling
#' assumptions are introduced. With `pre_periods = TRUE` the pre-treatment gaps
#' are returned as placebo / pre-trend points (a flat, near-zero pre-period
#' profile supports the design).
#'
#' Because the treated composition varies across bootstrap resamples and each
#' event time draws on relatively few cells, the per-period intervals are wider
#' than the overall ATT interval and cover each point individually (they are not
#' simultaneous bands).
#'
#' @param object A `trop` fit from [trop()].
#' @param se Standard-error method for the per-period effects: `"bootstrap"`
#'   (default; unit-level stratified block bootstrap), `"jackknife"`
#'   (leave-one-treated-unit-out; needs at least 2 treated units), or `"none"`.
#' @param pre_periods Logical; include pre-treatment event times as placebo /
#'   pre-trend points (default `TRUE`).
#' @param control A list of solver/CV/bootstrap settings from [trop_control()];
#'   defaults to `trop_control()`. Use it to set `n_boot`, `boot_ci`, `seed`, and
#'   `workers` for the resampling.
#' @param ... Unused.
#' @return An object of class `trop_event_study`: a list whose `estimates` is a
#'   `data.frame` with one row per event time (`event_time`, `estimate`,
#'   `std.error`, `conf.low`, `conf.high`, `n_cells`, `period`), plus the overall
#'   `att`, the `se.method`, and metadata.
#' @references Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2025).
#'   Triply Robust Panel Estimators. arXiv:2508.21536.
#' @seealso [autoplot.trop_event_study()] to plot the result.
#' @examples
#' \donttest{
#' df  <- sim_panel(N = 24, T = 12, n_treated = 6, t0 = 8, att = 3, seed = 1)
#' fit <- trop(df, "y", "w", "id", "t",
#'             control = trop_control(n_cv_cells = 10L, cv_cycles = 1L))
#' es  <- trop_event_study(fit, se = "bootstrap", pre_periods = TRUE,
#'                         control = trop_control(n_boot = 100L, seed = 1))
#' es
#' }
#' @export
trop_event_study <- function(object,
                             se = c("bootstrap", "jackknife", "none"),
                             pre_periods = TRUE, control = NULL, ...) {
  stopifnot(inherits(object, "trop"))
  se <- match.arg(se)
  cl <- match.call()
  if (is.null(control)) control <- trop_control()
  conf_level <- object$conf.level %||% control$conf_level
  Y <- object$panel$Y; W <- object$panel$W
  pat <- object$pattern; lam <- object$lambda
  X <- object$panel$X
  if (pat$n_treated_units < 1) stop("No treated units.", call. = FALSE)

  .with_workers(control$workers %||% 1L, {
    M <- .trop_pooled_M(Y, W, lam, control, pat, X = X)
    pe <- .trop_period_effects(Y, W, lam, control, pat, pre_periods, M = M,
                               X = X)
    point <- pe$effect
    if (!length(point)) stop("No event-time cells.", call. = FALSE)
    evt <- as.integer(names(point)); ord <- order(evt)
    point <- point[ord]; evt <- evt[ord]
    inf <- .trop_event_se(Y, W, lam, control, pat, point, se, conf_level,
                          pre_periods, X = X)
    res <- data.frame(
      event_time = evt,
      estimate   = as.numeric(point),
      std.error  = as.numeric(inf$se[as.character(evt)]),
      conf.low   = as.numeric(inf$conf.low[as.character(evt)]),
      conf.high  = as.numeric(inf$conf.high[as.character(evt)]),
      n_cells    = as.integer(pe$n[as.character(evt)]),
      period     = ifelse(evt < 0L, "pre", "post"),
      stringsAsFactors = FALSE)
    structure(
      list(estimates = res, att = object$estimate, se.method = se,
           n_boot = inf$n_boot, conf.level = conf_level,
           pre_periods = pre_periods, lambda = lam, pattern = pat,
           outcome = object$outcome, call = cl),
      class = "trop_event_study")
  })
}

#' @export
print.trop_event_study <- function(x, digits = 4, ...) {
  cat("TROP event study (per-period ATT)\n")
  cat(sprintf("  overall ATT = %.4f\n", x$att))
  est <- x$estimates
  n_post <- sum(est$period == "post")
  n_pre  <- sum(est$period == "pre")
  cat(sprintf("  event times: %d post-treatment", n_post))
  if (n_pre > 0) cat(sprintf(", %d pre-treatment (placebo)", n_pre))
  cat("\n")
  se_lab <- switch(x$se.method %||% "none",
    bootstrap = if (!is.null(x$n_boot))
                  sprintf("bootstrap (%d reps)", x$n_boot) else "bootstrap",
    jackknife = "jackknife", none = "none",
    x$se.method)
  cat(sprintf("  SE: %s; %.0f%% pointwise CI\n", se_lab, 100 * x$conf.level))
  print(format(est, digits = digits), row.names = FALSE)
  invisible(x)
}
