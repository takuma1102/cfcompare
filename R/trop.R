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
.trop_solve <- function(Y, W, wmat, lam, ctrl, L_init = NULL) {
  mask <- (W == 0) * 1
  .mcnnm_fit(Y, mask, wmat, lam$nn,
             max_iter = ctrl$max_iter, tol = ctrl$tol,
             L_init = L_init, svd_method = ctrl$svd %||% "truncated")
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
.trop_cv_Q <- function(Y, W, lam, ctrl, cv_cells, du_list = NULL) {
  Tt <- ncol(Y)
  one <- function(k) {
    i <- cv_cells[k, 1]; t <- cv_cells[k, 2]
    Wk <- W
    Wk[i, t] <- 1                       # hold this control cell out
    # unit distances do not depend on the penalties: reuse a cached value when
    # one is supplied (see .trop_select_lambda), else compute it here.
    du <- if (is.null(du_list)) .unit_distance_to(Y, Wk, i, t) else du_list[[k]]
    wmat <- .trop_weight_matrix(du, t, Tt, lam)
    fit <- .trop_solve(Y, Wk, wmat, lam, ctrl)
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
.trop_select_lambda <- function(Y, W, grids, ctrl, cv_cells, verbose = FALSE) {
  lam <- list(time = 0, unit = 0, nn = Inf)
  # Unit distances for each held-out CV cell do not depend on the penalties, so
  # compute them once and reuse across the whole penalty search (big saving:
  # the search evaluates the CV criterion dozens of times).
  du_cache <- lapply(seq_len(nrow(cv_cells)), function(k) {
    i <- cv_cells[k, 1]; t <- cv_cells[k, 2]
    Wk <- W; Wk[i, t] <- 1
    .unit_distance_to(Y, Wk, i, t)
  })
  eval_one <- function(lam) .trop_cv_Q(Y, W, lam, ctrl, cv_cells, du_list = du_cache)

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
.trop_default_grids <- function(Y, W, n_unit = 6L, n_time = 7L) {
  # nuclear-norm grid from singular values of the demeaned control matrix
  Yc <- Y
  Yc[W == 1 | is.na(Yc)] <- NA
  fill <- mean(Yc, na.rm = TRUE)
  Z <- Yc
  Z[is.na(Z)] <- fill
  Z <- Z - rowMeans(Z)
  Z <- t(t(Z) - colMeans(Z))
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
#' @param lambda Optional named list `list(time=, unit=, nn=)` to fix the
#'   penalties and skip cross-validation.
#' @param anchor How weights are anchored to treated cells: `"per_cell"`
#'   re-solves eq. (2) with cell-specific weights for every treated cell
#'   (faithful to the paper); `"pooled"` solves once with weights anchored to the
#'   treated set (fast); `"auto"` (default) uses `per_cell` when there are at
#'   most `max_cells` treated cells and `pooled` otherwise.
#' @param se Standard-error method: `"auto"`, `"jackknife"` (leave-one-treated-
#'   unit-out; needs >= 2 treated units), `"placebo"` (assign the treated
#'   pattern to controls; for a single treated unit), or `"none"`.
#' @param grids Optional list of penalty grids; see [trop_control()].
#' @param control A list of solver/CV settings from [trop_control()].
#' @param verbose Logical; print CV progress.
#' @return An object of class `trop`: a list with the ATT `estimate`,
#'   `std.error`, `conf.low`/`conf.high`, the selected penalties `lambda`,
#'   per-cell effects, the estimated counterfactual matrix, weights, and the
#'   reshaped panel.
#' @references Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2025).
#'   Triply Robust Panel Estimators. arXiv:2508.21536.
#' @examples
#' df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, seed = 1)
#' fit <- trop(df, "y", "w", "id", "t", se = "none",
#'             control = trop_control(n_cv_cells = 8L, cv_cycles = 1L))
#' fit
#' @export
trop <- function(data, outcome, treatment, unit, time,
                 lambda = NULL,
                 anchor = c("auto", "per_cell", "pooled"),
                 se = c("auto", "jackknife", "bootstrap", "placebo", "none"),
                 grids = NULL,
                 control = trop_control(),
                 verbose = FALSE) {
  anchor <- match.arg(anchor)
  se <- match.arg(se)
  m <- .panel_to_matrices(data, outcome, treatment, unit, time)
  Y <- m$Y; W <- m$W
  pat <- .assignment_pattern(W)

  n_cells <- pat$n_treated_cells
  if (anchor == "auto") {
    anchor <- if (n_cells <= control$max_cells) "per_cell" else "pooled"
  }

  if (is.null(grids)) grids <- .trop_default_grids(Y, W)

  eng <- .trop_engine(Y, W, pat, lambda, grids, control, anchor, se,
                      control$conf_level, verbose)

  structure(
    c(eng, list(
      conf.level = control$conf_level,
      pattern = pat,
      panel = m,
      outcome = outcome,
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
                         conf_level, verbose = FALSE) {
 .with_workers(control$workers %||% 1L, {
  # ---- choose penalties via CV (unless supplied) ----
  if (is.null(lambda)) {
    cv_cells <- .sample_control_cells(W, control$n_cv_cells, control$seed)
    lam <- .trop_select_lambda(Y, W, grids, control, cv_cells, verbose)
  } else {
    lam <- utils::modifyList(list(time = 0, unit = 0, nn = Inf), lambda)
    lam$Q <- NA_real_
  }

  est <- .trop_att(Y, W, lam, control, anchor, pat)

  if (identical(se, "auto")) {
    se <- if (pat$n_treated_units >= 2) "jackknife" else "placebo"
  }
  inf <- .trop_se(Y, W, lam, control, anchor, pat, est, se, conf_level)

  list(
    estimate = est$att,
    std.error = inf$se,
    conf.low = inf$conf.low,
    conf.high = inf$conf.high,
    se.method = se,
    lambda = lam,
    anchor = anchor,
    tau_cells = est$tau_cells,
    counterfactual = est$Mhat,
    rank = est$rank
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
.trop_att <- function(Y, W, lam, ctrl, anchor, pat) {
  N <- nrow(Y); Tt <- ncol(Y)
  treated_cells <- which(W == 1, arr.ind = TRUE)
  Mhat <- matrix(NA_real_, N, Tt, dimnames = dimnames(Y))
  rnk <- NA_integer_

  if (anchor == "pooled") {
    du <- .unit_distance_pooled(Y, W, pat$treated_units)
    t_anchor <- sort(unique(treated_cells[, 2]))
    wmat <- .trop_weight_matrix(du, t_anchor, Tt, lam)
    fit <- .trop_solve(Y, W, wmat, lam, ctrl)
    Mhat <- fit$M
    rnk <- fit$rank
  } else {
    ranks <- integer(nrow(treated_cells))
    L_prev <- NULL                      # warm start: reuse the previous cell's L
    for (k in seq_len(nrow(treated_cells))) {
      i <- treated_cells[k, 1]; t <- treated_cells[k, 2]
      du <- .unit_distance_to(Y, W, i, t)
      wmat <- .trop_weight_matrix(du, t, Tt, lam)
      fit <- .trop_solve(Y, W, wmat, lam, ctrl, L_init = L_prev)
      Mhat[i, t] <- fit$M[i, t]
      ranks[k] <- fit$rank
      L_prev <- fit$L                   # consecutive cells share dimensions
    }
    rnk <- round(mean(ranks))
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
  list(att = mean(tau), tau_cells = tau_cells, Mhat = Mhat, rank = rnk)
}

#' Standard errors for the TROP ATT (jackknife / placebo).
#' @keywords internal
#' @noRd
.trop_se <- function(Y, W, lam, ctrl, anchor, pat, est, method, conf_level) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  na_out <- list(se = NA_real_, conf.low = NA_real_, conf.high = NA_real_)
  if (method == "none") return(na_out)

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
      .trop_att(Yk, Wk, lam, ctrl, anchor, patk)$att
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
      tryCatch(.trop_att(Yb, Wb, lam, ctrl, anchor, patb)$att,
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
  } else { # placebo
    controls <- setdiff(seq_len(nrow(Y)), pat$treated_units)
    if (length(controls) < 2) return(na_out)
    treat_pattern <- W[pat$treated_units, , drop = FALSE]
    # replicate the treated cells' time pattern onto each control unit
    tcols <- which(colSums(treat_pattern) > 0)
    par <- (ctrl$workers %||% 1L) > 1L
    placebo <- unlist(.par_lapply(seq_along(controls), function(k) {
      ci <- controls[k]
      Wp <- matrix(0, nrow(Y), ncol(Y), dimnames = dimnames(Y))
      Wp[ci, tcols] <- 1
      patp <- .assignment_pattern(Wp)
      .trop_att(Y, Wp, lam, ctrl, anchor, patp)$att
    }, parallel = par), use.names = FALSE)
    se <- stats::sd(placebo)
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
#'   loops (cross-validation cells, and the bootstrap / jackknife / placebo
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
  invisible(x)
}
