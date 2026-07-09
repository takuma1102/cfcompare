# Triply RObust Panel (TROP) estimator.
# Athey, Imbens, Qu & Viviano (2026), Journal of Applied Econometrics, doi:10.1002/jae.70061.
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

#' Unit distances anchored to the AVERAGE treated trajectory (pooled weighting)
#'
#' RMS distance of each unit's outcome path to the mean path of the treated
#' units, over the non-target periods (columns in which no treated unit is
#' active). This is the pooled unit-distance convention of the official Python
#' package (`TROP_TWFE_average`) and of the Stata command's pooled mode
#' (`trop_unit_weights2(..., pooled = 1)`); on a complete block design the three
#' coincide exactly. Cells where the unit itself is treated (possible under
#' general assignment patterns) and non-finite cells are additionally excluded,
#' which is a no-op on block designs.
#' @keywords internal
#' @noRd
.unit_distance_to_avg <- function(Y, W, treated_units) {
  N <- nrow(Y); Tt <- ncol(Y)
  tcols <- which(colSums(W[treated_units, , drop = FALSE]) > 0)
  avg <- colMeans(Y[treated_units, , drop = FALSE], na.rm = TRUE)
  keep <- rep(TRUE, Tt)
  keep[tcols] <- FALSE                          # drop the target periods
  useM <- matrix(keep, N, Tt, byrow = TRUE) & (W == 0) & is.finite(Y) &
    matrix(is.finite(avg), N, Tt, byrow = TRUE)
  diff2 <- (Y - matrix(avg, N, Tt, byrow = TRUE))^2
  diff2[!useM] <- 0
  n <- rowSums(useM); ssq <- rowSums(diff2)
  d <- rep(NA_real_, N); pos <- n > 0
  d[pos] <- sqrt(ssq[pos] / n[pos])
  d
}

#' Centre (1-based, possibly half-integer) of the treated block of periods.
#'
#' Generalises the trailing-block centre `T - tp/2` (0-based; the official
#' Python convention) to treated blocks anywhere in the panel, following the
#' Stata command: centre = (min + max + 1) / 2 over the 0-based treated
#' columns, i.e. (min + max + 1) / 2 in 1-based indexing. For a trailing block
#' of width `tp` this equals the previous `(T - tp/2) + 1`, so results on
#' block designs are unchanged.
#' @keywords internal
#' @noRd
.treated_block_center <- function(treated_cols) {
  (min(treated_cols) + max(treated_cols) + 1) / 2
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
.trop_solve <- function(Y, W, wmat, lam, ctrl, state_init = NULL, X = NULL) {
  mask <- (W == 0) * 1
  .mcnnm_fit(Y, mask, wmat, lam$nn,
             max_iter = ctrl$max_iter, tol = ctrl$tol,
             state_init = state_init, svd_method = ctrl$svd %||% "truncated",
             X = X)
}

#' Extract the warm-startable state from a solver fit.
#' @keywords internal
#' @noRd
.solver_state <- function(fit) fit[c("L", "alpha", "beta", "grand", "phi")]

# ---- leave-one-out cross-validation (eq. 4-5) ------------------------------

#' LOO-CV criterion Q(lambda) on a subsample of control cells.
#'
#' For each sampled control cell (i, t) we treat it as if it were treated
#' (mask it out of the loss), build weights anchored to that cell, solve eq. (2),
#' and accumulate the squared prediction error -- a subsampled realisation of
#' eq. (5).
#'
#' Returns `list(Q, warm)`. `warm` is a per-cell list of solver states
#' (fixed effects, phi and low-rank part; see `.solver_state()`) when
#' `keep_state = TRUE` (else `NULL`); passing it back via `warm` on the next
#' call warm-starts each cell's solve from the previous penalty's solution.
#' Used by the lambda_nn grid sweep in `.trop_select_lambda()` (the
#' "warm-start nuclear path" of the official Stata command): the problem is
#' convex, so warm starts change only the iteration count, not the solution.
#' @keywords internal
#' @noRd
.trop_cv_Q <- function(Y, W, lam, ctrl, cv_cells, du_list = NULL, X = NULL,
                       warm = NULL, keep_state = FALSE) {
  Tt <- ncol(Y)
  one <- function(k) {
    i <- cv_cells[k, 1]; t <- cv_cells[k, 2]
    Wk <- W
    Wk[i, t] <- 1                       # hold this control cell out
    # unit distances do not depend on the penalties: reuse a cached value when
    # one is supplied (see .trop_select_lambda), else compute it here.
    du <- if (is.null(du_list)) .unit_distance_to(Y, Wk, i, t) else du_list[[k]]
    wmat <- .trop_weight_matrix(du, t, Tt, lam)
    fit <- .trop_solve(Y, Wk, wmat, lam, ctrl,
                       state_init = if (is.null(warm)) NULL else warm[[k]],
                       X = X)
    err <- (Y[i, t] - fit$M[i, t])^2
    if (keep_state) list(err = err, state = .solver_state(fit)) else err
  }
  par <- (ctrl$workers %||% 1L) > 1L
  res <- .par_lapply(seq_len(nrow(cv_cells)), one, parallel = par)
  if (keep_state) {
    list(Q = mean(vapply(res, function(r) r$err, numeric(1))),
         warm = lapply(res, function(r) r$state))
  } else {
    list(Q = mean(unlist(res, use.names = FALSE)), warm = NULL)
  }
}

# ---- placebo-RMSE cross-validation (option; Python/Stata-style) -------------

#' Placebo draws for placebo-RMSE penalty selection.
#'
#' Restricts the panel to control units (drops the real treated rows), then
#' draws `ctrl$n_placebo` placebo assignments. Each draw stamps the *actual*
#' treated-unit adoption patterns (the rows of `W` for the real treated units)
#' onto randomly chosen control units, one pattern per sampled unit -- the
#' pattern-bank scheme of the official Stata command (v0.2.4,
#' `trop_placebo_rmse_path()` / `__trop_Pat`). Under block adoption every
#' pattern is the same trailing block, so this reduces exactly to the previous
#' common-post-block behaviour; under staggered adoption each placebo draw
#' reproduces the real cohort structure (who adopts when).
#'
#' Placebo units are grouped into adoption cohorts (identical W rows). For each
#' cohort the draw caches its units, treated columns, block centre and unit
#' distances (none of which depend on the penalties), so the coordinate-descent
#' search reuses them across every penalty evaluation.
#' @keywords internal
#' @noRd
.trop_placebo_setup <- function(Y, W, pat, ctrl) {
  ctrl_units <- setdiff(seq_len(nrow(Y)), pat$treated_units)
  Yp <- Y[ctrl_units, , drop = FALSE]
  # Pattern bank: one adoption pattern per real treated unit.
  Pat <- W[pat$treated_units, , drop = FALSE]
  post_cols <- which(colSums(Pat) > 0)
  n_tr <- nrow(Pat)
  n_pl <- ctrl$n_placebo %||% 30L
  if (length(ctrl_units) <= n_tr || length(post_cols) == 0L ||
      length(post_cols) >= ncol(Yp)) {
    stop("Not enough control units / periods for placebo cross-validation.",
         call. = FALSE)
  }
  if (!is.null(ctrl$seed)) {
    old <- .Random.seed_safe()
    on.exit(.Random.seed_restore(old), add = TRUE)
    set.seed(ctrl$seed)
  }
  # Cohorts of the pattern bank: rows of Pat that are identical adopt together.
  pat_key <- apply(Pat, 1L, paste, collapse = "")
  cohorts <- split(seq_len(n_tr), pat_key)
  draws <- lapply(seq_len(n_pl), function(b) {
    tr <- sample.int(nrow(Yp), n_tr)
    Wp <- matrix(0, nrow(Yp), ncol(Yp))
    Wp[tr, ] <- Pat                       # stamp one real pattern per unit
    groups <- lapply(cohorts, function(idx) {
      us   <- tr[idx]
      cols <- which(Pat[idx[1L], ] == 1)
      list(us = us, cols = cols,
           size     = length(us) * length(cols),
           t_center = .treated_block_center(cols),
           du       = .unit_distance_to_avg(Yp, Wp, us))
    })
    names(groups) <- NULL
    list(Wp = Wp, groups = groups, post_cols = post_cols)
  })
  list(Yp = Yp, ctrl_units = ctrl_units, draws = draws)
}

#' Placebo-RMSE criterion: mean squared placebo ATT across the draws.
#'
#' Each placebo ATT is a pooled estimate on the control panel: within each
#' adoption cohort of the draw the weights are anchored to that cohort's
#' treated units and block centre, the cohort placebo ATT is the mean gap over
#' its placebo cells, and cohorts are combined cell-weighted -- mirroring the
#' per-spell, cell-weighted scoring of the official Stata command
#' (`trop_placebo_rmse_path()`), which "rehearses" the pooled estimand under
#' the real adoption pattern. With a single cohort (block adoption) this
#' collapses to the previous single-anchor pooled ATT, matching the
#' `TROP_TWFE_average` placebo tuning of the official Python package. The true
#' placebo effect is zero, so smaller is better. Penalty selection therefore
#' uses the pooled ATT even when the final fit uses `anchor = "per_cell"`.
#'
#' Returns `list(Q, warm)` like `.trop_cv_Q()`; `warm` holds, per placebo
#' draw, a list of fitted solver states (one per cohort) when
#' `keep_state = TRUE`, enabling the warm-start nuclear path in
#' `.trop_select_lambda()`.
#' @keywords internal
#' @noRd
.trop_cv_Q_placebo <- function(pb, lam, ctrl, X = NULL,
                               warm = NULL, keep_state = FALSE) {
  Yp <- pb$Yp; Tt <- ncol(Yp)
  one <- function(k) {
    d <- pb$draws[[k]]
    ng <- length(d$groups)
    states <- if (keep_state) vector("list", ng) else NULL
    att <- 0; wsum <- 0
    for (g in seq_len(ng)) {
      gr <- d$groups[[g]]
      wmat <- .trop_weight_matrix(gr$du, gr$t_center, Tt, lam)
      fit  <- .trop_solve(Yp, d$Wp, wmat, lam, ctrl,
                          state_init = if (is.null(warm)) NULL
                                       else warm[[k]][[g]],
                          X = X)
      cells <- matrix(FALSE, nrow(Yp), Tt)
      cells[gr$us, gr$cols] <- TRUE       # this cohort's placebo cells
      tau_g <- mean((Yp - fit$M)[cells])
      att   <- att + gr$size * tau_g
      wsum  <- wsum + gr$size
      if (keep_state) states[[g]] <- .solver_state(fit)
    }
    err <- (att / wsum)^2
    if (keep_state) list(err = err, state = states) else err
  }
  par <- (ctrl$workers %||% 1L) > 1L
  res <- .par_lapply(seq_along(pb$draws), one, parallel = par)
  if (keep_state) {
    list(Q = mean(vapply(res, function(r) r$err, numeric(1))),
         warm = lapply(res, function(r) r$state))
  } else {
    list(Q = mean(unlist(res, use.names = FALSE)), warm = NULL)
  }
}

#' Coordinate-descent search over (lambda_time, lambda_unit, lambda_nn).
#'
#' Follows the warm-start scheme of footnote 2 of the paper: start from
#' lambda_nn = Inf, lambda_unit = 0, optimise each penalty marginally, then
#' cycle.
#' @keywords internal
#' @noRd
.trop_select_lambda <- function(Y, W, grids, ctrl, cv_cells, verbose = FALSE,
                                X = NULL, pat = NULL) {
  lam <- list(time = 0, unit = 0, nn = Inf)
  if (identical(ctrl$cv_method %||% "loocv", "placebo")) {
    # Placebo-RMSE criterion (option): assign placebo blocks to control units and
    # minimise the mean squared placebo ATT (as in the official Python/Stata
    # packages). See .trop_placebo_setup() / .trop_cv_Q_placebo().
    if (is.null(pat)) pat <- .assignment_pattern(W)
    pb <- .trop_placebo_setup(Y, W, pat, ctrl)
    Xp <- if (is.null(X)) NULL else
      lapply(.as_cov_list(X, nrow(Y), ncol(Y)),
             function(M) M[pb$ctrl_units, , drop = FALSE])
    eval_one <- function(lam, warm = NULL, keep_state = FALSE)
      .trop_cv_Q_placebo(pb, lam, ctrl, X = Xp, warm = warm, keep_state = keep_state)
  } else {
    # Default: leave-one-control-cell-out prediction error (eq. (4)-(5)). Unit
    # distances for each held-out CV cell do not depend on the penalties, so
    # compute them once and reuse across the whole penalty search (big saving:
    # the search evaluates the CV criterion dozens of times).
    #
    # Under staggered / general adoption, holding out one control cell at a
    # time cannot mimic the real missingness pattern (whole post-periods
    # missing per adoption cohort), which is the known failure mode of LOOCV
    # in those designs (cf. the official Stata command, which restricts
    # cv(loocv) to per-cell mode and defaults to pattern-resampled placebo CV
    # for the pooled estimand). LOOCV remains the default here, but flag it.
    if (is.null(pat)) pat <- .assignment_pattern(W)
    if (!identical(pat$type, "block")) {
      warning("Treatment adoption is staggered / non-block; leave-one-cell-out ",
              "CV cannot mimic the design's missingness pattern and may select ",
              "unreliable penalties. Consider trop_control(cv_method = ",
              "\"placebo\"), which resamples the actual adoption patterns.",
              call. = FALSE)
    }
    du_cache <- lapply(seq_len(nrow(cv_cells)), function(k) {
      i <- cv_cells[k, 1]; t <- cv_cells[k, 2]
      Wk <- W; Wk[i, t] <- 1
      .unit_distance_to(Y, Wk, i, t)
    })
    eval_one <- function(lam, warm = NULL, keep_state = FALSE)
      .trop_cv_Q(Y, W, lam, ctrl, cv_cells, du_list = du_cache, X = X,
                 warm = warm, keep_state = keep_state)
  }

  pick <- function(lam, which, grid) {
    # Warm-start nuclear path (as in the official Stata command): sweep the
    # lambda_nn grid from the strongest penalty (Inf, where L = 0) down to the
    # weakest, initialising each held-out solve with the full solver state
    # (fixed effects, phi and low-rank part) fitted at the previous, larger
    # penalty. The problem in eq. (2) is convex, so this changes only the
    # iteration count, never the selected value. The time/unit sweeps solve
    # from a cold start as before.
    is_nn <- identical(which, "nn")
    if (is_nn) grid <- sort(grid, decreasing = TRUE)   # Inf first
    warm <- NULL
    qs <- numeric(length(grid))
    for (gi in seq_along(grid)) {
      l2 <- lam; l2[[which]] <- grid[gi]
      ev <- eval_one(l2, warm = if (is_nn) warm else NULL, keep_state = is_nn)
      qs[gi] <- ev$Q
      if (is_nn) warm <- ev$warm
    }
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

  # cycles; stop early once the selected triplet stops changing (the official
  # Python coordinate-descent iterates to a fixed point the same way).
  for (cyc in seq_len(ctrl$cv_cycles)) {
    prev <- lam[c("time", "unit", "nn")]
    lam$time <- pick(lam, "time", grids$time)
    lam$unit <- pick(lam, "unit", grids$unit)
    lam$nn   <- pick(lam, "nn",   grids$nn)
    if (identical(prev, lam[c("time", "unit", "nn")])) break
  }
  lam$Q <- eval_one(lam)$Q
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
  # The solver's soft-threshold is lambda_nn / (2 * max(weight)) (paper eq. (2)
  # loss, with no 1/2 factor), so with uniform weights the low-rank term is fully
  # shrunk to zero once lambda_nn >= 2 * s1. These multipliers spread the finite
  # grid over the non-degenerate range up to that point (the largest value nearly
  # kills L; Inf gives exactly L = 0), and keep lambda_nn on the official
  # Python/Stata scale. (Halving them would target the 1/2-loss variant and
  # under-cover the strongly-regularised region.)
  nn_grid <- c(Inf, s1 * c(1.0, 0.5, 0.2, 0.1, 0.04))

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
#' Fits the TROP estimator of Athey, Imbens, Qu & Viviano (2026) on a long
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
#'   penalties and skip cross-validation. When `standardize = TRUE`, supplied
#'   penalties are interpreted on the standardized-outcome scale.
#' @param standardize Logical; if `TRUE`, the outcome is standardized internally
#'   as `(Y - mean(Y)) / sd(Y)` (over all observed cells) before fitting, and
#'   the ATT, standard error, confidence interval, per-cell effects and
#'   counterfactual are mapped back to the raw outcome scale. This matches the
#'   convention of the official Stata command (v0.2.1+); as there, the selected
#'   `lambda` values are reported on the *standardized* scale, which makes them
#'   comparable across outcomes and across implementations. Default `FALSE`
#'   (raw outcome, the previous behaviour). The returned `panel$Y` is on the
#'   fitting (standardized) scale, with the centre/scale stored in `$scaling`.
#' @param anchor How weights are anchored to treated cells: `"per_cell"`
#'   re-solves eq. (2) with cell-specific weights for every treated cell
#'   (faithful to the paper); `"pooled"` solves once with weights anchored to the
#'   treated set (fast); `"auto"` (default) uses `per_cell` when there are at
#'   most `max_cells` treated cells and `pooled` otherwise. `trop_matrix()` has 
#'   no `anchor` argument; it is always pooled/block-center.
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
#'   counterfactual matrix, weights, the reshaped panel, and `scaling` (the
#'   outcome centre/scale used internally; the identity unless
#'   `standardize = TRUE`). All effect and outcome-level components are on the
#'   raw outcome scale; `panel$Y` and `lambda` are on the fitting scale.
#' @references Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026).
#'   Triply Robust Panel Estimators. \emph{Journal of Applied Econometrics}. \doi{10.1002/jae.70061}.
#' @examples
#' df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, seed = 1)
#' fit <- trop(df, "y", "w", "id", "t", se = "none",
#'             control = trop_control(n_cv_cells = 8L, cv_cycles = 1L))
#' fit
#' @export
trop <- function(data, outcome, treatment, unit, time,
                 covariates = NULL,
                 lambda = NULL,
                 standardize = FALSE,
                 anchor = c("auto", "per_cell", "pooled"),
                 se = c("bootstrap", "auto", "jackknife", "none"),
                 grids = NULL,
                 control = trop_control(),
                 verbose = FALSE) {
  anchor <- match.arg(anchor)
  se <- match.arg(se)
  m <- .panel_to_matrices(data, outcome, treatment, unit, time)
  # Optional outcome standardization (Stata convention): fit on
  # (Y - mean) / sd and map the results back to the raw scale below.
  sc <- .trop_scaling(m$Y, standardize)
  m$Y <- (m$Y - sc$center) / sc$scale
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
  eng <- .trop_unscale_engine(eng, sc)

  structure(
    c(eng, list(
      conf.level = control$conf_level,
      pattern = pat,
      panel = m,
      scaling = sc,
      outcome = outcome,
      covariates = covariates,
      call = match.call()
    )),
    class = "trop"
  )
}

#' Outcome scaling used by [trop()] (`standardize` option).
#'
#' Grand mean / SD over all observed (finite) outcome cells, as in the official
#' Stata command's `trop_standardize_outcome()`. A degenerate SD (zero, or not
#' finite) falls back to the identity transform.
#' @keywords internal
#' @noRd
.trop_scaling <- function(Y, standardize) {
  if (!isTRUE(standardize)) return(list(center = 0, scale = 1))
  v <- Y[is.finite(Y)]
  ctr <- mean(v)
  s <- stats::sd(v)
  if (!is.finite(ctr)) ctr <- 0
  if (!is.finite(s) || s <= 0) s <- 1
  list(center = ctr, scale = s)
}

#' Map a standardized-scale engine result back to the raw outcome scale.
#'
#' Effects (ATT, SE, CI half-widths, per-cell tau, covariate coefficients) scale
#' by `s`; levels (observed / counterfactual outcomes) map through
#' `center + s * y`. The selected penalties are left on the fitting
#' (standardized) scale, as in the Stata command.
#' @keywords internal
#' @noRd
.trop_unscale_engine <- function(eng, sc) {
  s <- sc$scale; ctr <- sc$center
  if (identical(s, 1) && identical(ctr, 0)) return(eng)
  eng$estimate  <- s * eng$estimate
  eng$std.error <- s * eng$std.error
  eng$conf.low  <- s * eng$conf.low
  eng$conf.high <- s * eng$conf.high
  if (!is.null(eng$tau_cells)) {
    eng$tau_cells$y      <- ctr + s * eng$tau_cells$y
    eng$tau_cells$y0_hat <- ctr + s * eng$tau_cells$y0_hat
    eng$tau_cells$tau    <- s * eng$tau_cells$tau
  }
  if (!is.null(eng$counterfactual))
    eng$counterfactual <- ctr + s * eng$counterfactual
  if (length(eng$phi)) eng$phi <- s * eng$phi
  eng
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
    lam <- .trop_select_lambda(Y, W, grids, control, cv_cells, verbose, X = X,
                               pat = pat)
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
  # n = Inf or n <= 0 means "all control cells": the paper's full eq. (5)
  # criterion, matching the Stata command's cells(0).
  if (!is.finite(n) || n <= 0) return(ctrl_idx)
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
    # Pooled weights: RMS distance to the AVERAGE treated trajectory, matching
    # the matrix-in trop_matrix()/.trop_reference_weights, the official Python
    # TROP_TWFE_average and the Stata pooled mode. (Averaging per-treated-unit
    # distances instead is not the reference convention.)
    du <- .unit_distance_to_avg(Y, W, pat$treated_units)
    # Anchor the time weights to the CENTRE of the treated block, matching the
    # references (for a trailing block, dist_time = |s - (T - tp/2)| 0-based).
    # .treated_block_center() generalises this to treated blocks anywhere in
    # the panel (Stata convention). Passing the vector of treated periods to
    # .trop_weight_matrix() would instead use the distance to the nearest
    # treated period, which flattens theta across the block and does not match
    # the reference.
    t_center <- .treated_block_center(unique(treated_cells[, 2]))
    wmat <- .trop_weight_matrix(du, t_center, Tt, lam)
    fit <- .trop_solve(Y, W, wmat, lam, ctrl, X = X)
    Mhat <- fit$M
    rnk <- fit$rank
    phi <- fit$phi
  } else {
    ranks <- integer(nrow(treated_cells))
    phis <- vector("list", nrow(treated_cells))
    st_prev <- NULL       # warm start: reuse the previous cell's full state
    for (k in seq_len(nrow(treated_cells))) {
      i <- treated_cells[k, 1]; t <- treated_cells[k, 2]
      du <- .unit_distance_to(Y, W, i, t)
      wmat <- .trop_weight_matrix(du, t, Tt, lam)
      fit <- .trop_solve(Y, W, wmat, lam, ctrl, state_init = st_prev, X = X)
      Mhat[i, t] <- fit$M[i, t]
      ranks[k] <- fit$rank
      phis[[k]] <- fit$phi
      st_prev <- .solver_state(fit)     # consecutive cells share dimensions
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
  rn <- rownames(Y) %||% as.character(seq_len(N))
  cn <- colnames(Y) %||% as.character(seq_len(Tt))
  tau_cells <- data.frame(
    unit = rn[treated_cells[, 1]],
    time = cn[treated_cells[, 2]],
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
#' @param max_iter Maximum solver iterations. Higher by default because small
#'   `lambda_nn` converges slowly under the soft-impute iteration.
#' @param tol Reference solver convergence tolerance. The applied tolerance is
#'   tightened automatically as `lambda_nn` weakens (small penalties converge
#'   slowly), keeping accuracy roughly constant across the penalty grid; a fixed
#'   tolerance otherwise under-converges at small `lambda_nn`. See [trop()].
#' @param n_cv_cells Number of control cells sampled for the CV criterion. Set
#'   to `Inf` (or `0`) to use *every* control cell, i.e. the paper's full
#'   eq. (5) criterion (the Stata command's `cells(0)`); slower but exact.
#' @param cv_cycles Number of coordinate-descent cycles in penalty selection.
#' @param cv_method Penalty cross-validation criterion. `"loocv"` (default)
#'   scores held-out control-cell prediction error (the paper's eq. (4)-(5));
#'   `"placebo"` instead stamps the *actual* treated-unit adoption patterns
#'   onto randomly resampled control units and minimises the mean squared
#'   placebo ATT (the placebo-RMSE / pattern-resampling criterion of the
#'   official Python and Stata packages). Under block adoption the two placebo
#'   variants coincide; under staggered adoption each placebo draw reproduces
#'   the real cohort structure, so `"placebo"` is the recommended choice there
#'   (a warning is emitted if `"loocv"` is used on a staggered design).
#'   Ignored when the penalties are supplied via `lambda` (CV is then
#'   skipped). The placebo criterion always uses the pooled (per-cohort
#'   block-centre, cell-weighted) ATT, even when the final fit uses
#'   `anchor = "per_cell"`.
#' @param n_placebo Number of placebo assignments drawn when
#'   `cv_method = "placebo"` (larger is more stable but slower).
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
trop_control <- function(max_iter = 2000L, tol = 1e-6,
                         n_cv_cells = 120L, cv_cycles = 2L,
                         max_cells = 60L, conf_level = 0.95,
                         n_boot = 200L, boot_ci = c("percentile", "normal"),
                         svd = c("truncated", "full"),
                         cv_method = c("loocv", "placebo"), n_placebo = 30L,
                         workers = 1L,
                         seed = NULL) {
  boot_ci <- match.arg(boot_ci)
  svd <- match.arg(svd)
  cv_method <- match.arg(cv_method)
  list(max_iter = max_iter, tol = tol, n_cv_cells = n_cv_cells,
       cv_cycles = cv_cycles, max_cells = max_cells,
       conf_level = conf_level, n_boot = n_boot, boot_ci = boot_ci,
       svd = svd, cv_method = cv_method, n_placebo = n_placebo,
       workers = workers,
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
  std <- !is.null(x$scaling) && (x$scaling$scale != 1 || x$scaling$center != 0)
  cat(sprintf("  penalties%s: lambda_time=%.3g, lambda_unit=%.3g, lambda_nn=%s\n",
              if (std) " (standardized-outcome scale)" else "",
              x$lambda$time, x$lambda$unit,
              ifelse(is.infinite(x$lambda$nn), "Inf",
                     sprintf("%.3g", x$lambda$nn))))
  if (std) {
    cat(sprintf("  outcome standardized internally (mean %.4g, sd %.4g); ATT/SE/CI on the raw scale\n",
                x$scaling$center, x$scaling$scale))
  }
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
  du <- .unit_distance_to_avg(Y, W, tu)
  # Block-centre time anchoring, consistent with .trop_att(anchor = "pooled"),
  # trop_matrix() and the official Python TROP_TWFE_average.
  tcells   <- which(W == 1, arr.ind = TRUE)
  t_center <- .treated_block_center(unique(tcells[, 2]))
  wmat <- .trop_weight_matrix(du, t_center, ncol(Y), lam)
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
#' @references Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026).
#'   Triply Robust Panel Estimators. \emph{Journal of Applied Econometrics}. \doi{10.1002/jae.70061}.
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
  # panel$Y (and lam) live on the fitting scale; effects computed below are
  # mapped back to the raw outcome scale via the fit's stored scaling.
  s_out <- (object$scaling %||% list(scale = 1))$scale
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
      estimate   = s_out * as.numeric(point),
      std.error  = s_out * as.numeric(inf$se[as.character(evt)]),
      conf.low   = s_out * as.numeric(inf$conf.low[as.character(evt)]),
      conf.high  = s_out * as.numeric(inf$conf.high[as.character(evt)]),
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
