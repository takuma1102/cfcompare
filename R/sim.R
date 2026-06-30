# Synthetic panel generator for examples, tests and quick experiments.

#' Simulate a panel from a low-rank factor model
#'
#' Generates a long panel whose untreated potential outcomes follow an
#' interactive-fixed-effects (factor) model on top of two-way fixed effects,
#' in the spirit of the data-generating processes in Athey, Imbens, Qu &
#' Viviano (2025). A block treatment is applied to `n_treated` units from period
#' `t0` onward, with a constant additive effect `att`.
#'
#' @param N Number of units.
#' @param T Number of periods.
#' @param n_treated Number of treated units.
#' @param t0 First treated period (block design).
#' @param rank Number of latent factors.
#' @param att True treatment effect added to treated cells.
#' @param noise Standard deviation of idiosyncratic noise.
#' @param seed Optional RNG seed.
#' @param n_cov Number of time-varying covariates to generate (default `0`,
#'   none). Generated covariates are cell-level standard-normal draws whose
#'   linear index is part of the untreated potential outcome `Y(0)`, matching the
#'   additive covariate model of [trop()] (`L = X.phi + R`); pass their names to
#'   `covariates` in [trop()] / [panel_compare()] to recover the coefficients.
#'   Defaults to `length(phi)` when `phi` is supplied.
#' @param phi Optional numeric coefficient vector for the covariates (length
#'   `n_cov`). When `NULL` and `n_cov > 0`, coefficients are drawn from a
#'   standard normal. The coefficients actually used are returned on the
#'   `"phi"` attribute of the result.
#' @param confounding Non-negative strength of selection-on-unobservables in the
#'   treatment assignment. With `confounding = 0` (default) the first
#'   `n_treated` units are treated independently of the latent factors, so
#'   treatment is ignorable and DID/TWFE is unbiased. With `confounding > 0`
#'   treated units are selected on their latent factor loadings (units loading
#'   heavily on the post-period factor direction are more likely to be treated),
#'   which biases DID/TWFE while factor-aware (MC, TROP) and weighting (SDID)
#'   estimators stay consistent. Larger values give stronger selection.
#' @return A long `data.frame` with columns `id`, `t`, `y`, `w`, and the noiseless
#'   counterfactual `y0` (useful for evaluating estimators). When `n_cov > 0` it
#'   additionally carries covariate columns `x1`, ..., `x{n_cov}`, and the true
#'   coefficients are stored on `attr(., "phi")`.
#' @examples
#' df <- sim_panel(N = 30, T = 15, n_treated = 5, t0 = 11, seed = 42)
#' head(df)
#'
#' # Two covariates with known coefficients, recoverable by trop(covariates=):
#' dc <- sim_panel(N = 30, T = 12, n_treated = 4, t0 = 9,
#'                 n_cov = 2, phi = c(1.2, -0.7), seed = 1)
#' attr(dc, "phi")
#' @export
sim_panel <- function(N = 30, T = 20, n_treated = 5, t0 = NULL,
                      rank = 3L, att = 1, noise = 0.5, seed = NULL,
                      n_cov = 0L, phi = NULL, confounding = 0) {
  if (!is.null(seed)) {
    old <- .Random.seed_safe(); on.exit(.Random.seed_restore(old), add = TRUE)
    set.seed(seed)
  }
  if (is.null(t0)) t0 <- floor(0.75 * T) + 1L
  stopifnot(n_treated < N, t0 >= 2, t0 <= T)

  # resolve covariate count / coefficients
  if (!is.null(phi)) {
    phi <- as.numeric(phi)
    if (missing(n_cov) || n_cov == 0L) n_cov <- length(phi)
    if (length(phi) != n_cov)
      stop("`phi` must have length `n_cov` (", n_cov, ").", call. = FALSE)
  }
  n_cov <- as.integer(n_cov)
  if (n_cov < 0L) stop("`n_cov` must be non-negative.", call. = FALSE)
  if (n_cov > 0L && is.null(phi)) phi <- stats::rnorm(n_cov)

  # two-way fixed effects
  alpha <- stats::rnorm(N, sd = 1)
  beta <- cumsum(stats::rnorm(T, sd = 0.3))          # smooth common trend

  # low-rank interactive component L = F %*% t(G)
  Fmat <- matrix(stats::rnorm(N * rank), N, rank)
  Gmat <- matrix(stats::rnorm(T * rank), T, rank)
  # give factors some temporal smoothness
  for (r in seq_len(rank)) Gmat[, r] <- cumsum(Gmat[, r]) / sqrt(T)
  L <- Fmat %*% t(Gmat)

  Y0 <- outer(alpha, rep(1, T)) + outer(rep(1, N), beta) + L
  noise_mat <- matrix(stats::rnorm(N * T, sd = noise), N, T)
  Y0_obs <- Y0 + noise_mat

  # covariates: cell-level regressors whose linear index X.phi is part of Y(0)
  # (both the observed outcome and the counterfactual y0 carry it), so trop()
  # with `covariates=` recovers `phi`.
  Xlist <- NULL
  if (n_cov > 0L) {
    Xlist <- lapply(seq_len(n_cov), function(k) matrix(stats::rnorm(N * T), N, T))
    Xsig <- Reduce(`+`, Map(function(Xk, b) b * Xk, Xlist, phi))
    Y0_obs <- Y0_obs + Xsig
  }

  # Treatment assignment. By default (confounding = 0) the first `n_treated`
  # units are treated, independent of the latent factors, so treatment is
  # ignorable and DID/TWFE is unbiased. With confounding > 0 the treated units
  # are instead selected on their latent factor loadings (paper-style selection
  # on unobservables): units whose loadings load heavily on the post-period
  # factor direction are more likely to be treated. This biases DID/TWFE (which
  # cannot model the interactive term) while factor-aware estimators (MC, TROP)
  # and weighting estimators (SDID) stay consistent. Larger values = stronger
  # selection.
  if (confounding > 0) {
    g_post <- colMeans(Gmat[t0:T, , drop = FALSE])        # post-period factor direction
    load   <- as.numeric(Fmat %*% g_post)
    load   <- (load - mean(load)) / stats::sd(load)
    # `confounding` sets the signal-to-noise ratio of selection: larger values
    # pick treated units more deterministically on their loadings (stronger
    # selection); as confounding -> 0 the score is dominated by the noise term,
    # i.e. selection becomes effectively random (ignorable), continuous with the
    # confounding = 0 branch below. (Note: the multiplier must scale signal
    # relative to a fixed-variance noise -- scaling the whole score would be a
    # no-op because order() is invariant to a positive multiplier.)
    score  <- confounding * load + stats::rnorm(N)
    treated_units <- sort(order(score, decreasing = TRUE)[seq_len(n_treated)])
  } else {
    treated_units <- seq_len(n_treated)
  }
  W <- matrix(0, N, T)
  W[treated_units, t0:T] <- 1
  Y <- Y0_obs + att * W

  df <- data.frame(
    id = rep(seq_len(N), times = T),
    t = rep(seq_len(T), each = N),
    y = as.numeric(Y),
    w = as.numeric(W),
    y0 = as.numeric(Y0_obs),
    stringsAsFactors = FALSE
  )
  if (n_cov > 0L)
    for (k in seq_len(n_cov)) df[[paste0("x", k)]] <- as.numeric(Xlist[[k]])
  df <- df[order(df$id, df$t), ]
  if (n_cov > 0L)
    attr(df, "phi") <- stats::setNames(phi, paste0("x", seq_len(n_cov)))
  df
}

#' Build a semi-synthetic panel from real data
#'
#' Takes a real long panel, uses its outcomes (optionally smoothed through a
#' low-rank-plus-two-way-fixed-effects fit) as the untreated potential outcomes
#' `Y(0)`, and imposes a *known* treatment effect on a chosen block of units and
#' periods. Because the baseline is real data but the effect is known, the result
#' is a ground-truth benchmark that "closely matches" the real setting -- the
#' style of semi-synthetic experiment used to evaluate panel estimators in
#' Athey, Imbens, Qu & Viviano (2026). Pair it with [panel_compare()] or
#' [panel_rmse()] to score estimators against the truth.
#'
#' @param data A real long `data.frame`, one row per unit-time.
#' @param outcome,unit,time Column names (strings).
#' @param n_treated Number of units to assign to the placebo treated group
#'   (sampled at random from all units).
#' @param t0 First treated period (block design). Defaults to about three
#'   quarters of the way through the panel.
#' @param att Constant additive treatment effect imposed on treated cells. Ignored
#'   if `effect` is supplied.
#' @param effect Optional per-treated-period effect: a single number, or a numeric
#'   vector of length `T - t0 + 1` giving a dynamic effect path.
#' @param baseline `"observed"` uses the real outcomes directly as `Y(0)`;
#'   `"lowrank"` replaces them with a low-rank + two-way-FE fit (optionally plus
#'   resampled residuals scaled by `noise`) for a smoother synthetic baseline.
#' @param lambda_nn Nuclear-norm penalty for the `"lowrank"` baseline fit.
#' @param noise For `"lowrank"`, standard-deviation multiplier on resampled
#'   residuals added back to the fit (0 = noiseless baseline).
#' @param seed Optional RNG seed.
#' @return A long `data.frame` with columns `id`, `t`, `y`, `w`, `y0` (the imposed
#'   untreated potential outcome) and `tau` (the true effect, 0 off treatment).
#' @seealso [sim_panel()], [panel_compare()], [panel_rmse()]
#' @export
#' @examples
#' real <- sim_panel(N = 40, T = 18, n_treated = 0L, att = 0, seed = 1)
#' ss <- sim_semisynthetic(real, "y", "id", "t",
#'                         n_treated = 6, t0 = 14, att = 3, seed = 2)
#' mean(ss$tau[ss$w == 1])   # true ATT = 3
sim_semisynthetic <- function(data, outcome, unit, time,
                              n_treated, t0 = NULL, att = 1, effect = NULL,
                              baseline = c("observed", "lowrank"),
                              lambda_nn = NULL, noise = 0, seed = NULL) {
  baseline <- match.arg(baseline)
  if (!is.null(seed)) {
    old <- .Random.seed_safe(); on.exit(.Random.seed_restore(old), add = TRUE)
    set.seed(seed)
  }
  units <- sort(unique(data[[unit]]))
  times <- sort(unique(data[[time]]))
  N <- length(units); Tt <- length(times)
  if (n_treated < 1 || n_treated >= N)
    stop("n_treated must be in 1..N-1.", call. = FALSE)
  if (is.null(t0)) t0 <- times[floor(0.75 * Tt) + 1L]
  ti0 <- match(t0, times)
  if (is.na(ti0) || ti0 < 2 || ti0 > Tt) stop("t0 not a valid period.", call. = FALSE)

  # outcome matrix (units x times)
  Y <- matrix(NA_real_, N, Tt, dimnames = list(as.character(units),
                                               as.character(times)))
  ui <- match(data[[unit]], units); tj <- match(data[[time]], times)
  Y[cbind(ui, tj)] <- data[[outcome]]
  if (anyNA(Y)) stop("Panel is unbalanced; sim_semisynthetic needs a full grid.",
                     call. = FALSE)

  # untreated potential outcome
  if (baseline == "lowrank") {
    mask <- matrix(1, N, Tt)
    lam_nn <- lambda_nn %||% (0.1 * mean(abs(Y)))
    fit <- .mcnnm_fit(Y, mask, matrix(1, N, Tt), lam_nn)
    Y0 <- fit$M
    if (noise > 0) {
      resid <- Y - fit$M
      Y0 <- Y0 + noise * matrix(sample(resid, N * Tt, replace = TRUE), N, Tt)
    }
  } else {
    Y0 <- Y
  }

  # block treatment on randomly chosen units
  treated_units <- sort(sample.int(N, n_treated))
  W <- matrix(0, N, Tt)
  W[treated_units, ti0:Tt] <- 1

  # effect path
  npost <- Tt - ti0 + 1L
  if (is.null(effect)) {
    eff_path <- rep(att, npost)
  } else if (length(effect) == 1L) {
    eff_path <- rep(effect, npost)
  } else if (length(effect) == npost) {
    eff_path <- effect
  } else {
    stop("effect must be length 1 or T - t0 + 1 (= ", npost, ").", call. = FALSE)
  }
  Tau <- matrix(0, N, Tt)
  Tau[treated_units, ti0:Tt] <- matrix(eff_path, n_treated, npost, byrow = TRUE)

  Yobs <- Y0 + Tau
  df <- data.frame(
    id = rep(units, times = Tt),
    t = rep(times, each = N),
    y = as.numeric(Yobs),
    w = as.numeric(W),
    y0 = as.numeric(Y0),
    tau = as.numeric(Tau),
    stringsAsFactors = FALSE
  )
  df[order(df$id, df$t), ]
}

#' True ATT of a (semi-)synthetic panel
#'
#' Convenience accessor for the ground-truth average treatment effect on the
#' treated of a simulated panel. It reads the per-cell true effect --- the `tau`
#' column produced by [sim_semisynthetic()], or `y - y0` for a [sim_panel()]
#' panel --- and averages it over the actively treated cells, giving the number
#' estimators are trying to recover. Pair it with [panel_compare()] or
#' [panel_rmse()] to score estimates against the truth without re-deriving it by
#' hand.
#'
#' @param data A long `data.frame` with a treatment indicator and either a `tau`
#'   column (true per-cell effect) or both `y` and `y0` columns.
#' @param treatment Name of the 0/1 treatment column (default `"w"`).
#' @param tau Name of the true-effect column (default `"tau"`). When that column
#'   is absent the per-cell effect is taken as `y - y0`.
#' @return A single numeric: the true ATT averaged over treated cells.
#' @seealso [sim_semisynthetic()], [sim_panel()], [panel_compare()]
#' @examples
#' real <- sim_panel(N = 40, T = 18, n_treated = 0L, att = 0, seed = 1)
#' ss <- sim_semisynthetic(real, "y", "id", "t",
#'                         n_treated = 6, t0 = 14, att = 3, seed = 2)
#' true_att(ss)   # 3
#' @export
true_att <- function(data, treatment = "w", tau = "tau") {
  data <- as.data.frame(data)
  if (!treatment %in% names(data))
    stop("Treatment column `", treatment, "` not found.", call. = FALSE)
  w <- data[[treatment]]
  if (tau %in% names(data)) {
    eff <- data[[tau]]
  } else if (all(c("y", "y0") %in% names(data))) {
    eff <- data[["y"]] - data[["y0"]]
  } else {
    stop("Need a `", tau, "` column, or both `y` and `y0`, to compute the ",
         "true ATT.", call. = FALSE)
  }
  treated <- w == 1 & is.finite(eff)
  if (!any(treated)) stop("No treated cells (", treatment, " == 1).", call. = FALSE)
  mean(eff[treated])
}
