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
#' @return A long `data.frame` with columns `id`, `t`, `y`, `w`, and the noiseless
#'   counterfactual `y0` (useful for evaluating estimators).
#' @examples
#' df <- sim_panel(N = 30, T = 15, n_treated = 5, t0 = 11, seed = 42)
#' head(df)
#' @export
sim_panel <- function(N = 30, T = 20, n_treated = 5, t0 = NULL,
                      rank = 3L, att = 1, noise = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(t0)) t0 <- floor(0.75 * T) + 1L
  stopifnot(n_treated < N, t0 >= 2, t0 <= T)

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

  W <- matrix(0, N, T)
  treated_units <- seq_len(n_treated)
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
  df[order(df$id, df$t), ]
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
  if (!is.null(seed)) set.seed(seed)
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
