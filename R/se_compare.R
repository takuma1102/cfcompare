# Compare TROP inference (standard-error) methods on a single fit: choose the
# penalties once by cross-validation, then re-run each requested SE method on
# that fixed estimate so the point estimate is identical across rows and only
# the uncertainty differs.

#' Compare TROP standard-error methods on one fit
#'
#' Refits the TROP estimator under several inference methods --- bootstrap,
#' jackknife and/or placebo --- and returns the estimate, standard error and
#' confidence interval from each on the common [`cf_att_tbl`][as_att] schema. The
#' cross-validated penalties are chosen once (a single `se = "none"` fit) and
#' reused for every inference method, exactly as [trop_ablation()] reuses one
#' fit's penalties, so the ATT is identical across rows and only the uncertainty
#' changes. That is the table to reach for when checking how sensitive the
#' confidence interval is to the choice of resampling scheme.
#'
#' @param data,outcome,treatment,unit,time,covariates Passed to [trop()].
#' @param se Character vector of inference methods to compare; any of
#'   `"bootstrap"`, `"jackknife"`, `"placebo"`. Defaults to all three. (To
#'   compare anchors or penalty constraints instead, see [trop()]'s `anchor`
#'   argument and [trop_ablation()].)
#' @param anchor Weight anchoring, as in [trop()] (default `"pooled"`).
#' @param control A [trop_control()] list, shared by every fit. Set `n_boot`
#'   here to control the bootstrap replications.
#' @param labels Optional `method` labels, one per `se`; defaults to
#'   `paste0("TROP_", se)`.
#' @return A `cf_att_tbl` (a `data.frame`) with one row per inference method,
#'   sharing a common point estimate. If a method does not apply to the design
#'   (for example jackknife with a single treated unit) it still contributes a
#'   row, carrying the shared estimate with `NA` uncertainty and an explanatory
#'   `note`.
#' @seealso [trop()], [trop_ablation()], [bind_att()]
#' @examples
#' df <- sim_panel(N = 24, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 1)
#' \donttest{
#' compare_se_modes(df, "y", "w", "id", "t",
#'                  se = c("jackknife", "placebo"),
#'                  control = trop_control(n_cv_cells = 10L, cv_cycles = 1L))
#' }
#' @export
compare_se_modes <- function(data, outcome, treatment, unit, time,
                             covariates = NULL,
                             se = c("bootstrap", "jackknife", "placebo"),
                             anchor = "pooled",
                             control = trop_control(),
                             labels = NULL) {
  se_choices <- c("bootstrap", "jackknife", "placebo")
  if (missing(se)) se <- se_choices
  se <- unique(match.arg(se, se_choices, several.ok = TRUE))
  labs <- labels %||% paste0("TROP_", se)
  if (length(labs) != length(se))
    stop("`labels` must have the same length as `se` (", length(se), ").",
         call. = FALSE)

  # One cross-validated fit supplies both the point estimate and the penalties
  # reused (fixed) by every inference method below.
  fit0 <- trop(data, outcome, treatment, unit, time,
               covariates = covariates, anchor = anchor,
               se = "none", control = control)
  lam0 <- fit0$lambda
  lam <- list(time = lam0$time, unit = lam0$unit, nn = lam0$nn)

  one <- function(s, lb) {
    f <- tryCatch(
      trop(data, outcome, treatment, unit, time,
           covariates = covariates, lambda = lam,
           anchor = anchor, se = s, control = control),
      error = function(e) e)
    if (inherits(f, "error")) {
      row <- as_att(fit0, method = lb)
      row$std.error <- NA_real_
      row$conf.low  <- NA_real_
      row$conf.high <- NA_real_
      row$note <- paste0("SE '", s, "' failed: ", conditionMessage(f))
      row
    } else {
      as_att(f, method = lb)
    }
  }

  out <- .rbind_att(Map(one, se, labs))
  attr(out, "outcome") <- outcome
  attr(out, "anchor")  <- anchor
  out
}
