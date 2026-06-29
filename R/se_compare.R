# Compare TROP inference (standard-error) methods on a single fit: choose the
# penalties once by cross-validation, then re-run each requested SE method on
# that fixed estimate so the point estimate is identical across rows and only
# the uncertainty differs.

#' Compare TROP standard-error methods on one fit
#'
#' Refits the TROP estimator under several inference methods --- the unit-level
#' stratified block bootstrap and/or the leave-one-treated-unit-out jackknife ---
#' and returns the estimate, standard error and confidence interval from each on
#' the common [`cf_att_tbl`][as_att] schema. The cross-validated penalties are
#' chosen once (a single `se = "none"` fit) and reused for every inference
#' method, exactly as [trop_ablation()] reuses one fit's penalties, so the ATT is
#' identical across rows and only the uncertainty changes. This is the table to
#' reach for when checking how sensitive the confidence interval is to the choice
#' of resampling scheme --- it is an *inference-mode* comparison, not an
#' estimator comparison, so all rows share one point estimate by construction.
#'
#' @param data,outcome,treatment,unit,time,covariates Passed to [trop()].
#' @param se Character vector of inference methods to compare; any of
#'   `"bootstrap"`, `"jackknife"`. Defaults to both. (To compare anchors or
#'   penalty constraints instead, see [trop()]'s `anchor` argument and
#'   [trop_ablation()].)
#' @param anchor Weight anchoring, as in [trop()] (default `"pooled"`).
#' @param control A [trop_control()] list, shared by every fit. Set `n_boot`
#'   here to control the bootstrap replications.
#' @param labels Optional `method` labels, one per `se`; defaults to
#'   `"Bootstrap SE"` / `"Jackknife SE"`.
#' @return An object of class `cf_se_comparison` (extending `cf_att_tbl`) with one
#'   row per inference method, sharing a common point estimate. If a method does
#'   not apply to the design (for example jackknife with a single treated unit) it
#'   still contributes a row, carrying the shared estimate with `NA` uncertainty
#'   and an explanatory `note`.
#' @seealso [trop()], [trop_ablation()], [bind_att()]
#' @examples
#' df <- sim_panel(N = 24, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 1)
#' \donttest{
#' compare_se_modes(df, "y", "w", "id", "t",
#'                  control = trop_control(n_cv_cells = 10L, cv_cycles = 1L))
#' }
#' @export
compare_se_modes <- function(data, outcome, treatment, unit, time,
                             covariates = NULL,
                             se = c("bootstrap", "jackknife"),
                             anchor = "pooled",
                             control = trop_control(),
                             labels = NULL) {
  se_choices <- c("bootstrap", "jackknife")
  if (missing(se)) se <- se_choices
  se <- unique(match.arg(se, se_choices, several.ok = TRUE))
  pretty <- c(bootstrap = "Bootstrap SE", jackknife = "Jackknife SE")
  labs <- labels %||% unname(pretty[se])
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
  class(out) <- c("cf_se_comparison", "cf_att_tbl", "data.frame")
  out
}

#' Plot a TROP standard-error comparison
#'
#' Forest plot of a [compare_se_modes()] result: every row carries the *same*
#' TROP point estimate (a vertical reference line marks it) and only the
#' confidence interval differs between resampling schemes. The title makes clear
#' this is an inference-mode comparison rather than an estimator comparison.
#'
#' @param object A `cf_se_comparison` from [compare_se_modes()].
#' @param ... Unused.
#' @return A \pkg{ggplot2} object.
#' @export
autoplot.cf_se_comparison <- function(object, ...) {
  .need_ggplot()
  df <- as.data.frame(object)
  df <- df[!is.na(df$estimate), , drop = FALSE]
  df$method <- factor(df$method, levels = rev(unique(df$method)))
  has_ci <- any(is.finite(df$conf.low) & is.finite(df$conf.high))
  est <- stats::median(df$estimate, na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$estimate, y = .data$method)) +
    ggplot2::geom_vline(xintercept = est, linetype = "dotted", colour = "grey55")
  if (has_ci) {
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$conf.low, xmax = .data$conf.high),
      height = 0.16, na.rm = TRUE)
  }
  p +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::labs(
      x = "ATT estimate", y = NULL,
      title = "Same TROP estimate under alternative SE methods",
      subtitle = "Point estimate fixed by cross-validation; only the interval differs") +
    ggplot2::theme_minimal(base_size = 12) +
    .center_titles()
}

#' @export
plot.cf_se_comparison <- function(x, ...) print(autoplot(x, ...))
