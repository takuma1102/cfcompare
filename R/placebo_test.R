# Placebo test for a TROP fit: relabel control units as pseudo-treated (true
# effect zero), matching the real design's treated-unit count and timing, refit
# at the fitted penalties, and compare the observed ATT to the placebo null.
# This is the "is the estimate plausibly zero?" check -- distinct from a standard
# error: the placebo SE option has been removed from the package in favour of
# this test plus the bootstrap / jackknife intervals.

#' Placebo (randomisation) test for a TROP ATT
#'
#' Assesses whether a fitted [trop()] ATT is distinguishable from zero by
#' comparing it to a placebo null distribution. Control units are repeatedly
#' relabelled as pseudo-treated --- using the *same number of treated units and
#' the same treatment timing* as the real design, so block, staggered and
#' non-absorbing patterns are all reproduced faithfully --- and the estimator is
#' refit at the fit's cross-validated penalties (the true effect is zero by
#' construction). The observed ATT is then placed against that null.
#'
#' This is a randomisation/permutation test, not a standard error: it answers
#' "could a TROP design like mine produce an effect this large when there is none?"
#' For interval estimates use the bootstrap or jackknife (see [trop()] and
#' [compare_se_modes()]).
#'
#' @param object A `trop` fit from [trop()].
#' @param B Number of placebo assignments to draw.
#' @param alternative Direction of the test: `"two.sided"` (default), `"greater"`
#'   or `"less"`.
#' @param control A [trop_control()] list for the placebo refits (penalties are
#'   held fixed at the fit's values; this mainly sets `seed` and `workers`).
#'   Defaults to [trop_control()]; set `seed` here for reproducible draws.
#' @return An object of class `trop_placebo_test` (a list) with the `observed`
#'   ATT, the vector of `placebo` ATTs, the `p.value`, the placebo-null `mean`
#'   and `sd`, a central null interval `null.low`/`null.high`, `alternative`,
#'   `B`, and `n_treated_units`. It has `print()` and [autoplot()] methods.
#' @seealso [trop()], [compare_se_modes()]
#' @examples
#' \donttest{
#' df  <- sim_panel(N = 40, T = 12, n_treated = 5, t0 = 9, att = 3, seed = 1)
#' fit <- trop(df, "y", "w", "id", "t", se = "none",
#'             control = trop_control(n_cv_cells = 10L, cv_cycles = 1L))
#' pt  <- trop_placebo_test(fit, B = 200)
#' pt
#' }
#' @export
trop_placebo_test <- function(object, B = 500L,
                              alternative = c("two.sided", "greater", "less"),
                              control = NULL) {
  stopifnot(inherits(object, "trop"))
  alternative <- match.arg(alternative)
  ctrl <- control %||% trop_control()

  Y <- object$panel$Y; W <- object$panel$W; X <- object$panel$X
  lam0 <- object$lambda
  lam <- list(time = lam0$time, unit = lam0$unit, nn = lam0$nn)
  anchor <- object$anchor %||% "pooled"
  N <- nrow(Y); Tt <- ncol(Y)
  tu <- object$pattern$treated_units
  G <- length(tu)
  if (G < 1L) stop("The fit has no treated units.", call. = FALSE)
  controls <- setdiff(seq_len(N), tu)
  if (length(controls) < G + 1L)
    stop("Not enough control units (", length(controls),
         ") to assign ", G, " pseudo-treated unit(s).", call. = FALSE)

  # control-only donor panel (real treated outcomes are contaminated, so drop
  # them) and the real treated rows' time patterns, reused as the placebo design.
  Yc <- Y[controls, , drop = FALSE]; ncy <- nrow(Yc)
  Xc <- if (is.null(X)) NULL else
    lapply(.as_cov_list(X, N, Tt), function(M) M[controls, , drop = FALSE])
  treat_rows <- W[tu, , drop = FALSE]            # G x T
  observed <- object$estimate

  if (!is.null(ctrl$seed)) {
    old <- .Random.seed_safe(); on.exit(.Random.seed_restore(old), add = TRUE)
    set.seed(ctrl$seed)
  }
  draws <- lapply(seq_len(B), function(b) sample.int(ncy, G))
  par <- (ctrl$workers %||% 1L) > 1L
  placebo <- unlist(.par_lapply(draws, function(idx) {
    Wp <- matrix(0, ncy, Tt, dimnames = dimnames(Yc))
    Wp[idx, ] <- treat_rows
    patp <- .assignment_pattern(Wp)
    tryCatch(.trop_att(Yc, Wp, lam, ctrl, anchor, patp, X = Xc)$att,
             error = function(e) NA_real_)
  }, parallel = par), use.names = FALSE)
  # panel$Y is on the fitting scale; the observed ATT is on the raw scale, so
  # map the placebo ATTs back before comparing (identity unless the fit was
  # standardized).
  placebo <- (object$scaling %||% list(scale = 1))$scale * placebo
  placebo <- placebo[is.finite(placebo)]
  if (length(placebo) < 2L)
    stop("Placebo draws did not produce enough finite estimates.", call. = FALSE)

  p.value <- switch(alternative,
    two.sided = mean(abs(placebo) >= abs(observed)),
    greater   = mean(placebo >= observed),
    less      = mean(placebo <= observed))
  q <- stats::quantile(placebo, c(0.025, 0.975), names = FALSE, type = 7)

  structure(
    list(observed = observed, placebo = placebo, p.value = p.value,
         mean = mean(placebo), sd = stats::sd(placebo),
         null.low = q[1], null.high = q[2],
         alternative = alternative, B = length(placebo),
         n_treated_units = G, outcome = object$outcome %||% NA_character_),
    class = "trop_placebo_test")
}

#' @rdname trop_placebo_test
#' @param x A `trop_placebo_test` object.
#' @param digits Number of significant digits to print.
#' @param ... Unused.
#' @export
print.trop_placebo_test <- function(x, digits = 4, ...) {
  cat("TROP placebo test\n")
  cat(sprintf("  observed ATT : %s\n", format(x$observed, digits = digits)))
  cat(sprintf("  placebo null : mean %s, sd %s  (%d draws, %d pseudo-treated unit%s)\n",
              format(x$mean, digits = digits), format(x$sd, digits = digits),
              x$B, x$n_treated_units, if (x$n_treated_units == 1L) "" else "s"))
  cat(sprintf("  null 95%% range: [%s, %s]\n",
              format(x$null.low, digits = digits),
              format(x$null.high, digits = digits)))
  cat(sprintf("  %s p-value : %.4f\n", x$alternative, x$p.value))
  invisible(x)
}

#' Plot a TROP placebo test
#'
#' Histogram of the placebo-null ATT distribution with the observed ATT drawn as
#' a solid reference line and zero as a dashed line.
#'
#' @param object A `trop_placebo_test` from [trop_placebo_test()].
#' @param ... Unused.
#' @return A \pkg{ggplot2} object.
#' @export
autoplot.trop_placebo_test <- function(object, ...) {
  .need_ggplot()
  d <- data.frame(att = object$placebo)
  obs <- object$observed
  rng <- range(c(object$placebo, obs, 0), na.rm = TRUE)
  pad <- 0.08 * diff(rng); rng <- c(rng[1] - pad, rng[2] + pad)
  bins <- max(15L, min(40L, ceiling(sqrt(length(object$placebo)))))
  lab_h <- if (obs > mean(rng)) 1.06 else -0.06   # keep the label inside the panel
  ggplot2::ggplot(d, ggplot2::aes(x = .data$att)) +
    ggplot2::geom_histogram(bins = bins, fill = "grey82", colour = "white") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_vline(xintercept = obs, colour = "#b2182b", linewidth = 1) +
    ggplot2::annotate("text", x = obs, y = Inf, vjust = 1.4, hjust = lab_h,
                      label = sprintf("observed = %.3g", obs),
                      colour = "#b2182b", size = 3.4) +
    ggplot2::labs(
      x = "Placebo ATT (true effect = 0)", y = "count",
      title = "Placebo test: observed ATT vs the placebo null",
      subtitle = sprintf("%d placebo draws \u00b7 %s p = %.3f",
                         object$B, object$alternative, object$p.value)) +
    ggplot2::coord_cartesian(xlim = rng) +
    ggplot2::theme_minimal(base_size = 12) +
    .center_titles()
}

#' @export
plot.trop_placebo_test <- function(x, ...) print(autoplot(x, ...))
