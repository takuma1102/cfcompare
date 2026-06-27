# Plotting for cfcompare comparisons. ggplot2 is an Import; methods degrade
# gracefully if it is somehow unavailable.

#' @keywords internal
#' @noRd
.need_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Install it with ",
         "install.packages('ggplot2').", call. = FALSE)
  }
}

#' Forest plot of ATT estimates across methods
#'
#' Point estimate and (where available) confidence interval for each method,
#' on a single axis -- the quick visual comparison this package is built for.
#'
#' @param object A `cf_comparison` (from [panel_compare()]) or a
#'   `cf_att_tbl`.
#' @param ... Unused.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \donttest{
#' df <- sim_panel(seed = 1)
#' cmp <- panel_compare(df, "y", "w", "id", "t",
#'                      methods = c("DID", "MC", "TROP"))
#' autoplot(cmp)
#' }
#' @importFrom ggplot2 autoplot
#' @export
autoplot.cf_comparison <- function(object, ...) {
  autoplot.cf_att_tbl(object$att, ...)
}

#' @rdname autoplot.cf_comparison
#' @export
autoplot.cf_att_tbl <- function(object, ...) {
  .need_ggplot()
  df <- as.data.frame(object)
  df <- df[!is.na(df$estimate), , drop = FALSE]
  df$method <- factor(df$method, levels = rev(unique(df$method)))
  has_ci <- any(is.finite(df$conf.low) & is.finite(df$conf.high))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$estimate, y = .data$method)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey60")
  if (has_ci) {
    p <- p + ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$conf.low, xmax = .data$conf.high),
      height = 0.18, na.rm = TRUE)
  }
  p +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::labs(x = "ATT estimate", y = NULL,
                  title = "Panel estimator comparison") +
    ggplot2::theme_minimal(base_size = 12)
}

#' @export
plot.cf_comparison <- function(x, ...) print(autoplot(x, ...))

#' Plot observed vs predicted counterfactual trajectories
#'
#' For each native engine in a comparison (DID, MC, TROP), draws the average
#' outcome over the treated units against the predicted control counterfactual,
#' to show how the methods extrapolate through the post-treatment period.
#'
#' @param x A `cf_comparison` from [panel_compare()].
#' @param methods Optional subset of methods to draw.
#' @return A \pkg{ggplot2} object.
#' @export
plot_counterfactual <- function(x, methods = NULL) {
  .need_ggplot()
  stopifnot(inherits(x, "cf_comparison"))
  Y <- x$panel$Y; W <- x$panel$W
  times <- x$panel$times
  treated_units <- x$pattern$treated_units
  cfs <- x$counterfactual
  if (!is.null(methods)) cfs <- cfs[intersect(methods, names(cfs))]
  cfs <- cfs[!vapply(cfs, is.null, logical(1))]
  if (!length(cfs)) stop("No counterfactual paths available to plot.",
                         call. = FALSE)

  obs <- colMeans(Y[treated_units, , drop = FALSE], na.rm = TRUE)
  base <- data.frame(time = times, value = obs, series = "observed",
                     stringsAsFactors = FALSE)
  pred <- do.call(rbind, lapply(names(cfs), function(nm) {
    M <- cfs[[nm]]
    val <- colMeans(M[treated_units, , drop = FALSE], na.rm = TRUE)
    data.frame(time = times, value = val, series = nm,
               stringsAsFactors = FALSE)
  }))
  dat <- rbind(base, pred)

  t0_time <- if (!is.na(x$pattern$block_t0)) times[x$pattern$block_t0] else NA

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$time, y = .data$value,
                                         colour = .data$series,
                                         linetype = .data$series)) +
    ggplot2::geom_line(linewidth = 0.9)
  if (!is.na(t0_time)) {
    p <- p + ggplot2::geom_vline(xintercept = t0_time, linetype = "dotted",
                                 colour = "grey50")
  }
  p +
    ggplot2::labs(x = "Time", y = "Outcome (treated-unit average)",
                  colour = NULL, linetype = NULL,
                  title = "Observed vs. predicted counterfactual") +
    ggplot2::theme_minimal(base_size = 12)
}

#' Synthetic-control-style trajectory plot for a single TROP fit
#'
#' Draws the treated-unit average observed path against the estimated untreated
#' (counterfactual) path, in the style of the \pkg{synthdid} plot: a dotted line
#' marks the first treated period, the post-treatment gap between the two lines is
#' the estimated effect, and -- since TROP carries explicit time weights -- the
#' time weights \eqn{\theta_s} are drawn as a ribbon along the bottom to show
#' which periods the counterfactual leans on.
#'
#' @param object A `trop` fit from [trop()].
#' @param show_weights Logical; draw the time-weight ribbon along the bottom.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
#' @examples
#' \donttest{
#' df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 3, seed = 1)
#' autoplot(trop(df, "y", "w", "id", "t",
#'               control = trop_control(n_cv_cells = 8L, cv_cycles = 1L)))
#' }
autoplot.trop <- function(object, show_weights = TRUE, ...) {
  .need_ggplot()
  Y <- object$panel$Y; W <- object$panel$W; pat <- object$pattern
  N <- nrow(Y); Tt <- ncol(Y)
  times <- suppressWarnings(as.numeric(colnames(Y)))
  if (anyNA(times)) times <- seq_len(Tt)
  tu <- pat$treated_units

  # full pooled counterfactual for the plot (regardless of the fit's anchor)
  lam <- object$lambda
  du <- .unit_distance_pooled(Y, W, tu)
  t_anchor <- sort(unique(which(W == 1, arr.ind = TRUE)[, 2]))
  wmat <- .trop_weight_matrix(du, t_anchor, Tt, lam)
  M <- .trop_solve(Y, W, wmat, lam, trop_control())$M

  obs <- colMeans(Y[tu, , drop = FALSE])
  cf  <- colMeans(M[tu, , drop = FALSE])
  dat <- rbind(
    data.frame(time = times, value = obs, series = "treated (observed)"),
    data.frame(time = times, value = cf,  series = "estimated Y(0)"))

  t0_time <- if (!is.na(pat$block_t0)) times[pat$block_t0] else NA

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$time, y = .data$value,
                                         colour = .data$series)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(
      values = c("treated (observed)" = "#1b9e9e",
                 "estimated Y(0)" = "#e8665a")) +
    ggplot2::labs(x = "Time", y = "Outcome (treated-unit average)",
                  colour = NULL,
                  title = "TROP: treated vs. estimated counterfactual",
                  subtitle = sprintf("ATT = %.3f%s", object$estimate,
                    if (is.finite(object$std.error))
                      sprintf("  (SE %.3f)", object$std.error) else ""))

  if (!is.na(t0_time)) {
    # shade the post-treatment gap (the estimated effect)
    post <- times >= t0_time
    gap <- data.frame(time = times[post], ymin = pmin(obs, cf)[post],
                      ymax = pmax(obs, cf)[post])
    p <- p +
      ggplot2::geom_ribbon(data = gap, inherit.aes = FALSE,
        ggplot2::aes(x = .data$time, ymin = .data$ymin, ymax = .data$ymax),
        fill = "grey50", alpha = 0.18) +
      ggplot2::geom_vline(xintercept = t0_time, linetype = "dotted",
                          colour = "grey50")
  }

  if (show_weights && lam$time > 0) {
    s <- seq_len(Tt)
    dtime <- vapply(s, function(ss) min(abs(ss - t_anchor)), numeric(1))
    theta <- exp(-lam$time * dtime)
    rng <- range(c(obs, cf), na.rm = TRUE)
    base <- rng[1] - 0.12 * diff(rng)
    h <- 0.10 * diff(rng)
    wr <- data.frame(time = times, ymin = base, ymax = base + h * theta)
    p <- p + ggplot2::geom_rect(data = wr, inherit.aes = FALSE,
      ggplot2::aes(xmin = .data$time - 0.45, xmax = .data$time + 0.45,
                   ymin = .data$ymin, ymax = .data$ymax),
      fill = "#1b9e9e", alpha = 0.35)
  }

  p + ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
}

#' @export
plot.trop <- function(x, ...) {
  print(autoplot.trop(x, ...))
  invisible(x)
}
