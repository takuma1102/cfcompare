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

# Centre the main title and subtitle. Appended to every ggplot the package
# returns, so titles default to centre alignment everywhere.
#' @keywords internal
#' @noRd
.center_titles <- function() {
  ggplot2::theme(
    plot.title    = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5)
  )
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
    ggplot2::theme_minimal(base_size = 12) +
    .center_titles()
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
                                         colour = .data$series)) +
    ggplot2::geom_line(linewidth = 0.9)
  if (!is.na(t0_time)) {
    p <- p + ggplot2::geom_vline(xintercept = t0_time, linetype = "dashed",
                                 colour = "grey40", linewidth = 0.6)
  }
  p +
    ggplot2::labs(x = "Time", y = "Outcome (treated-unit average)",
                  colour = NULL,
                  title = "Observed vs. predicted counterfactual") +
    ggplot2::theme_minimal(base_size = 12) +
    .center_titles()
}

#' Synthetic-control-style trajectory plot for a single TROP fit
#'
#' Draws the treated-unit average observed path against the estimated untreated
#' (counterfactual) path, in the style of the \pkg{synthdid} plot: a dashed line
#' marks the first treated period, the post-treatment gap between the two lines is
#' the estimated effect, and -- since TROP carries explicit time weights -- the
#' time weights \eqn{\theta_s = \exp(-\lambda_{time} |t - s|)} are drawn as a
#' filled band along the bottom (as in \pkg{synthdid}'s time-weight strip) to
#' show which periods the counterfactual leans on. The band shows weights, not
#' observation counts; heights are scaled to the largest weight, and with
#' \eqn{\lambda_{time} = 0} the weights are uniform and the band is flat.
#'
#' This is a *point-estimate* trajectory: it shows no standard errors or
#' confidence intervals (the post-treatment gap is left unshaded so it is not
#' mistaken for a confidence band). To visualise the fit's uncertainty, use
#' [trop_event_study()] and call [autoplot()] on the result, which draws
#' per-period estimates with pointwise confidence bars. Accordingly the subtitle
#' reports only the selected penalties by default; set `show_se = TRUE` to also
#' note which standard-error method the fit used.
#'
#' @param object A `trop` fit from [trop()].
#' @param show_weights Logical; draw the time-weight band along the bottom.
#' @param show_se Logical; append the fit's standard-error method to the subtitle.
#'   Defaults to `FALSE` because this plot does not display standard errors --
#'   see [trop_event_study()] for the uncertainty visualisation.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @seealso [trop_event_study()] and its `autoplot()` method for per-period
#'   effects with confidence bars.
#' @export
#' @examples
#' \donttest{
#' df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 3, seed = 1)
#' autoplot(trop(df, "y", "w", "id", "t",
#'               control = trop_control(n_cv_cells = 8L, cv_cycles = 1L)))
#' }
autoplot.trop <- function(object, show_weights = TRUE, show_se = FALSE, ...) {
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
  M <- .trop_solve(Y, W, wmat, lam, trop_control(), X = object$panel$X)$M

  obs <- colMeans(Y[tu, , drop = FALSE])
  cf  <- colMeans(M[tu, , drop = FALSE])
  dat <- rbind(
    data.frame(time = times, value = obs, series = "Treated (observed)"),
    data.frame(time = times, value = cf,  series = "Estimated Y(0)"))

  t0_time <- if (!is.na(pat$block_t0)) times[pat$block_t0] else NA

  # subtitle: the three penalties, and (only when show_se = TRUE) how the SE was
  # obtained. By default the SE method is omitted, because this plot does not
  # display standard errors -- see trop_event_study() for the uncertainty view.
  # Built as a plotmath expression so the lambda and centre-dot glyphs are drawn
  # from the symbol font. Embedding the raw Unicode characters in device text
  # triggers an "mbcsToSbcs" conversion failure under non-UTF-8 locales (e.g.
  # the C locale used by R CMD check on some platforms).
  fmt <- function(z) if (is.infinite(z)) "Inf" else sprintf("%.2f", z)
  lam_txt <- sprintf("(unit %s, time %s, nn %s)",
                     fmt(lam$unit), fmt(lam$time), fmt(lam$nn))
  if (isTRUE(show_se)) {
    se_txt <- switch(object$se.method %||% "none",
      bootstrap = if (!is.null(object$n_boot))
                    sprintf("Bootstrap SE (%d reps)", object$n_boot) else "Bootstrap SE",
      jackknife = "Jackknife SE",
      none      = "no SE",
      sprintf("%s SE", object$se.method))
    sub_txt <- substitute("Penalties " * lambda == lt %.% st,
                          list(lt = lam_txt, st = se_txt))
  } else {
    sub_txt <- substitute("Penalties " * lambda == lt, list(lt = lam_txt))
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$time, y = .data$value,
                                         colour = .data$series)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(
      values = c("Treated (observed)" = "#1b9e9e",
                 "Estimated Y(0)" = "#e8665a")) +
    ggplot2::labs(x = "Time", y = "Outcome",
                  colour = NULL,
                  title = "TROP: treated vs. estimated counterfactual",
                  subtitle = sub_txt)

  if (!is.na(t0_time)) {
    # mark the first treated period. The post-treatment gap between the two
    # lines is the estimated effect; it is left unshaded so it is not mistaken
    # for a confidence band.
    p <- p +
      ggplot2::geom_vline(xintercept = t0_time, linetype = "dashed",
                          colour = "grey40", linewidth = 0.6)
  }

  if (show_weights && length(t_anchor)) {
    # TROP time weights theta_s = exp(-lambda_time * |t - s|): how much each
    # period is leaned on. Drawn as a filled band along the bottom in the style
    # of synthdid's time-weight strip (not observation counts). Heights are
    # scaled to the largest weight so the shape is readable; with lambda_time = 0
    # the weights are uniform and the band is flat.
    s <- seq_len(Tt)
    dtime <- vapply(s, function(ss) min(abs(ss - t_anchor)), numeric(1))
    theta <- exp(-lam$time * dtime)
    theta <- theta / max(theta)
    rng  <- range(c(obs, cf), na.rm = TRUE)
    base <- rng[1] - 0.14 * diff(rng)      # baseline below the trajectories
    h    <- 0.10 * diff(rng)               # band height at full weight
    wr <- data.frame(time = times, ymin = base, ymax = base + h * theta)
    p <- p +
      ggplot2::geom_ribbon(data = wr, inherit.aes = FALSE,
        ggplot2::aes(x = .data$time, ymin = .data$ymin, ymax = .data$ymax),
        fill = "#1b9e9e", alpha = 0.30) +
      ggplot2::geom_line(data = wr, inherit.aes = FALSE,
        ggplot2::aes(x = .data$time, y = .data$ymax),
        colour = "#1b9e9e", linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = base, colour = "grey75", linewidth = 0.3) +
      ggplot2::labs(caption = "Bottom band: TROP time weights")
  }

  p + ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top") +
    .center_titles()
}

#' @export
plot.trop <- function(x, ...) {
  print(autoplot.trop(x, ...))
  invisible(x)
}

#' Event-study plot for a TROP fit
#'
#' Plots the per-period treatment effects from [trop_event_study()] against event
#' time, with pointwise confidence-interval bars, a solid line at zero, and a
#' dotted line marking treatment onset. When pre-treatment periods are present
#' they are drawn as distinct placebo / pre-trend points (open markers) so the
#' pre-period profile can be read at a glance. The SE method is named at the end
#' of the subtitle.
#'
#' @param object A `trop_event_study` object from [trop_event_study()].
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
#' @examples
#' \donttest{
#' df  <- sim_panel(N = 24, T = 12, n_treated = 6, t0 = 8, att = 3, seed = 1)
#' fit <- trop(df, "y", "w", "id", "t",
#'             control = trop_control(n_cv_cells = 10L, cv_cycles = 1L))
#' autoplot(trop_event_study(fit, se = "bootstrap",
#'                           control = trop_control(n_boot = 100L, seed = 1)))
#' }
autoplot.trop_event_study <- function(object, ...) {
  .need_ggplot()
  est <- object$estimates
  has_ci <- any(is.finite(est$conf.low) & is.finite(est$conf.high))
  has_pre <- any(est$period == "pre")

  se_txt <- switch(object$se.method %||% "none",
    bootstrap = if (!is.null(object$n_boot))
                  sprintf("Bootstrap SE (%d reps)", object$n_boot) else "Bootstrap SE",
    jackknife = "Jackknife SE",
    none      = "no SE",
    sprintf("%s SE", object$se.method))
  ci_txt <- sprintf("%.0f%% pointwise CI", 100 * object$conf.level)
  sub_txt <- sprintf("Per-period ATT, %s; %s", ci_txt, se_txt)

  est$Period <- factor(ifelse(est$period == "pre",
                              "Pre-treatment (placebo)", "Post-treatment"),
                       levels = c("Pre-treatment (placebo)", "Post-treatment"))

  p <- ggplot2::ggplot(est, ggplot2::aes(x = .data$event_time,
                                         y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", colour = "grey60") +
    ggplot2::geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50")

  if (has_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high,
                   colour = .data$Period),
      width = 0.18, na.rm = TRUE)
  }

  p <- p +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Period,
                                     shape = .data$Period), size = 2.4) +
    ggplot2::scale_colour_manual(
      values = c("Pre-treatment (placebo)" = "grey55",
                 "Post-treatment" = "#1b9e9e"), drop = FALSE) +
    ggplot2::scale_shape_manual(
      values = c("Pre-treatment (placebo)" = 1, "Post-treatment" = 16),
      drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = est$event_time) +
    ggplot2::labs(
      x = "Event time (periods relative to treatment)",
      y = "ATT", colour = NULL, shape = NULL,
      title = "TROP event study", subtitle = sub_txt) +
    ggplot2::theme_minimal(base_size = 12) +
    .center_titles()

  if (has_pre) {
    p <- p + ggplot2::theme(legend.position = "top")
  } else {
    p <- p + ggplot2::guides(colour = "none", shape = "none")
  }
  p
}

#' @export
plot.trop_event_study <- function(x, ...) {
  print(autoplot.trop_event_study(x, ...))
  invisible(x)
}
