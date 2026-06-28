# Penalty-sensitivity grid for TROP: how the ATT estimate and the
# cross-validation loss vary across a grid of penalties. Written independently.

# ---- internal helpers -------------------------------------------------------

# The three TROP penalties, in canonical order.
.trop_penalties <- c("time", "unit", "nn")

# Which two penalties form the heatmap axes (x, y), and which is held fixed.
# Falls back to the historical default (x = nn, y = time, fixed = unit) for
# grid objects created before the `axes` attribute existed.
.grid_axes <- function(g) {
  ax <- attr(g, "axes")
  if (is.null(ax)) ax <- c(x = "nn", y = "time")
  ax
}
.grid_fixed <- function(g) {
  fx <- attr(g, "fixed")
  if (is.null(fx)) fx <- setdiff(.trop_penalties, .grid_axes(g))
  unname(fx)
}

# Plotmath axis label for a penalty, e.g. lambda[nn].
.lambda_lab <- function(pen) switch(pen,
  time = expression(lambda[time]),
  unit = expression(lambda[unit]),
  nn   = expression(lambda[nn]))

# Plotmath symbol (language object) for a penalty, for use inside substitute().
.lambda_sym <- function(pen) switch(pen,
  time = quote(lambda[time]),
  unit = quote(lambda[unit]),
  nn   = quote(lambda[nn]))

#' TROP penalty-sensitivity grid (heatmap data)
#'
#' Sweeps any two of the three TROP penalties --- \eqn{\lambda_{time}},
#' \eqn{\lambda_{unit}} and \eqn{\lambda_{nn}} --- against each other while
#' holding the third fixed, and records at each grid point both the resulting
#' ATT estimate and the leave-one-out cross-validation loss. This is the data
#' behind the diagnostic heatmap: colour the cells by CV loss to see where the
#' data-driven choice lands, and read the ATT off each cell to judge how
#' sensitive the estimate is to the penalties. Use `axes` to choose which two
#' penalties go on the x and y axes; the remaining (fixed) penalty is reported
#' in the heatmap subtitle.
#'
#' @param data A long `data.frame`, one row per unit-time.
#' @param outcome,treatment,unit,time Column names (strings).
#' @param lambda_time,lambda_nn,lambda_unit Penalty inputs. For a penalty that
#'   is on an axis (see `axes`), supply a numeric grid to sweep (or `NULL` for a
#'   data-driven default grid; the `lambda_nn` default keeps only the finite
#'   part, as the heatmap axis needs finite values). For the penalty held fixed,
#'   supply a single value, or `NULL` to choose it once by cross-validation and
#'   hold it fixed across the sweep.
#' @param axes Length-2 character vector naming the two penalties to sweep, as
#'   `c(x, y)`: the first is the x axis, the second the y axis. Each must be one
#'   of `"time"`, `"unit"`, `"nn"`. Defaults to `c("nn", "time")` (the
#'   historical layout). The third penalty is held fixed.
#' @param anchor Estimation anchor passed to the ATT computation
#'   (`"pooled"` by default for speed).
#' @param control A list of solver/CV controls from [trop_control()].
#' @param seed Optional integer seed for CV-cell sampling.
#' @param verbose Logical; print progress.
#' @return A `cf_trop_grid` (a `data.frame`) with columns `lambda_time`,
#'   `lambda_nn`, `lambda_unit`, `att`, `cv_loss`. The two swept penalties vary
#'   across rows; the fixed penalty is constant. The grid point minimising
#'   `cv_loss` (the data-driven choice) is stored in `attr(., "selected")`, and
#'   the chosen layout in `attr(., "axes")` (named `x`/`y`) and
#'   `attr(., "fixed")`.
#' @seealso [trop()], [autoplot.cf_trop_grid()]
#' @export
#' @examples
#' \donttest{
#' df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 1)
#' g <- trop_sensitivity(df, "y", "w", "id", "t",
#'                       lambda_time = c(0, 0.1, 0.5), lambda_nn = c(2, 5),
#'                       control = trop_control(n_cv_cells = 8L, cv_cycles = 1L))
#' autoplot(g)
#'
#' # Sweep the unit penalty against the nuclear-norm penalty instead, holding
#' # the time penalty fixed (chosen by CV):
#' g2 <- trop_sensitivity(df, "y", "w", "id", "t", axes = c("nn", "unit"),
#'                        control = trop_control(n_cv_cells = 8L, cv_cycles = 1L))
#' autoplot(g2)
#' }
trop_sensitivity <- function(data, outcome, treatment, unit, time,
                             lambda_time = NULL, lambda_nn = NULL,
                             lambda_unit = NULL, axes = c("nn", "time"),
                             anchor = "pooled",
                             control = trop_control(), seed = NULL,
                             verbose = FALSE) {
  pens <- .trop_penalties
  axes <- as.character(axes)
  if (length(axes) != 2L || !all(axes %in% pens) || anyDuplicated(axes))
    stop('`axes` must name two distinct penalties from "time", "unit", "nn" ',
         "(given as c(x, y)).")
  x_pen <- axes[1L]; y_pen <- axes[2L]
  fixed_pen <- setdiff(pens, axes)

  m <- .panel_to_matrices(data, outcome, treatment, unit, time)
  Y <- m$Y; W <- m$W
  pat <- .assignment_pattern(W)
  grids <- .trop_default_grids(Y, W)

  user_grids <- list(time = lambda_time, unit = lambda_unit, nn = lambda_nn)
  axis_vals <- function(pen) {
    g <- user_grids[[pen]]
    if (is.null(g)) {
      g <- grids[[pen]]
      if (pen == "nn") g <- g[is.finite(g)]  # heatmap axis needs finite values
    }
    sort(unique(g))
  }
  x_vals <- axis_vals(x_pen)
  y_vals <- axis_vals(y_pen)

  if (!is.null(seed)) control$seed <- seed
  cv_cells <- .sample_control_cells(W, control$n_cv_cells, control$seed)

  # value of the held-fixed penalty: user scalar, or CV-chosen once
  fixed_val <- user_grids[[fixed_pen]]
  if (is.null(fixed_val)) {
    sel <- .trop_select_lambda(Y, W, grids, control, cv_cells, FALSE)
    fixed_val <- sel[[fixed_pen]]
  } else {
    fixed_val <- fixed_val[1L]
  }

  rows <- list(); k <- 0L
  for (yv in y_vals) {
    for (xv in x_vals) {
      k <- k + 1L
      lam <- list(time = NA_real_, unit = NA_real_, nn = NA_real_)
      lam[[x_pen]] <- xv
      lam[[y_pen]] <- yv
      lam[[fixed_pen]] <- fixed_val
      cvl <- .trop_cv_Q(Y, W, lam, control, cv_cells)
      att <- .trop_att(Y, W, lam, control, anchor, pat)$att
      rows[[k]] <- data.frame(lambda_time = lam$time, lambda_nn = lam$nn,
                              lambda_unit = lam$unit, att = att,
                              cv_loss = cvl, stringsAsFactors = FALSE)
      if (verbose) message(sprintf("  %s=%.4g %s=%.4g  att=%.3f cv=%.5g",
                                   x_pen, xv, y_pen, yv, att, cvl))
    }
  }
  out <- do.call(rbind, rows)
  sel_row <- out[which.min(out$cv_loss), , drop = FALSE]
  attr(out, "selected") <- sel_row
  attr(out, "outcome") <- outcome
  attr(out, "anchor") <- anchor
  attr(out, "axes") <- c(x = x_pen, y = y_pen)
  attr(out, "fixed") <- fixed_pen
  class(out) <- c("cf_trop_grid", "data.frame")
  out
}

#' @export
print.cf_trop_grid <- function(x, digits = 4, ...) {
  s <- attr(x, "selected")
  ax <- .grid_axes(x); x_pen <- ax[["x"]]; y_pen <- ax[["y"]]
  fixed_pen <- .grid_fixed(x)
  xcol <- paste0("lambda_", x_pen); ycol <- paste0("lambda_", y_pen)
  fcol <- paste0("lambda_", fixed_pen)
  cat("TROP penalty-sensitivity grid\n")
  cat(sprintf("  %d x %d grid (lambda_%s x lambda_%s), lambda_%s = %.4g (fixed)\n",
              length(unique(x[[xcol]])), length(unique(x[[ycol]])),
              x_pen, y_pen, fixed_pen, x[[fcol]][1]))
  cat(sprintf("  CV-selected: lambda_%s=%.4g, lambda_%s=%.4g -> ATT=%.4f (CV=%.5g)\n",
              x_pen, s[[xcol]], y_pen, s[[ycol]], s$att, s$cv_loss))
  cat(sprintf("  ATT range over grid: [%.4f, %.4f]\n",
              min(x$att), max(x$att)))
  invisible(x)
}

#' Heatmap of the TROP penalty-sensitivity grid
#'
#' Cells are coloured by cross-validation loss (darker = better out-of-sample
#' fit) and annotated with the ATT estimate at that penalty pair; the
#' CV-selected cell is outlined. The two swept penalties (see the `axes`
#' argument of [trop_sensitivity()]) form the x and y axes, and the fixed
#' penalty is reported in the subtitle.
#'
#' @param object A `cf_trop_grid` from [trop_sensitivity()].
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
autoplot.cf_trop_grid <- function(object, ...) {
  .need_ggplot()
  d <- as.data.frame(object)
  ax <- .grid_axes(object); x_pen <- ax[["x"]]; y_pen <- ax[["y"]]
  fixed_pen <- .grid_fixed(object)
  xcol <- paste0("lambda_", x_pen); ycol <- paste0("lambda_", y_pen)
  fcol <- paste0("lambda_", fixed_pen)
  d$xf <- factor(round(d[[xcol]], 4))
  d$yf <- factor(round(d[[ycol]], 4))
  s <- attr(object, "selected")

  # report the penalty held fixed across the sweep. Built as a plotmath
  # expression so the lambda glyph renders via the symbol font (locale-portable).
  fv <- d[[fcol]][1L]
  fv_txt <- if (is.infinite(fv)) "Inf" else sprintf("%.3g", fv)
  sub_txt <- substitute(
    "fixed " * FS == v * ";  cell values = ATT; red outline = CV-selected",
    list(FS = .lambda_sym(fixed_pen), v = fv_txt))

  ggplot2::ggplot(d, ggplot2::aes(x = .data$xf, y = .data$yf,
                                  fill = .data$cv_loss)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$att)),
                       size = 3, colour = "grey15") +
    ggplot2::geom_tile(
      data = data.frame(xf = factor(round(s[[xcol]], 4), levels = levels(d$xf)),
                        yf = factor(round(s[[ycol]], 4), levels = levels(d$yf))),
      ggplot2::aes(x = .data$xf, y = .data$yf), inherit.aes = FALSE,
      fill = NA, colour = "red", linewidth = 1) +
    ggplot2::scale_fill_viridis_c(option = "D", direction = -1, name = "CV loss") +
    ggplot2::labs(
      x = .lambda_lab(x_pen), y = .lambda_lab(y_pen),
      title = "TROP Penalty Sensitivity",
      subtitle = sub_txt) +
    ggplot2::theme_minimal(base_size = 12) +
    .center_titles()
}

#' @export
plot.cf_trop_grid <- function(x, ...) {
  print(autoplot.cf_trop_grid(x, ...))
  invisible(x)
}

#' Base-R surface views of a TROP penalty-sensitivity grid
#'
#' Draws the cross-validation-loss surface and the ATT surface from a
#' [trop_sensitivity()] grid as two separate, full-width base-graphics plots,
#' each with a colour key (legend) on the right. Unlike [autoplot.cf_trop_grid()]
#' --- a single `ggplot2` heatmap that packs both into one panel --- this view
#' keeps the surfaces apart so neither is squeezed, and needs no `ggplot2`. The
#' surface axes follow the swept penalties recorded on the grid (see the `axes`
#' argument of [trop_sensitivity()]).
#'
#' @param grid A `cf_trop_grid` returned by [trop_sensitivity()].
#' @param which Which surface(s) to draw: `"both"` (default), `"cv_loss"`, or
#'   `"att"`. With `"both"` the two surfaces are drawn as separate figures; on an
#'   interactive screen device you are prompted between them (use a multi-page
#'   device such as `pdf()`, or call once per surface, to keep both).
#' @param ask Logical; when drawing both on an interactive device, prompt before
#'   the second figure. Defaults to `interactive()`.
#' @param ... Further arguments passed to [graphics::filled.contour()].
#' @return Invisibly, a list with the `cv_loss` and `att` surface matrices
#'   (rows = y-axis penalty, columns = x-axis penalty).
#' @seealso [trop_sensitivity()], [autoplot.cf_trop_grid()]
#' @examples
#' \donttest{
#' df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 1)
#' g <- trop_sensitivity(df, "y", "w", "id", "t",
#'                       lambda_time = c(0, 0.25, 1),
#'                       control = trop_control(n_cv_cells = 12L, cv_cycles = 1L))
#' plot_trop_surfaces(g, which = "cv_loss")
#' }
#' @export
plot_trop_surfaces <- function(grid, which = c("both", "cv_loss", "att"),
                               ask = interactive(), ...) {
  stopifnot(inherits(grid, "cf_trop_grid"))
  which <- match.arg(which)

  ax <- .grid_axes(grid); x_pen <- ax[["x"]]; y_pen <- ax[["y"]]
  xcol <- paste0("lambda_", x_pen); ycol <- paste0("lambda_", y_pen)

  # Reshape the long grid into (y-axis x x-axis) matrices without a formula, so
  # there are no undefined-global NOTEs.
  yv <- sort(unique(grid[[ycol]]))
  xv <- sort(unique(grid[[xcol]]))
  cv_surface  <- matrix(NA_real_, length(yv), length(xv), dimnames = list(yv, xv))
  att_surface <- cv_surface
  ix <- cbind(match(grid[[ycol]], yv), match(grid[[xcol]], xv))
  cv_surface[ix]  <- grid$cv_loss
  att_surface[ix] <- grid$att

  viridis <- function(n) grDevices::hcl.colors(n, "Viridis", rev = TRUE)
  bluered <- function(n) grDevices::hcl.colors(n, "Blue-Red")
  draw <- function(z, main, key, pal) {
    graphics::filled.contour(xv, yv, t(z), color.palette = pal,
      xlab = .lambda_lab(x_pen), ylab = .lambda_lab(y_pen),
      main = main,
      key.title = graphics::title(main = key, cex.main = 0.9, font.main = 1),
      ...)
  }

  if (which %in% c("both", "cv_loss")) {
    draw(cv_surface, "CV loss surface", "CV loss", viridis)
  }
  if (which == "both" && isTRUE(ask)) {
    op <- graphics::par(ask = TRUE)            # prompt before the 2nd figure only
    on.exit(graphics::par(op), add = TRUE)
  }
  if (which %in% c("both", "att")) {
    draw(att_surface, "ATT surface", "ATT", bluered)
  }
  invisible(list(cv_loss = cv_surface, att = att_surface))
}
