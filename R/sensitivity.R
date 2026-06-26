# Penalty-sensitivity grid for TROP: how the ATT estimate and the
# cross-validation loss vary across a grid of penalties. Written independently.

#' TROP penalty-sensitivity grid (heatmap data)
#'
#' Sweeps the time penalty \eqn{\lambda_{time}} against the nuclear-norm penalty
#' \eqn{\lambda_{nn}} (holding the unit penalty fixed) and records, at each grid
#' point, both the resulting ATT estimate and the leave-one-out cross-validation
#' loss. This is the data behind the diagnostic heatmap: colour the cells by CV
#' loss to see where the data-driven choice lands, and read the ATT off each cell
#' to judge how sensitive the estimate is to the penalties.
#'
#' @param data A long `data.frame`, one row per unit-time.
#' @param outcome,treatment,unit,time Column names (strings).
#' @param lambda_time,lambda_nn Numeric grids for the two swept penalties.
#'   Defaults are derived from the data scale; `lambda_nn` defaults to the finite
#'   part of the default grid (the heatmap axis needs finite values).
#' @param lambda_unit The fixed unit penalty. If `NULL`, it is chosen once by
#'   cross-validation and held fixed across the sweep.
#' @param anchor Estimation anchor passed to the ATT computation
#'   (`"pooled"` by default for speed).
#' @param control A list of solver/CV controls from [trop_control()].
#' @param seed Optional integer seed for CV-cell sampling.
#' @param verbose Logical; print progress.
#' @return A `cf_trop_grid` (a `data.frame`) with columns `lambda_time`,
#'   `lambda_nn`, `lambda_unit`, `att`, `cv_loss`. The grid point minimising
#'   `cv_loss` (the data-driven choice) is stored in `attr(., "selected")`.
#' @seealso [trop()], [autoplot.cf_trop_grid()]
#' @export
#' @examples
#' \donttest{
#' df <- sim_panel(N = 25, T = 14, n_treated = 4, t0 = 11, att = 2, seed = 1)
#' g <- trop_sensitivity(df, "y", "w", "id", "t",
#'                       control = trop_control(n_cv_cells = 40L, cv_cycles = 1L))
#' autoplot(g)
#' }
trop_sensitivity <- function(data, outcome, treatment, unit, time,
                             lambda_time = NULL, lambda_nn = NULL,
                             lambda_unit = NULL, anchor = "pooled",
                             control = trop_control(), seed = NULL,
                             verbose = FALSE) {
  m <- .panel_to_matrices(data, outcome, treatment, unit, time)
  Y <- m$Y; W <- m$W
  pat <- .assignment_pattern(W)
  grids <- .trop_default_grids(Y, W)
  if (is.null(lambda_time)) lambda_time <- grids$time
  if (is.null(lambda_nn)) lambda_nn <- grids$nn[is.finite(grids$nn)]
  lambda_time <- sort(unique(lambda_time))
  lambda_nn <- sort(unique(lambda_nn))

  if (!is.null(seed)) control$seed <- seed
  cv_cells <- .sample_control_cells(W, control$n_cv_cells, control$seed)

  if (is.null(lambda_unit)) {
    sel <- .trop_select_lambda(Y, W, grids, control, cv_cells, FALSE)
    lambda_unit <- sel$unit
  }

  rows <- list(); k <- 0L
  for (ln in lambda_nn) {
    for (lt in lambda_time) {
      k <- k + 1L
      lam <- list(time = lt, unit = lambda_unit, nn = ln)
      cvl <- .trop_cv_Q(Y, W, lam, control, cv_cells)
      att <- .trop_att(Y, W, lam, control, anchor, pat)$att
      rows[[k]] <- data.frame(lambda_time = lt, lambda_nn = ln,
                              lambda_unit = lambda_unit, att = att,
                              cv_loss = cvl, stringsAsFactors = FALSE)
      if (verbose) message(sprintf("  nn=%.4g time=%.4g  att=%.3f cv=%.5g",
                                   ln, lt, att, cvl))
    }
  }
  out <- do.call(rbind, rows)
  sel_row <- out[which.min(out$cv_loss), , drop = FALSE]
  attr(out, "selected") <- sel_row
  attr(out, "outcome") <- outcome
  attr(out, "anchor") <- anchor
  class(out) <- c("cf_trop_grid", "data.frame")
  out
}

#' @export
print.cf_trop_grid <- function(x, digits = 4, ...) {
  s <- attr(x, "selected")
  cat("TROP penalty-sensitivity grid\n")
  cat(sprintf("  %d x %d grid (lambda_time x lambda_nn), lambda_unit = %.4g\n",
              length(unique(x$lambda_time)), length(unique(x$lambda_nn)),
              x$lambda_unit[1]))
  cat(sprintf("  CV-selected: lambda_time=%.4g, lambda_nn=%.4g -> ATT=%.4f (CV=%.5g)\n",
              s$lambda_time, s$lambda_nn, s$att, s$cv_loss))
  cat(sprintf("  ATT range over grid: [%.4f, %.4f]\n",
              min(x$att), max(x$att)))
  invisible(x)
}

#' Heatmap of the TROP penalty-sensitivity grid
#'
#' Cells are coloured by cross-validation loss (darker = better out-of-sample
#' fit) and annotated with the ATT estimate at that penalty pair; the
#' CV-selected cell is outlined.
#'
#' @param object A `cf_trop_grid` from [trop_sensitivity()].
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
autoplot.cf_trop_grid <- function(object, ...) {
  .need_ggplot()
  d <- as.data.frame(object)
  d$ltf <- factor(round(d$lambda_time, 4))
  d$lnf <- factor(round(d$lambda_nn, 4))
  s <- attr(object, "selected")
  ggplot2::ggplot(d, ggplot2::aes(x = .data$lnf, y = .data$ltf,
                                  fill = .data$cv_loss)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$att)),
                       size = 3, colour = "grey15") +
    ggplot2::geom_tile(
      data = data.frame(lnf = factor(round(s$lambda_nn, 4), levels = levels(d$lnf)),
                        ltf = factor(round(s$lambda_time, 4), levels = levels(d$ltf))),
      ggplot2::aes(x = .data$lnf, y = .data$ltf), inherit.aes = FALSE,
      fill = NA, colour = "red", linewidth = 1) +
    ggplot2::scale_fill_viridis_c(option = "D", direction = -1, name = "CV loss") +
    ggplot2::labs(
      x = expression(lambda[nn]), y = expression(lambda[time]),
      title = "TROP penalty sensitivity",
      subtitle = "colour = CV loss; numbers = ATT; red = CV-selected") +
    ggplot2::theme_minimal(base_size = 12)
}

#' @export
plot.cf_trop_grid <- function(x, ...) {
  print(autoplot.cf_trop_grid(x, ...))
  invisible(x)
}
