# Estimation-RMSE curves as a function of a design dimension (number of control
# units, or number of pre-treatment periods), one line per estimator. Reproduces
# the style of the paper's "RMSE vs. N_control / T_pre" diagnostic figures.
# Written independently; uses a semi-synthetic factor-model DGP with interactive
# confounding so that plain DID is biased and the robust estimators are not.

#' Sequential or parallel lapply. Uses future.apply::future_lapply when
#' `parallel = TRUE` and the package is installed (honouring the future plan set
#' by the user, e.g. future::plan(future::multisession)); otherwise base lapply.
#' Results are identical either way because each task is self-seeded.
#' @keywords internal
#' @noRd
.par_lapply <- function(X, FUN, parallel = FALSE) {
  if (isTRUE(parallel) && requireNamespace("future.apply", quietly = TRUE)) {
    return(future.apply::future_lapply(X, FUN, future.seed = TRUE))
  }
  lapply(X, FUN)
}

#' Evaluate `expr` with a transient parallel `future` plan.
#'
#' When `workers > 1` and the `future`/`future.apply` packages are installed, a
#' `multisession` plan with that many workers is set up for the duration of
#' `expr` and the previously active plan is restored afterwards. This lets the
#' embarrassingly parallel inner loops (which call `.par_lapply()`) share one
#' worker pool instead of paying start-up costs on every call. When `workers` is
#' `1`, or the packages are unavailable, `expr` runs unchanged in the current
#' process, so the code path still works with only base R installed.
#' @keywords internal
#' @noRd
.with_workers <- function(workers, expr) {
  workers <- if (is.null(workers)) 1L else workers
  if (workers > 1L &&
      requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)) {
    oplan <- future::plan()
    on.exit(future::plan(oplan), add = TRUE)
    future::plan(future::multisession, workers = workers)
  }
  force(expr)
}

#' Format a duration in seconds as a short human string for progress messages.
#' @keywords internal
#' @noRd
.fmt_dur <- function(s) {
  if (!is.finite(s)) return("?")
  if (s < 60) return(sprintf("%.0fs", s))
  if (s < 3600) return(sprintf("%dm%02ds", as.integer(s) %/% 60L,
                                as.integer(round(s)) %% 60L))
  sprintf("%dh%02dm", as.integer(s) %/% 3600L, (as.integer(s) %% 3600L) %/% 60L)
}

#' Default semi-synthetic generator: a richer factor model with interactive
#' confounding. Y(0) combines unit/time fixed effects, several random-walk latent
#' factors, heterogeneous unit-specific linear trends, and AR(1) idiosyncratic
#' errors. Treatment is selected on BOTH the factor loadings and the trend slope,
#' so plain DID/TWFE (and, partly, SC) are biased while the robust estimators are
#' not. Written independently.
#' @keywords internal
#' @noRd
.rmse_curve_gen <- function(seed, n_control, n_treated, n_pre, n_post,
                            rank, att, noise, ar = 0.4, trend_sd = 0.05) {
  set.seed(seed)
  N <- n_control + n_treated
  T <- n_pre + n_post
  t0 <- n_pre + 1L

  alpha <- stats::rnorm(N, sd = 0.5)              # unit fixed effects
  beta  <- cumsum(stats::rnorm(T, sd = 0.3))      # common time path

  # several random-walk latent factors with heterogeneous loadings
  F <- matrix(stats::rnorm(N * rank), N, rank)
  G <- matrix(stats::rnorm(T * rank), T, rank)
  for (r in seq_len(rank)) G[, r] <- cumsum(G[, r]) / sqrt(T)
  L <- F %*% t(G)

  # heterogeneous unit-specific linear trends (non-parallel trends)
  gamma <- stats::rnorm(N, sd = trend_sd)
  trend <- outer(gamma, seq_len(T))

  # AR(1) idiosyncratic errors
  E <- matrix(0, N, T)
  innov <- matrix(stats::rnorm(N * T), N, T) * noise
  E[, 1] <- innov[, 1]
  if (T > 1) for (s in 2:T) E[, s] <- ar * E[, s - 1] + sqrt(1 - ar^2) * innov[, s]

  Y0 <- outer(alpha, rep(1, T)) + outer(rep(1, N), beta) + L + trend + E

  # confounding: selection on loadings AND on the trend slope
  score <- F[, 1] + 0.5 * F[, min(2, rank)] + 8 * gamma + stats::rnorm(N, sd = 0.3)
  treated <- order(score, decreasing = TRUE)[seq_len(n_treated)]

  W <- matrix(0, N, T)
  W[treated, t0:T] <- 1
  data.frame(id = rep(seq_len(N), T), t = rep(seq_len(T), each = N),
             y = as.numeric(Y0 + att * W), w = as.numeric(W))
}

#' Estimation-RMSE curve over one design dimension
#'
#' Sweeps one design dimension -- the number of control units (`"n_control"`),
#' pre-treatment periods (`"n_pre"`), treated units (`"n_treated"`) or
#' post-treatment periods (`"n_post"`) -- and, at each value, runs a
#' semi-synthetic Monte Carlo: panels are drawn from a latent factor model in
#' which treatment is selected on the factor loadings (so plain DID/TWFE is
#' biased), a known constant effect `att` is imposed, every estimator is run, and
#' the estimation RMSE against the known truth is recorded. The result is the
#' data behind the paper-style "RMSE vs. N_control / T_pre" line plot, one line
#' per estimator. With a dense `values` grid and enough `n_runs` the curves are
#' smooth, as in the paper.
#'
#' This is a simulation diagnostic on synthetic data; it does not take a user
#' panel. To run a curve on your own data, draw panels with [sim_semisynthetic()]
#' and call [panel_compare()] in a loop.
#'
#' @param vary Which dimension to sweep: `"n_control"`, `"n_pre"`,
#'   `"n_treated"` or `"n_post"`.
#' @param values Integer vector of values for the swept dimension. Defaults to
#'   `seq(20, 45, by = 5)` for control units, `seq(5, 38, by = 2)` for
#'   pre-periods, and `seq(2, 12, by = 2)` for treated units or post-periods.
#'   Each value is an actual measured point (its own Monte Carlo); pass a
#'   denser/wider grid for smoother or longer curves.
#' @param n_runs Monte Carlo replications per value (higher = smoother; the paper
#'   uses ~1000). Default 500. **Note:** with the dense default grid and six
#'   estimators this is a large simulation (thousands of fits per panel); reduce
#'   `n_runs` and/or `values` for a quick look, or parallelise the outer loop.
#' @param methods Estimators to include; any subset of the [panel_compare()]
#'   methods. Default `c("DID", "SDID", "SC", "MC", "DIFP", "TROP")`, the six
#'   estimators compared in the paper. `gsynth`, `augsynth` and `CS` can be added
#'   if those packages are installed.
#' @param exclude Optional character vector of methods to drop from `methods`
#'   (e.g. `exclude = "DIFP"`). Unknown names are ignored with a warning.
#' @param n_control,n_treated,n_pre,n_post Base design; the dimension named in
#'   `vary` is overridden by `values`, the rest are held fixed. Defaults give a
#'   sizeable panel (60 controls, 8 treated, 16 pre- and 6 post-periods).
#' @param rank,att,noise Number of latent factors, the imposed (true) ATT, and
#'   the idiosyncratic-noise scale.
#' @param ar AR(1) coefficient of the idiosyncratic errors (serial correlation).
#' @param trend_sd SD of the heterogeneous unit-specific linear trend slopes
#'   (non-parallel trends); treatment is also selected on this slope.
#' @param anchor Estimation anchor for the native ATT (`"pooled"` by default for
#'   speed across the many fits).
#' @param control Solver/CV controls from [trop_control()].
#' @param seed Base integer seed (each replication uses a distinct offset).
#' @param parallel Logical; if `TRUE` and the `future.apply` package is
#'   installed, the Monte Carlo tasks are run in parallel honouring the active
#'   `future` plan (e.g. `future::plan(future::multisession, workers = 6)`).
#'   Results are identical to the sequential default because each task is
#'   self-seeded. Default `FALSE`.
#' @param progress Logical; print a rough progress indicator (completed fits, a
#'   percentage, elapsed time and an ETA) periodically while the curve runs.
#'   Helpful for the heavy large-`n_runs` jobs. Defaults to `interactive()`, so
#'   it shows in an interactive session and stays silent in scripts and tests.
#'   Messages go to `stderr()`; wrap in `suppressMessages()` to mute.
#' @param verbose Logical; additionally print a line as each swept value
#'   finishes. Independent of `progress`.
#' @return A `cf_rmse_curve` (a `data.frame`) with columns `method`, `x`
#'   (the swept value), `rmse`, `bias`, `n_runs`. The swept-dimension label is in
#'   `attr(., "vary")`.
#' @seealso [rmse_curves()], [panel_rmse()], [autoplot.cf_rmse_curve()]
#' @export
#' @examples
#' \donttest{
#' # quick look (small grid + few reps); raise n_runs/values for paper quality
#' cc <- rmse_curve("n_control", values = c(12, 16), n_runs = 2, n_pre = 6, n_post = 3,
#'                  methods = c("DID", "TROP"),
#'                  control = trop_control(n_cv_cells = 5L, cv_cycles = 1L))
#' autoplot(cc)
#' }
rmse_curve <- function(vary = c("n_control", "n_pre", "n_treated", "n_post"),
                       values = NULL, n_runs = 500L,
                       methods = c("DID", "SDID", "SC", "MC", "DIFP", "TROP"),
                       exclude = NULL,
                       n_control = 60L, n_treated = 8L,
                       n_pre = 16L, n_post = 6L,
                       rank = 4L, att = 2, noise = 1,
                       ar = 0.4, trend_sd = 0.05,
                       anchor = "pooled", control = trop_control(),
                       seed = 1L, parallel = FALSE, progress = interactive(),
                       verbose = FALSE) {
  vary <- match.arg(vary)
  # rmse_curve() seeds each Monte Carlo replication internally (via
  # .rmse_curve_gen) so its results are reproducible; save and restore the
  # caller's RNG state so the function leaves the global stream untouched.
  old_seed <- .Random.seed_safe()
  on.exit(.Random.seed_restore(old_seed), add = TRUE)
  methods <- .resolve_methods(
    methods, exclude,
    c("DID", "SDID", "SC", "MC", "DIFP", "TROP", "gsynth", "augsynth", "CS"))
  if (is.null(values)) {
    values <- switch(vary,
      n_control = seq(20L, 45L, by = 5L),
      n_pre     = seq(5L, 38L, by = 2L),
      n_treated = seq(2L, 12L, by = 2L),
      n_post    = seq(2L, 12L, by = 2L))
  }
  values <- sort(unique(as.integer(values)))

  # One Monte Carlo replication for a given (value, run). Each task is fully
  # determined by its seed (.rmse_curve_gen re-seeds), so results are identical
  # whether tasks run sequentially or in parallel.
  one_task <- function(task) {
    v <- task$v; r <- task$r
    nc  <- if (vary == "n_control") v else n_control
    np  <- if (vary == "n_pre")     v else n_pre
    nt  <- if (vary == "n_treated") v else n_treated
    npo <- if (vary == "n_post")    v else n_post
    df <- .rmse_curve_gen(seed * 100000L + v * 100L + r,
                          n_control = nc, n_treated = nt,
                          n_pre = np, n_post = npo,
                          rank = rank, att = att, noise = noise,
                          ar = ar, trend_sd = trend_sd)
    cmp <- tryCatch(
      suppressMessages(panel_compare(df, "y", "w", "id", "t",
        methods = methods, se = "none", control = control, anchor = anchor)),
      error = function(e) NULL)
    est <- stats::setNames(rep(NA_real_, length(methods)), methods)
    if (!is.null(cmp)) est[cmp$att$method] <- cmp$att$estimate
    est
  }

  grid <- expand.grid(r = seq_len(n_runs), vi = seq_along(values),
                      KEEP.OUT.ATTRS = FALSE)
  tasks <- lapply(seq_len(nrow(grid)),
                  function(i) list(v = values[grid$vi[i]], r = grid$r[i]))

  # Run the tasks in chunks so a rough progress indicator (completed count, %,
  # elapsed time, ETA) can be printed periodically -- these runs get heavy when
  # n_runs is large. Chunking does not change results: every task is self-seeded.
  total <- length(tasks)
  res <- vector("list", total)
  if (isTRUE(progress)) {
    message(sprintf("rmse_curve [%s]: %d runs x %d value(s) = %d fits%s",
                    vary, n_runs, length(values), total,
                    if (isTRUE(parallel)) " (parallel)" else ""))
  }
  n_chunks <- max(1L, min(20L, total))
  bounds <- floor(seq(0, total, length.out = n_chunks + 1L))
  t0 <- Sys.time()
  for (ci in seq_len(n_chunks)) {
    lo <- bounds[ci] + 1L; hi <- bounds[ci + 1L]
    if (hi < lo) next                                # skip any empty chunk
    res[lo:hi] <- .par_lapply(tasks[lo:hi], one_task, parallel = parallel)
    if (isTRUE(progress)) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      eta <- if (hi < total) el / hi * (total - hi) else 0
      message(sprintf("  %d/%d fits (%2.0f%%) | elapsed %s | ETA %s",
                      hi, total, 100 * hi / total, .fmt_dur(el), .fmt_dur(eta)))
    }
  }

  rows <- list(); k <- 0L
  for (vi in seq_along(values)) {
    est <- do.call(rbind, res[grid$vi == vi])      # n_runs x methods, run order
    rmse <- sqrt(colMeans((est - att)^2, na.rm = TRUE))
    bias <- colMeans(est - att, na.rm = TRUE)
    for (m in methods) {
      k <- k + 1L
      rows[[k]] <- data.frame(method = m, x = values[vi], rmse = unname(rmse[m]),
                              bias = unname(bias[m]), n_runs = n_runs,
                              stringsAsFactors = FALSE)
    }
    if (verbose) message(sprintf("  %s = %d done", vary, values[vi]))
  }
  out <- do.call(rbind, rows)
  attr(out, "vary") <- vary
  attr(out, "xlab") <- switch(vary,
    n_control = "Number of control units",
    n_pre     = "Number of pre-treatment periods",
    n_treated = "Number of treated units",
    n_post    = "Number of post-treatment periods")
  attr(out, "att") <- att
  class(out) <- c("cf_rmse_curve", "data.frame")
  out
}

#' Estimation-RMSE curves over both design dimensions
#'
#' Convenience wrapper that runs [rmse_curve()] for both the number of control
#' units and the number of pre-treatment periods and bundles the two curves, so
#' they can be plotted individually (default) or side by side (`combined = TRUE`),
#' as in the paper's two-panel figure.
#'
#' @param values_control,values_pre Optional grids for each sweep (see
#'   [rmse_curve()] defaults).
#' @param ... Passed to [rmse_curve()] (`n_runs`, `methods`, base design, `att`,
#'   `noise`, `control`, `seed`, `parallel`, `verbose`, ...). Set
#'   `parallel = TRUE` (plus a `future::plan`) to run both sweeps in parallel.
#' @return A `cf_rmse_curves` object: a list with `$n_control` and `$n_pre`,
#'   each a `cf_rmse_curve`.
#' @seealso [rmse_curve()], [autoplot.cf_rmse_curves()]
#' @export
#' @examples
#' \donttest{
#' # defaults are a large simulation; use a small grid + few reps for a quick look
#' g <- rmse_curves(values_control = c(12, 16), values_pre = c(5, 8), n_post = 3,
#'                  n_runs = 2, methods = c("DID", "TROP"),
#'                  control = trop_control(n_cv_cells = 5L, cv_cycles = 1L))
#' plot(g)                  # two separate figures (default)
#' # combined = TRUE needs the optional 'patchwork' package:
#' if (requireNamespace("patchwork", quietly = TRUE))
#'   plot(g, combined = TRUE) # side-by-side, paper-style
#' }
rmse_curves <- function(values_control = NULL, values_pre = NULL, ...) {
  cc <- rmse_curve("n_control", values = values_control, ...)
  cp <- rmse_curve("n_pre", values = values_pre, ...)
  structure(list(n_control = cc, n_pre = cp), class = "cf_rmse_curves")
}

#' @export
print.cf_rmse_curve <- function(x, ...) {
  vary <- attr(x, "vary")
  cat("RMSE curve over", attr(x, "xlab"), "\n")
  cat(sprintf("  %d values of %s x %d methods, %d MC reps each\n",
              length(unique(x$x)), vary, length(unique(x$method)), x$n_runs[1]))
  best <- tapply(x$rmse, x$method, mean, na.rm = TRUE)
  cat("  mean RMSE by method (low = good):\n")
  ord <- sort(best)
  for (m in names(ord)) cat(sprintf("    %-6s %.3f\n", m, ord[[m]]))
  invisible(x)
}

#' @export
print.cf_rmse_curves <- function(x, ...) {
  cat("Two estimation-RMSE curves (control units; pre-treatment periods)\n\n")
  print(x$n_control); cat("\n"); print(x$n_pre)
  invisible(x)
}

# Shared method colour palette so both panels match.
#' @keywords internal
#' @noRd
.curve_palette <- function(methods) {
  base <- c(DID = "#d1495b", SDID = "#edae49", SC = "#66a61e", MC = "#7570b3",
            DIFP = "#8c613c", TROP = "#2c7fb8", gsynth = "#e377c2",
            augsynth = "#17becf", CS = "#999999")
  miss <- setdiff(methods, names(base))
  if (length(miss)) base[miss] <- grDevices::hcl.colors(length(miss), "Dark 3")
  base[methods]
}

#' Line plot of estimation RMSE versus a design dimension
#'
#' One line per estimator; the y axis is on a log scale (as in the paper).
#' Lower is better.
#'
#' @param object A `cf_rmse_curve` from [rmse_curve()].
#' @param log_y Logical; log10 y axis (default `TRUE`).
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
autoplot.cf_rmse_curve <- function(object, log_y = TRUE, ...) {
  .need_ggplot()
  d <- as.data.frame(object)
  d <- d[is.finite(d$rmse), , drop = FALSE]
  pal <- .curve_palette(unique(d$method))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$rmse,
                                       colour = .data$method)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.7) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::labs(
      x = attr(object, "xlab"), y = "Estimation RMSE of the ATT",
      title = sprintf("RMSE vs. %s", tolower(attr(object, "xlab"))),
      subtitle = sprintf("%d Monte Carlo reps per point", d$n_runs[1])) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "right") +
    .center_titles()
  if (isTRUE(log_y)) p <- p + ggplot2::scale_y_log10()
  p
}

#' Plot one or both estimation-RMSE curves
#'
#' @param x A `cf_rmse_curve` or `cf_rmse_curves`.
#' @param combined For `cf_rmse_curves`: `FALSE` (default) draws the two panels
#'   separately; `TRUE` draws them side by side (paper-style).
#' @param log_y Logical; log10 y axis (default `TRUE`).
#' @param file Optional path. If given, the figure is written to a PNG at
#'   `width` x `height` inches; for separate panels two files are written with
#'   `-control`/`-pre` suffixes. If `NULL`, draws to the current device.
#' @param width,height Figure size in inches. Defaults are deliberately generous
#'   (single 8x5; combined 12x5).
#' @param ... Unused.
#' @return The input, invisibly.
#' @export
plot.cf_rmse_curve <- function(x, log_y = TRUE, file = NULL,
                               width = 8, height = 5, ...) {
  .render_curve(x, log_y = log_y, file = file, width = width, height = height)
  invisible(x)
}

#' @rdname plot.cf_rmse_curve
#' @export
plot.cf_rmse_curves <- function(x, combined = FALSE, log_y = TRUE, file = NULL,
                                width = NULL, height = NULL, ...) {
  if (combined) {
    if (is.null(width)) width <- 12; if (is.null(height)) height <- 5
    .render_curve(x, combined = TRUE, log_y = log_y, file = file,
                  width = width, height = height)
  } else {
    if (is.null(width)) width <- 8; if (is.null(height)) height <- 5
    f1 <- f2 <- NULL
    if (!is.null(file)) {
      stem <- sub("\\.png$", "", file)
      f1 <- paste0(stem, "-control.png"); f2 <- paste0(stem, "-pre.png")
    }
    .render_curve(x$n_control, log_y = log_y, file = f1, width = width, height = height)
    .render_curve(x$n_pre, log_y = log_y, file = f2, width = width, height = height)
  }
  invisible(x)
}

#' Autoplot for paired RMSE curves
#'
#' @param object A `cf_rmse_curves`.
#' @param combined `FALSE` (default) returns a named list of two `ggplot`s;
#'   `TRUE` returns a single side-by-side `ggplot` (needs the `patchwork`
#'   package).
#' @param log_y Logical; log10 y axis.
#' @param ... Unused.
#' @return A list of two `ggplot`s, or one combined `ggplot`.
#' @export
autoplot.cf_rmse_curves <- function(object, combined = FALSE, log_y = TRUE, ...) {
  p1 <- autoplot.cf_rmse_curve(object$n_control, log_y = log_y)
  p2 <- autoplot.cf_rmse_curve(object$n_pre, log_y = log_y)
  if (!combined) return(list(n_control = p1, n_pre = p2))
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("install 'patchwork' for combined = TRUE, or use combined = FALSE.",
         call. = FALSE)
  }
  patchwork::wrap_plots(p1, p2, ncol = 2, guides = "collect")
}

# Internal renderer: ggplot when available, base-R otherwise; optional PNG file.
#' @keywords internal
#' @noRd
.render_curve <- function(obj, combined = FALSE, log_y = TRUE, file = NULL,
                          width = 8, height = 5) {
  have_gg <- requireNamespace("ggplot2", quietly = TRUE)
  open_dev <- !is.null(file)
  if (open_dev) grDevices::png(file, width = width, height = height,
                               units = "in", res = 130)
  on.exit(if (open_dev) grDevices::dev.off(), add = TRUE)

  if (have_gg) {
    if (inherits(obj, "cf_rmse_curves")) {
      pr <- autoplot.cf_rmse_curves(obj, combined = TRUE, log_y = log_y)
    } else {
      pr <- autoplot.cf_rmse_curve(obj, log_y = log_y)
    }
    print(pr); return(invisible())
  }
  # ---- base-R fallback ----
  draw_one <- function(d) {
    xl <- attr(d, "xlab") %||% ""
    d <- d[is.finite(d$rmse), , drop = FALSE]
    ms <- unique(d$method); pal <- .curve_palette(ms)
    yl <- range(d$rmse, na.rm = TRUE)
    plot(NA, xlim = range(d$x), ylim = yl, log = if (log_y) "y" else "",
         xlab = xl, ylab = "Estimation RMSE of the ATT",
         main = sprintf("RMSE vs. %s", tolower(xl)))
    for (m in ms) {
      di <- d[d$method == m, ]; di <- di[order(di$x), ]
      graphics::lines(di$x, di$rmse, col = pal[[m]], lwd = 2)
      graphics::points(di$x, di$rmse, col = pal[[m]], pch = 19, cex = 0.6)
    }
    graphics::legend("topright", legend = ms, col = pal[ms], lwd = 2,
                     bty = "n", cex = 0.85)
  }
  if (inherits(obj, "cf_rmse_curves")) {
    op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op), add = TRUE)
    draw_one(obj$n_control); draw_one(obj$n_pre)
  } else {
    draw_one(obj)
  }
  invisible()
}
