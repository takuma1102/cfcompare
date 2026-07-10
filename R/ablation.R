# Penalty ablation for TROP: refit the estimator under a sequence of penalty
# constraints (the "Table 5"-style robustness check of the paper -- constraining
# lambda_time and/or lambda_unit to zero, and/or lambda_nn to infinity) and
# return the results as a publication-ready table rather than a plot.

#' Penalty ablation for the TROP estimator
#'
#' Refits TROP under a sequence of penalty constraints to show how the ATT moves
#' as the estimator is stripped back towards matrix completion and
#' difference-in-differences. This is the robustness exercise behind Table 5 of
#' Athey, Imbens, Qu & Viviano (2026): constraining `lambda_unit` and/or
#' `lambda_time` to zero (dropping unit/time weights) and/or `lambda_nn` to
#' infinity (dropping the low-rank regression adjustment).
#'
#' The full specification's penalties are taken from a single cross-validated
#' [trop()] fit; every constrained specification then reuses those values where
#' they are not constrained, so the rows are directly comparable. The result is a
#' table object: print it for a clean console table, or [format()] it to
#' paste-ready LaTeX (booktabs) or Markdown. It is also an ordinary
#' `data.frame`, so it can be passed to `knitr::kable()`, `gt::gt()`, etc.
#'
#' @param data,outcome,treatment,unit,time,covariates Passed to [trop()].
#' @param control A [trop_control()] list (used for every fit).
#' @param anchor Weight anchoring, as in [trop()] (default `"pooled"`).
#' @param se Standard-error method for every specification (default `"none"`;
#'   any [trop()] method is allowed, e.g. `"jackknife"`). When not `"none"`, the
#'   standard error and confidence interval appear in the table.
#' @param lambda_full Optional named list `list(time=, unit=, nn=)` to fix the
#'   full specification's penalties and skip cross-validation.
#' @param rmse Logical; if `TRUE`, also score each specification by a placebo
#'   RMSE. A pseudo block-treatment (true effect zero) is assigned to control
#'   units, every specification is refit at its penalties, and the RMSE is
#'   `sqrt(mean(estimate^2))` over `n_runs` draws. The same draws are scored under
#'   every specification, so the relative RMSE with TROP (full) = 1 is the
#'   quantity reported in Table 5 of the paper. Adds `rmse` and `rmse_rel`
#'   columns, and makes `plot()` show the (relative) RMSE instead of the ATT.
#' @param horizon,n_pseudo Placebo design used when `rmse = TRUE`: the number of
#'   final periods to pseudo-treat and the number of pseudo-treated control
#'   units. Default to the real design's post length and treated-unit count.
#' @param n_runs Number of placebo draws when `rmse = TRUE`.
#' @return An object of class `trop_ablation` (a `data.frame`) with one row per
#'   specification: `spec`, the three penalties `lambda_time`/`lambda_unit`/
#'   `lambda_nn`, `estimate`, `std.error`, `conf.low`/`conf.high`, and `rank`.
#'   When `rmse = TRUE` it also has `rmse` and `rmse_rel` (RMSE relative to the
#'   full TROP).
#' @references Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026).
#'   Triply Robust Panel Estimators. \emph{Journal of Applied Econometrics}. \doi{10.1002/jae.70061}.
#' @seealso [trop()]; `format()` for paste-ready LaTeX/Markdown output
#' @examples
#' df <- sim_panel(N = 14, T = 8, n_treated = 3, t0 = 5, att = 2, seed = 1)
#' # fixed base penalties keep the example fast; omit `lambda_full` to select
#' # them by cross-validation
#' ab <- trop_ablation(df, "y", "w", "id", "t",
#'                     lambda_full = list(time = 0.1, unit = 0.5, nn = 2))
#' ab                                   # clean console table
#' \donttest{
#' cat(format(ab, "latex"), sep = "\n") # paste-ready LaTeX (booktabs)
#' }
#' @export
trop_ablation <- function(data, outcome, treatment, unit, time,
                          covariates = NULL,
                          control = trop_control(),
                          anchor = "pooled",
                          se = "none",
                          lambda_full = NULL,
                          rmse = FALSE,
                          horizon = NULL,
                          n_pseudo = NULL,
                          n_runs = 200L) {
  fit_full <- trop(data, outcome, treatment, unit, time,
                   covariates = covariates, lambda = lambda_full,
                   anchor = anchor, se = se, control = control)
  lam <- fit_full$lambda
  nn_finite <- if (is.finite(lam$nn)) lam$nn else
    stats::sd(data[[outcome]], na.rm = TRUE)
  # Every spec carries the base fit's resolved nuclear-norm stabilising floor
  # (0 when disabled), so the trop() refits and the placebo-RMSE .trop_att()
  # solves below all use the same effective lambda_nn convention as the full
  # fit. It is a no-op for the lambda_nn = Inf rows.
  nn_fl <- lam$nn_floor %||% 0

  specs <- list(
    list(key = "full", spec = "TROP (full)",
         lambda = list(time = lam$time, unit = lam$unit, nn = lam$nn,
                       nn_floor = nn_fl)),
    list(key = "no_nn", spec = "No regression adjustment",
         lambda = list(time = lam$time, unit = lam$unit, nn = Inf,
                       nn_floor = nn_fl)),
    list(key = "no_unit", spec = "No unit weights",
         lambda = list(time = lam$time, unit = 0, nn = lam$nn,
                       nn_floor = nn_fl)),
    list(key = "no_time", spec = "No time weights",
         lambda = list(time = 0, unit = lam$unit, nn = lam$nn,
                       nn_floor = nn_fl)),
    list(key = "mc", spec = "Matrix completion",
         lambda = list(time = 0, unit = 0, nn = nn_finite,
                       nn_floor = nn_fl)),
    list(key = "did", spec = "Difference-in-differences",
         lambda = list(time = 0, unit = 0, nn = Inf,
                       nn_floor = nn_fl))
  )

  one <- function(s) {
    f <- if (identical(s$key, "full")) fit_full else
      trop(data, outcome, treatment, unit, time,
           covariates = covariates, lambda = s$lambda,
           anchor = anchor, se = se, control = control)
    data.frame(
      spec        = s$spec,
      lambda_time = s$lambda$time,
      lambda_unit = s$lambda$unit,
      lambda_nn   = s$lambda$nn,
      estimate    = f$estimate,
      std.error   = f$std.error %||% NA_real_,
      conf.low    = f$conf.low  %||% NA_real_,
      conf.high   = f$conf.high %||% NA_real_,
      rank        = f$rank      %||% NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, lapply(specs, one))
  rownames(out) <- NULL

  if (isTRUE(rmse)) {
    pm <- .panel_to_matrices(data, outcome, treatment, unit, time)
    Ym <- pm$Y; Wm <- pm$W
    controls <- which(rowSums(Wm) == 0)
    treated  <- which(rowSums(Wm) > 0)
    if (length(controls) < 2L)
      stop("rmse = TRUE needs at least two control units.", call. = FALSE)
    Tt <- ncol(Ym)
    if (is.null(horizon))
      horizon <- max(1L, sum(colSums(Wm[treated, , drop = FALSE]) > 0))
    horizon <- min(horizon, Tt - 1L)
    if (is.null(n_pseudo))
      n_pseudo <- max(1L, min(length(treated), length(controls) - 1L))
    n_pseudo <- min(n_pseudo, length(controls) - 1L)
    held <- utils::tail(seq_len(Tt), horizon)
    Yc <- Ym[controls, , drop = FALSE]; ncy <- nrow(Yc)
    # placebo-ATT RMSE (true effect 0): assign a pseudo block-treatment to control
    # units and score each penalty spec by sqrt(mean(estimate^2)). The SAME pseudo
    # draws are scored under every spec, so the ratios are paired -- the relative
    # RMSE with TROP (full) = 1 is the Table-5 quantity of the paper.
    if (!is.null(control$seed)) set.seed(control$seed)
    ps_list <- lapply(seq_len(n_runs), function(r) sample(seq_len(ncy), n_pseudo))
    lam_list <- lapply(specs, `[[`, "lambda")
    err <- matrix(NA_real_, nrow = n_runs, ncol = length(specs))
    for (r in seq_len(n_runs)) {
      Wp <- matrix(0, ncy, Tt, dimnames = dimnames(Yc))
      Wp[ps_list[[r]], held] <- 1
      patp <- .assignment_pattern(Wp)
      for (j in seq_along(lam_list))
        err[r, j] <- tryCatch(
          .trop_att(Yc, Wp, lam_list[[j]], control, "pooled", patp)$att,
          error = function(e) NA_real_)
    }
    rmse_vec <- sqrt(colMeans(err^2, na.rm = TRUE))
    full_idx <- which(vapply(specs, function(s) identical(s$key, "full"),
                             logical(1)))
    out$rmse     <- rmse_vec
    out$rmse_rel <- rmse_vec / rmse_vec[full_idx]
    attr(out, "rmse_runs")     <- n_runs
    attr(out, "rmse_horizon")  <- horizon
    attr(out, "rmse_n_pseudo") <- n_pseudo
  }
  attr(out, "outcome") <- outcome
  attr(out, "anchor")  <- anchor
  attr(out, "se")      <- fit_full$se.method
  attr(out, "conf.level") <- fit_full$conf.level %||% 0.95
  class(out) <- c("trop_ablation", "data.frame")
  out
}

# ---- display helpers --------------------------------------------------------

# Format a numeric vector for display; Inf becomes a format-specific token and
# NA becomes an empty cell.
.abl_num <- function(x, digits, inf = "Inf") {
  vapply(x, function(z) {
    if (is.na(z)) ""
    else if (is.infinite(z)) if (z > 0) inf else paste0("-", inf)
    else formatC(z, format = "f", digits = digits)
  }, character(1))
}

# Assemble the display columns shared by every output format. Standard-error and
# CI columns are included only when an SE was actually computed.
.abl_columns <- function(x, digits, inf, ci_fmt) {
  has_se <- any(is.finite(x$std.error))
  cols <- list(
    Specification = x$spec,
    lt = .abl_num(x$lambda_time, digits, inf),
    lu = .abl_num(x$lambda_unit, digits, inf),
    ln = .abl_num(x$lambda_nn,   digits, inf),
    ATT = .abl_num(x$estimate, digits, inf)
  )
  if (has_se) {
    cols$SE <- .abl_num(x$std.error, digits, inf)
    lo <- .abl_num(x$conf.low,  digits, inf)
    hi <- .abl_num(x$conf.high, digits, inf)
    cols[["95% CI"]] <- ifelse(lo == "" | hi == "", "", ci_fmt(lo, hi))
  }
  if (!is.null(x$rmse)) {
    cols$RMSE <- .abl_num(x$rmse, digits, inf)
    cols[["vs TROP"]] <- .abl_num(x$rmse_rel, digits, inf)
  }
  cols$`rank(L)` <- ifelse(is.na(x$rank), "", as.character(x$rank))
  cols
}

#' @rdname trop_ablation
#' @param x A `trop_ablation` object.
#' @param digits Number of decimal places for the penalties and estimates.
#' @param ... Unused.
#' @export
print.trop_ablation <- function(x, digits = 3, ...) {
  cols <- .abl_columns(x, digits, inf = "Inf",
                       ci_fmt = function(lo, hi) sprintf("[%s, %s]", lo, hi))
  disp <- as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE)
  names(disp)[names(disp) %in% c("lt", "lu", "ln")] <-
    c("lambda_time", "lambda_unit", "lambda_nn")
  cat(sprintf("TROP penalty ablation  (outcome: %s; anchor: %s; SE: %s)\n\n",
              attr(x, "outcome") %||% "?", attr(x, "anchor") %||% "?",
              attr(x, "se") %||% "none"))
  print(disp, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' Render a TROP ablation table as paste-ready LaTeX or Markdown
#'
#' Produces a publication-quality table from a `trop_ablation()` result. The
#' `"latex"` output uses \pkg{booktabs} rules and math-mode penalty headers,
#' matching the style of the paper's tables (load the booktabs LaTeX package).
#'
#' @param x A `trop_ablation` object from `trop_ablation()`.
#' @param output `"latex"` (booktabs) or `"markdown"` (GitHub pipe table).
#' @param digits Number of decimal places.
#' @param caption,label Table caption and cross-reference label (LaTeX only).
#' @param ... Unused.
#' @return A character vector of lines (paste-ready). Combine with
#'   `cat(..., sep = "\n")` or `writeLines()`.
#' @seealso `trop_ablation()`
#' @examples
#' df <- sim_panel(N = 14, T = 8, n_treated = 3, t0 = 5, att = 2, seed = 1)
#' ab <- trop_ablation(df, "y", "w", "id", "t",
#'                     lambda_full = list(time = 0.1, unit = 0.5, nn = 2))
#' writeLines(format(ab, "markdown"))
#' @export
format.trop_ablation <- function(x, output = c("latex", "markdown"),
                                 digits = 3,
                                 caption = "TROP penalty ablation",
                                 label = "tab:trop-ablation", ...) {
  output <- match.arg(output)
  conf <- round(100 * (attr(x, "conf.level") %||% 0.95))

  if (output == "markdown") {
    cols <- .abl_columns(x, digits, inf = "Inf",
                         ci_fmt = function(lo, hi) sprintf("[%s, %s]", lo, hi))
    nm <- names(cols)
    head_map <- c(Specification = "Specification",
                  lt = "&lambda;_time", lu = "&lambda;_unit", ln = "&lambda;_nn",
                  ATT = "ATT", SE = "SE",
                  "95% CI" = paste0(round(100 * (attr(x, "conf.level") %||% 0.95)),
                                    "% CI"),
                  RMSE = "RMSE", "vs TROP" = "RMSE/TROP",
                  "rank(L)" = "rank(L)")
    hdr <- ifelse(nm %in% names(head_map), head_map[nm], nm)
    align <- ifelse(nm == "Specification", ":--", "--:")
    body <- do.call(paste, c(lapply(cols, identity), sep = " | "))
    c(paste0("| ", paste(hdr, collapse = " | "), " |"),
      paste0("| ", paste(align, collapse = " | "), " |"),
      paste0("| ", body, " |"))
  } else {
    inf <- "$\\infty$"
    cols <- .abl_columns(x, digits, inf = inf,
                         ci_fmt = function(lo, hi) sprintf("[%s, %s]", lo, hi))
    nm <- names(cols)
    head_map <- c(
      "Specification" = "Specification",
      "lt" = "$\\lambda_{\\mathrm{time}}$",
      "lu" = "$\\lambda_{\\mathrm{unit}}$",
      "ln" = "$\\lambda_{\\mathrm{nn}}$",
      "ATT" = "ATT",
      "SE" = "SE",
      "95% CI" = paste0(conf, "\\% CI"),
      "RMSE" = "RMSE",
      "vs TROP" = "RMSE$/$TROP",
      "rank(L)" = "rank($L$)")
    headers <- ifelse(nm %in% names(head_map), head_map[nm], nm)
    # left-align the label column, right-align the rest
    coldef <- paste0("l", paste(rep("r", length(nm) - 1L), collapse = ""))
    esc <- function(v) gsub("%", "\\\\%", v)               # escape stray %
    rows <- do.call(paste, c(lapply(cols, esc), sep = " & "))
    rows <- paste0(rows, " \\\\")
    c("\\begin{table}[!ht]",
      "\\centering",
      sprintf("\\caption{%s}", caption),
      sprintf("\\label{%s}", label),
      sprintf("\\begin{tabular}{%s}", coldef),
      "\\toprule",
      paste0(paste(headers, collapse = " & "), " \\\\"),
      "\\midrule",
      rows,
      "\\bottomrule",
      "\\end{tabular}",
      "\\end{table}")
  }
}

# ---- figure output ----------------------------------------------------------

# Row label: human name + the penalty constraint that defines the specification
# (so "DID-like"/"MC-like" are explained the way Table 5 of the paper does it).
# Returns a plotmath expression so lambda renders as a proper Greek symbol.
#' @keywords internal
#' @noRd
.abl_label_expr <- function(spec) {
  switch(spec,
    "TROP (full)"               = "TROP (full)",
    "No regression adjustment"  =
      bquote("No regression adj.  (" * lambda[nn] * " = " * infinity * ")"),
    "No unit weights"           =
      bquote("No unit weights  (" * lambda[unit] * " = 0)"),
    "No time weights"           =
      bquote("No time weights  (" * lambda[time] * " = 0)"),
    "Matrix completion"         =
      bquote("MC-like  (" * lambda[unit] * " = " * lambda[time] * " = 0)"),
    "Difference-in-differences" =
      bquote("DID-like  (" * lambda[unit] * " = " * lambda[time] *
             " = 0, " * lambda[nn] * " = " * infinity * ")"),
    spec)
}

# ASCII stand-in used only to size the label column.
#' @keywords internal
#' @noRd
.abl_label_plain <- function(spec) {
  switch(spec,
    "No regression adjustment"  = "No regression adj.  (l_nn = Inf)",
    "No unit weights"           = "No unit weights  (l_unit = 0)",
    "No time weights"           = "No time weights  (l_time = 0)",
    "Matrix completion"         = "MC-like  (l_unit = l_time = 0)",
    "Difference-in-differences" = "DID-like  (l_unit = l_time = 0, l_nn = Inf)",
    spec)
}

# Shared column spec used by both the renderer and the auto-sizer: a left label
# column plus the value columns (RMSE when available, otherwise the ATT block).
#' @keywords internal
#' @noRd
.abl_fig_cols <- function(x, digits = 3) {
  num <- function(v) .abl_num(v, digits, inf = "Inf")
  vcols <- list()
  if (!is.null(x$rmse)) {
    vcols[["RMSE"]]    <- num(x$rmse)
    vcols[["vs TROP"]] <- num(x$rmse_rel)
  } else {
    vcols[["ATT"]] <- num(x$estimate)
    if (any(is.finite(x$std.error))) {
      vcols[["SE"]] <- num(x$std.error)
      lo <- num(x$conf.low); hi <- num(x$conf.high)
      vcols[["95% CI"]] <- ifelse(lo == "" | hi == "", "",
                                  sprintf("[%s, %s]", lo, hi))
    }
  }
  list(
    labels_expr  = lapply(x$spec, .abl_label_expr),
    labels_plain = vapply(x$spec, .abl_label_plain, character(1)),
    vnames = names(vcols),
    vcols  = vcols
  )
}

# Vertical layout in inches (top-down), shared by the renderer and the
# auto-sizer so the row pitch and the title gaps stay at a normal, fixed level
# regardless of the device height (rather than stretching to fill a tall pane).
#' @keywords internal
#' @noRd
.abl_layout <- function(nb) {
  c_title <- 0.13
  c_sub   <- c_title + 0.24
  y_top   <- c_sub   + 0.18
  c_head  <- y_top   + 0.15
  y_mid   <- c_head  + 0.15
  ROW     <- 0.21                        # single-line row pitch (normal table)
  c_rows  <- (y_mid + 0.16) + (seq_len(nb) - 1L) * ROW
  y_bot   <- c_rows[nb] + 0.13
  list(c_title = c_title, c_sub = c_sub, y_top = y_top, c_head = c_head,
       y_mid = y_mid, c_rows = c_rows, y_bot = y_bot, content_h = y_bot + 0.07)
}

# Pure-grid renderer (no gridExtra/gt/Chrome, so it works headless). Booktabs
# style after the paper: serif type, centred title/subtitle, no shading, three
# rules. Draws onto the current device; .render_ablation() opens a file device.
#' @keywords internal
#' @noRd
.draw_trop_ablation <- function(x, digits = 3) {
  fc <- .abl_fig_cols(x, digits)
  labels_expr <- fc$labels_expr; vnames <- fc$vnames; vcols <- fc$vcols
  nb <- nrow(x); nvc <- length(vcols)
  has_rmse <- !is.null(x$rmse)

  label_w <- max(nchar(fc$labels_plain), nchar("Specification"))
  vw <- vapply(seq_len(nvc), function(j)
    max(nchar(vnames[j]), max(nchar(vcols[[j]]), 1L)), integer(1))
  wrel  <- c(label_w, vw) + 3L
  wfrac <- wrel / sum(wrel)
  xr <- cumsum(wfrac); xl <- xr - wfrac

  title <- "TROP penalty ablation"
  meta <- if (has_rmse)
    sprintf("Placebo RMSE over %s runs (TROP = 1.00)",
            attr(x, "rmse_runs") %||% NA)
  else
    sprintf("SE: %s", attr(x, "se") %||% "none")

  c_rule <- "#222222"; c_text <- "#111111"; c_sub <- "#555555"; ff <- "serif"
  lay <- .abl_layout(nb)

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = 0.92, height = 1))
  on.exit(grid::popViewport(), add = TRUE)

  # Centre the fixed-height content block vertically; map inch offsets from its
  # top, so the row pitch is constant whatever the device/pane height is.
  dev_h <- grid::convertHeight(grid::unit(1, "npc"), "inches", valueOnly = TRUE)
  ytop  <- 0.5 + (lay$content_h / 2) / max(dev_h, lay$content_h)
  at <- function(off) grid::unit(ytop, "npc") - grid::unit(off, "inches")

  xpos <- function(j) if (j == 1L) grid::unit(xl[1], "npc") else
    grid::unit(xr[j], "npc") - grid::unit(2, "pt")
  jjust <- function(j) if (j == 1L) c("left", "centre") else c("right", "centre")
  put <- function(lbl, j, off, bold = FALSE)
    grid::grid.text(lbl, x = xpos(j), y = at(off), just = jjust(j),
                    gp = grid::gpar(fontsize = 11, fontfamily = ff, col = c_text,
                                    fontface = if (bold) "bold" else "plain"))
  rule <- function(off, lwd) grid::grid.lines(
    x = grid::unit(c(0, 1), "npc"),
    y = grid::unit(c(ytop, ytop), "npc") - grid::unit(c(off, off), "inches"),
    gp = grid::gpar(col = c_rule, lwd = lwd))

  grid::grid.text(title, x = 0.5, y = at(lay$c_title), just = c("centre", "centre"),
                  gp = grid::gpar(fontface = "bold", fontsize = 15,
                                  fontfamily = ff, col = c_text))
  grid::grid.text(meta, x = 0.5, y = at(lay$c_sub), just = c("centre", "centre"),
                  gp = grid::gpar(fontsize = 10.5, fontfamily = ff, col = c_sub))

  rule(lay$y_top, 1.6)                                  # toprule
  put("Specification", 1L, lay$c_head, bold = TRUE)
  for (j in seq_len(nvc)) put(vnames[j], j + 1L, lay$c_head, bold = TRUE)
  rule(lay$y_mid, 1.0)                                  # midrule
  for (i in seq_len(nb)) {
    put(labels_expr[[i]], 1L, lay$c_rows[i])
    for (j in seq_len(nvc)) put(vcols[[j]][i], j + 1L, lay$c_rows[i])
  }
  rule(lay$y_bot, 1.6)                                  # bottomrule

  invisible()
}

# Open a file device by extension (png default, pdf for vector) and draw.
#' @keywords internal
#' @noRd
.render_ablation <- function(x, file = NULL, width = NULL, height = NULL,
                             res = 200, digits = 3) {
  fc <- .abl_fig_cols(x, digits)
  nb <- nrow(x)
  label_chars <- max(nchar(fc$labels_plain), nchar("Specification"))
  val_chars <- sum(vapply(seq_along(fc$vcols), function(j)
    max(nchar(fc$vnames[j]), max(nchar(fc$vcols[[j]]), 1L)) + 3L, integer(1)))
  tot_chars <- label_chars + 3L + val_chars
  if (is.null(width))  width  <- max(6.5, 0.092 * tot_chars + 1)
  if (is.null(height)) height <- .abl_layout(nb)$content_h + 0.22

  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "") { file <- paste0(file, ".png"); ext <- "png" }
    if (ext == "pdf") {
      grDevices::pdf(file, width = width, height = height)
    } else if (ext == "png") {
      grDevices::png(file, width = width, height = height, units = "in", res = res)
    } else {
      stop("plot() on a trop_ablation writes .png or .pdf; ",
           "for .tex/.md use format(x, \"latex\"/\"markdown\").", call. = FALSE)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  .draw_trop_ablation(x, digits = digits)
  invisible(file)
}

#' @rdname trop_ablation
#' @param file Optional output path. When supplied, the table is rendered to that
#'   file instead of the active graphics device; the format is taken from the
#'   extension (`.png`, the default, or `.pdf` for vector output). For LaTeX or
#'   Markdown source use [format()] instead.
#' @param width,height Figure size in inches. Defaults adapt to the number of
#'   rows and columns.
#' @param res Resolution in PPI for the `.png` device (ignored for `.pdf`).
#' @export
plot.trop_ablation <- function(x, file = NULL, width = NULL, height = NULL,
                               res = 200, digits = 3, ...) {
  .render_ablation(x, file = file, width = width, height = height,
                   res = res, digits = digits)
  invisible(x)
}
