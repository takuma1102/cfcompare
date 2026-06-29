# Penalty ablation for TROP: refit the estimator under a sequence of penalty
# constraints (the "Table 5"-style robustness check of the paper -- constraining
# lambda_time and/or lambda_unit to zero, and/or lambda_nn to infinity) and
# return the results as a publication-ready table rather than a plot.

#' Penalty ablation for the TROP estimator
#'
#' Refits TROP under a sequence of penalty constraints to show how the ATT moves
#' as the estimator is stripped back towards matrix completion and
#' difference-in-differences. This is the robustness exercise behind Table 5 of
#' Athey, Imbens, Qu & Viviano (2025): constraining `lambda_unit` and/or
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
#' @return An object of class `trop_ablation` (a `data.frame`) with one row per
#'   specification: `spec`, the three penalties `lambda_time`/`lambda_unit`/
#'   `lambda_nn`, `estimate`, `std.error`, `conf.low`/`conf.high`, and `rank`.
#' @references Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2025).
#'   Triply Robust Panel Estimators. arXiv:2508.21536.
#' @seealso [trop()]; `format()` for paste-ready LaTeX/Markdown output
#' @examples
#' df <- sim_panel(N = 30, T = 10, n_treated = 4, t0 = 5, att = 2, seed = 1)
#' ab <- trop_ablation(df, "y", "w", "id", "t",
#'                     control = trop_control(n_cv_cells = 10L, cv_cycles = 1L))
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
                          lambda_full = NULL) {
  fit_full <- trop(data, outcome, treatment, unit, time,
                   covariates = covariates, lambda = lambda_full,
                   anchor = anchor, se = se, control = control)
  lam <- fit_full$lambda
  nn_finite <- if (is.finite(lam$nn)) lam$nn else
    stats::sd(data[[outcome]], na.rm = TRUE)

  specs <- list(
    list(key = "full", spec = "TROP (full)",
         lambda = list(time = lam$time, unit = lam$unit, nn = lam$nn)),
    list(key = "no_nn", spec = "No regression adjustment",
         lambda = list(time = lam$time, unit = lam$unit, nn = Inf)),
    list(key = "no_unit", spec = "No unit weights",
         lambda = list(time = lam$time, unit = 0, nn = lam$nn)),
    list(key = "no_time", spec = "No time weights",
         lambda = list(time = 0, unit = lam$unit, nn = lam$nn)),
    list(key = "mc", spec = "Matrix completion",
         lambda = list(time = 0, unit = 0, nn = nn_finite)),
    list(key = "did", spec = "Difference-in-differences",
         lambda = list(time = 0, unit = 0, nn = Inf))
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
#' df <- sim_panel(N = 30, T = 10, n_treated = 4, t0 = 5, att = 2, seed = 1)
#' ab <- trop_ablation(df, "y", "w", "id", "t",
#'                     control = trop_control(n_cv_cells = 10L, cv_cycles = 1L))
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

# Pure-grid renderer for the ablation table (no gridExtra/gt/Chrome needed, so
# it works on any R install and on a headless machine). Draws onto the current
# device/page; .render_ablation() opens a file device around it when asked.
#' @keywords internal
#' @noRd
.draw_trop_ablation <- function(x, digits = 3) {
  cols <- .abl_columns(x, digits, inf = "Inf",
                       ci_fmt = function(lo, hi) sprintf("[%s, %s]", lo, hi))
  key  <- names(cols)
  hdr  <- key
  hdr[key == "lt"] <- "lambda_time"
  hdr[key == "lu"] <- "lambda_unit"
  hdr[key == "ln"] <- "lambda_nn"
  align <- ifelse(key == "Specification", "l", "r")
  body  <- lapply(cols, as.character)
  nc <- length(body); nb <- length(body[[1]])

  charw <- vapply(seq_len(nc), function(j)
    max(nchar(hdr[j]), max(nchar(body[[j]]), 1L)), integer(1)) + 2L
  wfrac <- charw / sum(charw)
  xr <- cumsum(wfrac); xl <- xr - wfrac

  title <- "TROP penalty ablation"
  subtitle <- sprintf("outcome: %s    anchor: %s    SE: %s",
                      attr(x, "outcome") %||% "?", attr(x, "anchor") %||% "?",
                      attr(x, "se") %||% "none")

  c_head_bg <- "#2b3a55"; c_head_fg <- "white"
  c_stripe  <- "#eef1f6"; c_rule <- "#9aa6bd"; c_text <- "#1d2740"

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = 0.94, height = 0.92))
  on.exit(grid::popViewport(), add = TRUE)

  top_title <- 0.17
  table_top <- 1 - top_title
  n_tot <- nb + 1L
  rh <- table_top / n_tot

  grid::grid.text(title, x = 0, y = 0.98, just = c("left", "top"),
                  gp = grid::gpar(fontface = "bold", fontsize = 14, col = c_text))
  grid::grid.text(subtitle, x = 0, y = 0.98 - 0.075, just = c("left", "top"),
                  gp = grid::gpar(fontsize = 9.5, col = "#5a6477"))

  cell_x <- function(j) if (align[j] == "l")
    grid::unit(xl[j], "npc") + grid::unit(4, "pt") else
    grid::unit(xr[j], "npc") - grid::unit(4, "pt")
  cell_just <- function(j) c(if (align[j] == "l") "left" else "right", "centre")

  hy <- table_top - rh / 2
  grid::grid.rect(x = 0.5, y = hy, width = 1, height = rh,
                  gp = grid::gpar(fill = c_head_bg, col = NA))
  for (j in seq_len(nc))
    grid::grid.text(hdr[j], x = cell_x(j), y = hy, just = cell_just(j),
                    gp = grid::gpar(col = c_head_fg, fontface = "bold", fontsize = 10))

  for (i in seq_len(nb)) {
    ry <- table_top - rh * (i + 0.5)
    if (i %% 2L == 0L)
      grid::grid.rect(x = 0.5, y = ry, width = 1, height = rh,
                      gp = grid::gpar(fill = c_stripe, col = NA))
    for (j in seq_len(nc))
      grid::grid.text(body[[j]][i], x = cell_x(j), y = ry, just = cell_just(j),
                      gp = grid::gpar(fontsize = 9.5, col = c_text,
                                      fontfamily = if (align[j] == "r") "mono" else ""))
  }

  for (yy in c(table_top, table_top - rh, table_top - rh * n_tot))
    grid::grid.lines(x = c(0, 1), y = yy, gp = grid::gpar(col = c_rule, lwd = 1.2))

  invisible()
}

# Open a file device by extension (png default, pdf for vector) and draw.
#' @keywords internal
#' @noRd
.render_ablation <- function(x, file = NULL, width = NULL, height = NULL,
                             res = 200, digits = 3) {
  cols <- .abl_columns(x, digits, inf = "Inf",
                       ci_fmt = function(lo, hi) sprintf("[%s, %s]", lo, hi))
  nb <- nrow(x)
  tot_chars <- sum(vapply(cols, function(z) max(nchar(z), 6L), integer(1))) +
    2L * length(cols)
  if (is.null(width))  width  <- max(6, 0.10 * tot_chars + 1)
  if (is.null(height)) height <- 0.95 + 0.34 * (nb + 1)

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
