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
