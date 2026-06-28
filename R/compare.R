# Comparison harness: run several panel estimators on the SAME data and return
# them on one tidy schema (class 'cf_att_tbl').

#' Compare DID, SDID, MC and TROP on the same panel
#'
#' Runs a set of panel estimators on a single long panel and returns their
#' average-treatment-effect-on-the-treated (ATT) estimates on a common tidy
#' schema, so applied researchers can compare them at a glance. DID, MC and TROP
#' are computed natively (no external dependencies); SDID and SC are routed
#' through the \pkg{synthdid} package when available, and an alternative MC/IFE
#' can be routed through \pkg{gsynth}. Methods whose optional package is missing,
#' or that do not apply to the design, are skipped with a message.
#'
#' @inheritParams trop
#' @param methods Character vector of methods to run. Any of `"DID"`, `"SDID"`,
#'   `"SC"`, `"MC"`, `"DIFP"`, `"TROP"`, `"gsynth"`, `"augsynth"`, `"CS"`.
#'   Defaults to the native engines plus SDID: `c("DID", "SDID", "MC", "TROP",
#'   "DIFP")`.
#' @param exclude Optional character vector of methods to drop from `methods`
#'   (after defaults are applied). Convenient for running "everything except one",
#'   e.g. `exclude = "DIFP"`. Unknown names are ignored with a warning.
#' @param anchor Weight anchoring for TROP; see [trop()].
#' @param se Standard-error method for the native engines; see [trop()].
#' @return An object of class `cf_comparison`: a list with the tidy table
#'   `att` (class `cf_att_tbl`), per-method counterfactual matrices
#'   `counterfactual`, the native fit objects `fits`, and the reshaped `panel`.
#' @seealso [trop()], [autoplot.cf_comparison()], [plot_counterfactual()]
#' @examples
#' df <- sim_panel(N = 25, T = 14, n_treated = 4, t0 = 10, seed = 3)
#' cmp <- panel_compare(df, "y", "w", "id", "t",
#'                      methods = c("DID", "MC", "TROP"), se = "none",
#'                      control = trop_control(n_cv_cells = 8L, cv_cycles = 1L))
#' cmp$att
#' @export
panel_compare <- function(data, outcome, treatment, unit, time,
                          methods = c("DID", "SDID", "MC", "TROP", "DIFP"),
                          exclude = NULL,
                          anchor = "auto",
                          se = c("auto", "jackknife", "bootstrap", "placebo", "none"),
                          control = trop_control(),
                          verbose = FALSE) {
  se <- match.arg(se)
  known <- c("DID", "SDID", "SC", "MC", "DIFP", "TROP", "gsynth", "augsynth", "CS")
  methods <- .resolve_methods(methods, exclude, known)

  m <- .panel_to_matrices(data, outcome, treatment, unit, time)
  Y <- m$Y; W <- m$W
  pat <- .assignment_pattern(W)
  cl <- control$conf_level

  rows <- list()
  cfs <- list()
  fits <- list()

  add_native <- function(name, engine_name, eng) {
    rows[[name]] <<- .engine_row(name, engine_name, eng, pat, outcome)
    cfs[[name]] <<- eng$counterfactual
    fits[[name]] <<- eng
  }

  for (meth in methods) {
    if (verbose) message("Running ", meth, " ...")
    if (meth == "DID") {
      add_native("DID", "cfcompare", .engine_did(Y, W, pat, control, se, cl))
    } else if (meth == "MC") {
      add_native("MC", "cfcompare", .engine_mc(Y, W, pat, control, se, cl, verbose))
    } else if (meth == "DIFP") {
      eng <- .engine_difp(Y, W, pat, control, cl, se)
      if (inherits(eng, "skip")) {
        message("  skipping DIFP: ", eng$note)
        rows[["DIFP"]] <- .engine_row("DIFP", "cfcompare",
          list(estimate = NA_real_), pat, outcome, note = eng$note)
      } else {
        add_native("DIFP", "cfcompare", eng)
      }
    } else if (meth == "TROP") {
      add_native("TROP", "cfcompare",
                 .engine_trop(Y, W, pat, control, anchor, se, cl, verbose))
    } else if (meth %in% c("SDID", "SC")) {
      eng <- .engine_synthdid(Y, W, pat,
                              which = if (meth == "SDID") "sdid" else "sc", cl,
                              se = se)
      if (inherits(eng, "skip")) {
        message("  skipping ", meth, ": ", eng$note)
        rows[[meth]] <- .engine_row(meth, "synthdid",
          list(estimate = NA_real_), pat, outcome, note = eng$note)
      } else {
        rows[[meth]] <- .engine_row(meth, "synthdid", eng, pat, outcome)
        fits[[meth]] <- eng$fit
        if (!is.null(eng$counterfactual)) cfs[[meth]] <- eng$counterfactual
      }
    } else if (meth == "gsynth") {
      eng <- .engine_gsynth(data, outcome, treatment, unit, time, "mc", cl)
      if (inherits(eng, "skip")) {
        message("  skipping gsynth: ", eng$note)
        rows[["gsynth"]] <- .engine_row("gsynth", "gsynth",
          list(estimate = NA_real_), pat, outcome, note = eng$note)
      } else {
        rows[["gsynth"]] <- .engine_row("gsynth", "gsynth", eng, pat, outcome)
        fits[["gsynth"]] <- eng$fit
        if (!is.null(eng$counterfactual)) cfs[["gsynth"]] <- eng$counterfactual
      }
    } else if (meth == "augsynth") {
      eng <- .engine_augsynth(data, outcome, treatment, unit, time, pat, cl)
      if (inherits(eng, "skip")) {
        message("  skipping augsynth: ", eng$note)
        rows[["augsynth"]] <- .engine_row("augsynth", "augsynth",
          list(estimate = NA_real_), pat, outcome, note = eng$note)
      } else {
        rows[["augsynth"]] <- .engine_row("augsynth", "augsynth", eng, pat, outcome)
        fits[["augsynth"]] <- eng$fit
      }
    } else if (meth == "CS") {
      eng <- .engine_did_cs(data, outcome, treatment, unit, time, pat, cl)
      if (inherits(eng, "skip")) {
        message("  skipping CS: ", eng$note)
        rows[["CS"]] <- .engine_row("CS", "did",
          list(estimate = NA_real_), pat, outcome, note = eng$note)
      } else {
        rows[["CS"]] <- .engine_row("CS", "did", eng, pat, outcome)
        fits[["CS"]] <- eng$fit
      }
    }
  }

  att <- do.call(rbind, rows[methods[methods %in% names(rows)]])
  rownames(att) <- NULL
  class(att) <- c("cf_att_tbl", "data.frame")

  structure(
    list(att = att, counterfactual = cfs, fits = fits, panel = m,
         pattern = pat, outcome = outcome, call = match.call()),
    class = "cf_comparison"
  )
}

#' @export
print.cf_comparison <- function(x, ...) {
  cat("Panel estimator comparison (ATT)\n")
  cat(sprintf("  design: %s | treated cells: %d across %d unit(s) | outcome: %s\n",
              x$pattern$type, x$pattern$n_treated_cells,
              x$pattern$n_treated_units, x$outcome))
  cat("\n")
  print(x$att, ...)
  invisible(x)
}

#' @export
print.cf_att_tbl <- function(x, digits = 4, ...) {
  df <- as.data.frame(x)
  num <- c("estimate", "std.error", "conf.low", "conf.high")
  for (nm in num) if (nm %in% names(df)) df[[nm]] <- round(df[[nm]], digits)
  print.data.frame(df, row.names = FALSE)
  invisible(x)
}

# ---- tidiers ----------------------------------------------------------------

# Canonical column order of the cf_att_tbl schema.
.att_schema <- c("method", "estimate", "std.error", "conf.low", "conf.high",
                 "n_treated_cells", "n_treated_units", "outcome", "engine",
                 "rank", "note")

# Row-bind a list of ATT data frames into one `cf_att_tbl`, aligning columns
# (any missing column is filled with NA) and keeping the canonical schema order
# for known columns. NULL / empty parts are dropped. This is the shared backbone
# of `as_att.list()` and `bind_att()`, so stacking tables is robust to inputs
# that carry extra or missing columns.
.rbind_att <- function(parts) {
  parts <- lapply(parts, function(p) if (is.null(p)) NULL else as.data.frame(p))
  parts <- Filter(function(p) !is.null(p) && nrow(p) > 0L, parts)
  if (!length(parts)) return(NULL)
  all_cols <- unique(unlist(lapply(parts, names), use.names = FALSE))
  ordered  <- c(intersect(.att_schema, all_cols), setdiff(all_cols, .att_schema))
  parts <- lapply(parts, function(d) {
    for (m in setdiff(ordered, names(d))) d[[m]] <- NA
    d[ordered]
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  class(out) <- c("cf_att_tbl", "data.frame")
  out
}

#' Coerce estimator output to the cfcompare ATT schema
#'
#' Brings the output of [trop()], [panel_compare()], a \pkg{synthdid} estimate,
#' or any object with an `estimate` into the common `cf_att_tbl` schema, so
#' results computed elsewhere can be slotted into the same comparison and plots.
#'
#' @param x An object to tidy.
#' @param ... Passed to methods (e.g. `method`, `outcome`, `conf_level`).
#' @return A `cf_att_tbl` (a `data.frame`).
#' @export
as_att <- function(x, ...) UseMethod("as_att")

#' @export
as_att.cf_comparison <- function(x, ...) x$att

#' @export
as_att.cf_att_tbl <- function(x, ...) x

#' @export
as_att.trop <- function(x, method = "TROP", ...) {
  row <- data.frame(
    method = method, estimate = x$estimate, std.error = x$std.error,
    conf.low = x$conf.low, conf.high = x$conf.high,
    n_treated_cells = x$pattern$n_treated_cells,
    n_treated_units = x$pattern$n_treated_units,
    outcome = x$outcome, engine = "cfcompare", rank = x$rank,
    note = NA_character_, stringsAsFactors = FALSE
  )
  class(row) <- c("cf_att_tbl", "data.frame")
  row
}

#' @export
as_att.synthdid_estimate <- function(x, method = "SDID",
                                           outcome = NA_character_,
                                           conf_level = 0.95, ...) {
  att <- as.numeric(x)
  v <- tryCatch(stats::vcov(x, method = "jackknife"), error = function(e) NA_real_)
  se <- if (is.matrix(v)) sqrt(v[1, 1]) else as.numeric(v)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  row <- data.frame(
    method = method, estimate = att, std.error = se,
    conf.low = att - z * se, conf.high = att + z * se,
    n_treated_cells = NA_integer_, n_treated_units = NA_integer_,
    outcome = outcome, engine = "synthdid", rank = NA_integer_,
    note = NA_character_, stringsAsFactors = FALSE
  )
  class(row) <- c("cf_att_tbl", "data.frame")
  row
}

#' @export
as_att.augsynth <- function(x, method = "augsynth",
                                  outcome = NA_character_,
                                  conf_level = 0.95, ...) {
  s <- summary(x, inf_type = "jackknife")
  avg <- s$average_att
  att <- as.numeric(avg$Estimate[1])
  se <- tryCatch(as.numeric(avg$Std.Error[1]), error = function(e) NA_real_)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  row <- data.frame(
    method = method, estimate = att, std.error = se,
    conf.low = att - z * se, conf.high = att + z * se,
    n_treated_cells = NA_integer_, n_treated_units = NA_integer_,
    outcome = outcome, engine = "augsynth", rank = NA_integer_,
    note = NA_character_, stringsAsFactors = FALSE
  )
  class(row) <- c("cf_att_tbl", "data.frame")
  row
}

#' @export
as_att.list <- function(x, ...) {
  parts <- lapply(x, as_att, ...)
  nms <- names(x)
  if (!is.null(nms)) {
    # use list names as method labels, but only for single-row parts (a name
    # cannot sensibly relabel a multi-row table such as a panel_compare result).
    parts <- Map(function(p, nm) {
      if (!is.na(nm) && nzchar(nm) && nrow(p) == 1L) p$method <- nm
      p
    }, parts, nms)
  }
  .rbind_att(parts)
}

#' @export
as_att.default <- function(x, method = NA_character_,
                                 outcome = NA_character_, ...) {
  est <- tryCatch(as.numeric(x[["estimate"]]), error = function(e) NULL)
  if (is.null(est)) {
    est <- tryCatch(as.numeric(x), error = function(e) NA_real_)
  }
  row <- data.frame(
    method = method, estimate = est[1], std.error = NA_real_,
    conf.low = NA_real_, conf.high = NA_real_,
    n_treated_cells = NA_integer_, n_treated_units = NA_integer_,
    outcome = outcome, engine = "external", rank = NA_integer_,
    note = NA_character_, stringsAsFactors = FALSE
  )
  class(row) <- c("cf_att_tbl", "data.frame")
  row
}

#' Stack ATT results into one comparison table
#'
#' Row-binds several ATT results --- [trop()] fits, [panel_compare()] results,
#' \pkg{synthdid} estimates, existing `cf_att_tbl`s, or anything [as_att()]
#' understands --- into a single `cf_att_tbl` ready for [autoplot()] or
#' [plot_counterfactual()]. Each argument is passed through [as_att()] first, so
#' multi-row inputs (e.g. a `panel_compare()` result) keep all their rows, and
#' columns are aligned automatically. This replaces the manual
#' `rbind(as.data.frame(...), ...)` + `class(...) <- "cf_att_tbl"` dance.
#'
#' Names given to the arguments become `method` labels, overriding the label
#' already on a *single-row* result; a name on a multi-row input is ignored with
#' a warning (it cannot sensibly relabel several methods at once). For per-object
#' control of other fields such as `outcome`, coerce that object with [as_att()]
#' first and pass the result in unnamed.
#'
#' @param ... Objects coercible by [as_att()], optionally named to set the
#'   `method` label.
#' @return A `cf_att_tbl` (a `data.frame`) with the inputs stacked row-wise.
#' @seealso [as_att()], [panel_compare()], [compare_se_modes()]
#' @examples
#' df <- sim_panel(N = 25, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 1)
#' ctrl <- trop_control(n_cv_cells = 8L, cv_cycles = 1L)
#' f_pooled   <- trop(df, "y", "w", "id", "t", anchor = "pooled",
#'                    se = "none", control = ctrl)
#' f_per_cell <- trop(df, "y", "w", "id", "t", anchor = "per_cell",
#'                    se = "none", control = ctrl)
#' bind_att(pooled = f_pooled, per_cell = f_per_cell)
#' @export
bind_att <- function(...) {
  args <- list(...)
  if (!length(args)) stop("bind_att() needs at least one object.", call. = FALSE)
  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  parts <- vector("list", length(args))
  for (i in seq_along(args)) {
    a <- as_att(args[[i]])
    if (nzchar(nms[i])) {
      if (nrow(a) == 1L) {
        a$method <- nms[i]
      } else {
        warning("bind_att(): name \"", nms[i],
                "\" ignored for a multi-row result (", nrow(a), " rows).",
                call. = FALSE)
      }
    }
    parts[[i]] <- a
  }
  .rbind_att(parts)
}
