# Out-of-sample placebo validation across estimators, written independently
# from the paper. The idea follows the paper's simulation logic: hold out a
# placebo block of control cells (or assign a placebo treatment to control
# units), have each method predict / estimate it while the true effect is
# zero, and score by root-mean-square error. Lower is better. No code from any
# official TROP implementation is used.

#' Honest out-of-sample prediction for a native (unified-program) method.
#'
#' Fits the working model on the cells flagged available in `fit_mask`
#' (real treated cells and the held-out placebo cells are both excluded), then
#' returns the predicted untreated outcome matrix.
#' @keywords internal
#' @noRd
.oos_predict_native <- function(Y, fit_mask, wmat, lambda_nn, control) {
  fit <- .mcnnm_fit(Y, fit_mask, wmat, lambda_nn,
                    max_iter = control$max_iter, tol = control$tol,
                    svd_method = control$svd %||% "truncated")
  fit$M
}

#' Build a placebo weight matrix anchored to pseudo-treated units / held periods.
#' @keywords internal
#' @noRd
.oos_weight_matrix <- function(Y, excl, pseudo_units, held_periods, lam) {
  # excl: N x T 0/1, 1 = cell excluded from fitting (real treated or held out).
  # Unit distances use only available control cells (excl == 0).
  du <- .unit_distance_pooled(Y, excl, pseudo_units)
  .trop_weight_matrix(du, held_periods, ncol(Y), lam)
}

#' Out-of-sample RMSE across panel estimators
#'
#' Compares estimators by how well each does on held-out control data, following
#' the "random blocks" placebo idea of the doubly/triply robust panel estimator
#' paper. Two scoring rules are available:
#'
#' * `metric = "placebo"` (default): in each run a random set of control units is
#'   given a *placebo* block treatment in the final `horizon` periods (true effect
#'   zero), every method is estimated on that control-only panel, and the score is
#'   `sqrt(mean(ATT^2))` across runs. This works for **every** method -- native
#'   (DID, MC, TROP) and wrapped (`SDID`/`SC` via \pkg{synthdid}, `gsynth`,
#'   `augsynth`, `CS` = Callaway & Sant'Anna via \pkg{did}).
#' * `metric = "prediction"`: per-cell one-step-ahead held-out RMSE (the paper's
#'   Table-1 unit-period counterfactual prediction error). A block of control
#'   cells is held out and each method predicts the untreated outcome there via
#'   its counterfactual matrix ([counterfactual_matrix()]). Supported for the
#'   native methods (`DID`, `MC`, `TROP`), for `DIFP`, and -- when \pkg{synthdid}
#'   is installed -- for `SDID` and `SC`. The remaining wrapped methods
#'   (`gsynth`, `augsynth`, `CS`) are reported as `NA` for this metric; use
#'   `metric = "placebo"` for them.
#'
#' Lower is better. Native methods are tuned once on the real data (DID is
#' parameter-free; MC selects `lambda_nn`; TROP selects the full triplet by
#' cross-validation) and reused across runs. Wrapped methods are skipped with a
#' note when their package is missing or the design does not apply.
#'
#' @param data A long `data.frame`, one row per unit-time.
#' @param outcome,treatment,unit,time Column names (strings).
#' @param methods Methods to compare; subset of `"DID"`, `"MC"`, `"TROP"`,
#'   `"DIFP"`, `"SDID"`, `"SC"`, `"gsynth"`, `"augsynth"`, `"CS"`. Defaults to
#'   `c("DID", "SC", "SDID", "MC", "TROP", "DIFP")`.
#' @param exclude Optional character vector of methods to drop from `methods`
#'   (e.g. `exclude = "DIFP"`). Unknown names are ignored with a warning.
#' @param metric `"placebo"` (placebo-ATT RMSE, all methods) or `"prediction"`
#'   (per-cell held-out counterfactual RMSE for `"DID"`/`"MC"`/`"TROP"`/`"DIFP"`,
#'   and `"SDID"`/`"SC"` when \pkg{synthdid} is installed; `"gsynth"`,
#'   `"augsynth"` and `"CS"` are `NA` for this metric).
#' @param horizon Number of final periods held out per placebo cohort.
#' @param n_pseudo Number of placebo (pseudo-treated) control units per run.
#' @param n_runs Number of placebo runs to average over.
#' @param control A list of solver/CV controls from [trop_control()].
#' @param seed Optional integer seed for reproducible placebo draws.
#' @param verbose Logical; print progress.
#' @return A `cf_rmse_tbl` (a `data.frame`) with one row per method and columns
#'   `method`, `rmse`, `rmse_se`, `n_runs`, `engine`, `note`.
#' @references Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2025/2026).
#'   Triply/Doubly Robust Panel Estimators.
#' @seealso [panel_compare()], [autoplot.cf_rmse_tbl()]
#' @export
#' @examples
#' \donttest{
#' df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 1)
#' r <- panel_rmse(df, "y", "w", "id", "t",
#'                 methods = c("DID", "TROP"),
#'                 horizon = 2, n_pseudo = 3, n_runs = 2,
#'                 control = trop_control(n_cv_cells = 8L, cv_cycles = 1L),
#'                 seed = 1)
#' r
#' autoplot(r)
#' }
panel_rmse <- function(data, outcome, treatment, unit, time,
                       methods = c("DID", "SC", "SDID", "MC", "TROP", "DIFP"),
                       exclude = NULL,
                       metric = c("placebo", "prediction"),
                       horizon = 10L, n_pseudo = 10L, n_runs = 10L,
                       control = trop_control(), seed = NULL, verbose = FALSE) {
  metric <- match.arg(metric)
  known <- c("DID", "SC", "SDID", "MC", "DIFP", "TROP", "gsynth", "augsynth", "CS")
  methods <- .resolve_methods(methods, exclude, known)

  m <- .panel_to_matrices(data, outcome, treatment, unit, time)
  Y <- m$Y; W <- m$W
  pat <- .assignment_pattern(W)
  N <- nrow(Y); Tt <- ncol(Y); cl <- control$conf_level
  if (horizon < 1L || horizon >= Tt) stop("horizon must be in 1..T-1.")

  controls <- setdiff(seq_len(N), pat$treated_units)
  if (length(controls) < n_pseudo + 1L) {
    stop("Not enough control units (", length(controls),
         ") for n_pseudo = ", n_pseudo, ".", call. = FALSE)
  }
  held_periods <- (Tt - horizon + 1L):Tt

  eng_label <- function(meth) switch(meth,
    DID = , MC = , TROP = , DIFP = "cfcompare", SDID = , SC = "synthdid",
    gsynth = "gsynth", augsynth = "augsynth", CS = "did")

  # ---- tune native methods once on the real data --------------------------
  native_cfg <- list(DID = list(lam = list(time = 0, unit = 0, nn = Inf)))
  if ("MC" %in% methods) {
    native_cfg$MC <- list(lam = .engine_mc(Y, W, pat, control, "none", cl)$lambda)
  }
  if ("TROP" %in% methods) {
    native_cfg$TROP <- list(lam = .engine_trop(Y, W, pat, control, "pooled",
                                               "none", cl)$lambda)
  }

  if (!is.null(seed)) {
    old <- .Random.seed_safe(); on.exit(.Random.seed_restore(old), add = TRUE)
    set.seed(seed)
  }

  rows <- list()

  if (metric == "prediction") {
    # Per-cell held-out prediction RMSE (the paper's Table-1 quantity): hold out
    # a block of control cells, have each method predict the untreated outcome
    # there, and score sqrt(mean((Y - Yhat)^2)). Native engines (DID/MC/TROP) fit
    # on the full panel with the real-treated and held-out cells masked. The
    # block estimators (DIFP, and SDID/SC via synthdid) only define a
    # counterfactual for a clean treated block, so they run on the control-only
    # panel with the held units pseudo-treated; both score the same held cells.
    pred_block <- c("DIFP", "SDID", "SC")          # via a counterfactual matrix
    pred_native <- c("DID", "MC", "TROP")
    pred_unsupported <- c("gsynth", "augsynth", "CS")
    use_methods <- intersect(c(pred_native, pred_block), methods)
    # draw all pseudo-cohorts up front (reproducible given seed); the per-run
    # work below is deterministic and runs in parallel when workers > 1.
    pseudo_list <- lapply(seq_len(n_runs), function(r) sample(controls, n_pseudo))
    par <- (control$workers %||% 1L) > 1L
    one_pred_run <- function(pseudo) {
      out <- stats::setNames(rep(NA_real_, length(use_methods)), use_methods)
      # --- native path: full panel, mask real-treated + held cells ----------
      nat <- intersect(pred_native, use_methods)
      if (length(nat)) {
        excl <- W; excl[pseudo, held_periods] <- 1
        fit_mask <- 1 - excl
        H <- cbind(rep(pseudo, each = length(held_periods)),
                   rep(held_periods, times = length(pseudo)))
        for (meth in nat) {
          lam <- native_cfg[[meth]]$lam
          wmat <- if (meth == "TROP")
            .oos_weight_matrix(Y, excl, pseudo, held_periods, lam) else
            matrix(1, N, Tt)
          Mhat <- .oos_predict_native(Y, fit_mask, wmat, lam$nn, control)
          out[[meth]] <- sqrt(mean((Y[H] - Mhat[H])^2))
        }
      }
      # --- block path: control-only panel, held units pseudo-treated --------
      blk <- intersect(pred_block, use_methods)
      if (length(blk)) {
        Yc <- Y[controls, , drop = FALSE]
        pc <- match(pseudo, controls)
        Wp <- matrix(0, nrow(Yc), Tt, dimnames = dimnames(Yc))
        Wp[pc, held_periods] <- 1
        patp <- .assignment_pattern(Wp)
        Hc <- cbind(rep(pc, each = length(held_periods)),
                    rep(held_periods, times = length(pc)))
        cf_of <- function(meth) switch(meth,
          DIFP = .difp_counterfactual(Yc, Wp, patp),
          SDID = { e <- .engine_synthdid(Yc, Wp, patp, "sdid", cl, se = "none")
                   if (inherits(e, "skip")) NULL else e$counterfactual },
          SC   = { e <- .engine_synthdid(Yc, Wp, patp, "sc", cl, se = "none")
                   if (inherits(e, "skip")) NULL else e$counterfactual })
        for (meth in blk) {
          Mh <- tryCatch(cf_of(meth), error = function(e) NULL)
          if (!is.null(Mh) && all(is.finite(Mh[Hc])))
            out[[meth]] <- sqrt(mean((Yc[Hc] - Mh[Hc])^2))
        }
      }
      out
    }
    if (verbose) message("RMSE: ", n_runs, " prediction run(s)")
    runs <- .with_workers(control$workers %||% 1L,
                          .par_lapply(pseudo_list, one_pred_run, parallel = par))
    run_rmse <- do.call(rbind, runs)
    if (is.null(dim(run_rmse)))
      run_rmse <- matrix(run_rmse, nrow = n_runs,
                         dimnames = list(NULL, use_methods))
    for (meth in use_methods) {
      v <- run_rmse[, meth]; v <- v[is.finite(v)]
      if (length(v) >= 1L) {
        rows[[meth]] <- data.frame(method = meth, rmse = sqrt(mean(v^2)),
          rmse_se = if (length(v) >= 2L) stats::sd(v) / sqrt(length(v)) else NA_real_,
          n_runs = length(v), engine = eng_label(meth), note = NA_character_,
          stringsAsFactors = FALSE)
      } else {
        rows[[meth]] <- data.frame(method = meth, rmse = NA_real_,
          rmse_se = NA_real_, n_runs = 0L, engine = eng_label(meth),
          note = if (meth %in% pred_block)
            "no finite predictions (needs a block design / the 'synthdid' package)"
            else "no finite predictions",
          stringsAsFactors = FALSE)
      }
    }
    for (meth in intersect(pred_unsupported, methods)) {
      rows[[meth]] <- data.frame(method = meth, rmse = NA_real_,
        rmse_se = NA_real_, n_runs = 0L, engine = eng_label(meth),
        note = "metric='prediction' not implemented for this method; use metric='placebo'",
        stringsAsFactors = FALSE)
    }

  } else {
    # placebo-ATT RMSE: assign placebo treatment to control units (true effect
    # 0) and score each method by sqrt(mean(ATT^2)). Works for every method.
    Yc <- Y[controls, , drop = FALSE]; ncy <- nrow(Yc)
    ctrl_ids <- rownames(Y)[controls]; tlabs <- suppressWarnings(as.numeric(colnames(Y)))
    if (anyNA(tlabs)) tlabs <- seq_len(Tt)
    long_unit <- rep(ctrl_ids, times = Tt)
    long_time <- rep(tlabs, each = ncy)
    long_y <- as.numeric(Yc)
    # draw all placebo cohorts up front (reproducible given seed); each run is
    # otherwise self-contained and runs in parallel when workers > 1.
    ps_list <- lapply(seq_len(n_runs), function(r) sample(seq_len(ncy), n_pseudo))
    par <- (control$workers %||% 1L) > 1L
    one_placebo_run <- function(ps) {
      Wp <- matrix(0, ncy, Tt, dimnames = dimnames(Yc))
      Wp[ps, held_periods] <- 1
      patp <- .assignment_pattern(Wp)
      dcp <- data.frame(.u = long_unit, .t = long_time, .y = long_y,
                        .w = as.numeric(Wp), stringsAsFactors = FALSE)
      take <- function(e) if (inherits(e, "skip")) NA_real_ else e$estimate
      att <- stats::setNames(rep(NA_real_, length(methods)), methods)
      notes <- stats::setNames(rep(NA_character_, length(methods)), methods)
      for (meth in methods) {
        val <- tryCatch(switch(meth,
          DID  = .trop_att(Yc, Wp, native_cfg$DID$lam,  control, "pooled", patp)$att,
          MC   = .trop_att(Yc, Wp, native_cfg$MC$lam,   control, "pooled", patp)$att,
          TROP = .trop_att(Yc, Wp, native_cfg$TROP$lam, control, "pooled", patp)$att,
          SDID = { e <- .engine_synthdid(Yc, Wp, patp, "sdid", cl, se = "none")
                   if (inherits(e, "skip")) notes[[meth]] <- e$note; take(e) },
          SC   = { e <- .engine_synthdid(Yc, Wp, patp, "sc", cl, se = "none")
                   if (inherits(e, "skip")) notes[[meth]] <- e$note; take(e) },
          gsynth   = { e <- .engine_gsynth(dcp, ".y", ".w", ".u", ".t", "mc", cl)
                   if (inherits(e, "skip")) notes[[meth]] <- e$note; take(e) },
          augsynth = { e <- .engine_augsynth(dcp, ".y", ".w", ".u", ".t", patp, cl)
                   if (inherits(e, "skip")) notes[[meth]] <- e$note; take(e) },
          CS   = { e <- .engine_did_cs(dcp, ".y", ".w", ".u", ".t", patp, cl)
                   if (inherits(e, "skip")) notes[[meth]] <- e$note; take(e) },
          DIFP = .difp_att(Yc, Wp, patp)
        ), error = function(err) NA_real_)
        att[[meth]] <- val
      }
      list(att = att, notes = notes)
    }
    if (verbose) message("placebo: ", n_runs, " run(s)")
    runs <- .with_workers(control$workers %||% 1L,
                          .par_lapply(ps_list, one_placebo_run, parallel = par))
    run_att <- do.call(rbind, lapply(runs, `[[`, "att"))
    if (is.null(dim(run_att)))
      run_att <- matrix(run_att, nrow = n_runs, dimnames = list(NULL, methods))
    note_mat <- do.call(rbind, lapply(runs, `[[`, "notes"))
    skip_notes <- vapply(methods, function(meth) {
      v <- note_mat[, meth]; v <- v[!is.na(v)]
      if (length(v)) v[1] else NA_character_
    }, character(1))
    for (meth in methods) {
      a <- run_att[, meth]; a <- a[is.finite(a)]
      if (length(a) >= 1) {
        rmse <- sqrt(mean(a^2))
        rse <- if (length(a) >= 2) stats::sd(a^2) / (2 * rmse * sqrt(length(a))) else NA_real_
        rows[[meth]] <- data.frame(method = meth, rmse = rmse, rmse_se = rse,
          n_runs = length(a), engine = eng_label(meth), note = NA_character_,
          stringsAsFactors = FALSE)
      } else {
        rows[[meth]] <- data.frame(method = meth, rmse = NA_real_,
          rmse_se = NA_real_, n_runs = 0L, engine = eng_label(meth),
          note = if (is.na(skip_notes[[meth]])) "no finite placebo estimates" else skip_notes[[meth]],
          stringsAsFactors = FALSE)
      }
    }
  }

  out <- do.call(rbind, rows[methods[methods %in% names(rows)]])
  rownames(out) <- NULL
  attr(out, "horizon") <- horizon
  attr(out, "n_pseudo") <- n_pseudo
  attr(out, "outcome") <- outcome
  attr(out, "metric") <- metric
  class(out) <- c("cf_rmse_tbl", "data.frame")
  out
}

#' @export
print.cf_rmse_tbl <- function(x, digits = 4, ...) {
  cat("Out-of-sample RMSE by method",
      sprintf("(horizon = %s, n_pseudo = %s)\n",
              attr(x, "horizon"), attr(x, "n_pseudo")))
  df <- as.data.frame(x)
  for (nm in c("rmse", "rmse_se")) if (nm %in% names(df))
    df[[nm]] <- round(df[[nm]], digits)
  ord <- order(df$rmse, na.last = TRUE)
  print.data.frame(df[ord, , drop = FALSE], row.names = FALSE)
  invisible(x)
}

#' Bar chart of out-of-sample RMSE across methods
#'
#' Visualises a [panel_rmse()] result as a ranked bar chart (lowest RMSE first),
#' with +/- 1 standard-error whiskers across placebo runs. This is the
#' cross-model RMSE comparison from the paper.
#'
#' @param object A `cf_rmse_tbl` from [panel_rmse()].
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
autoplot.cf_rmse_tbl <- function(object, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("install 'ggplot2' to plot.", call. = FALSE)
  }
  df <- as.data.frame(object)
  df <- df[!is.na(df$rmse), , drop = FALSE]
  if (!nrow(df)) stop("No methods with a computed RMSE to plot.", call. = FALSE)
  df$method <- stats::reorder(df$method, df$rmse)
  best <- df$method[which.min(df$rmse)]
  df$is_best <- df$method == best
  lo <- df$rmse - df$rmse_se; hi <- df$rmse + df$rmse_se
  metric <- attr(object, "metric") %||% "placebo"
  if (identical(metric, "prediction")) {
    rmse_lab  <- "Out-of-sample prediction RMSE"
    title_lab <- "Out-of-sample Prediction RMSE by Method"
    sub_lab   <- sprintf("Horizon = %s periods, %s runs",
                         attr(object, "horizon"), df$n_runs[1])
  } else {
    rmse_lab  <- "Placebo RMSE"
    title_lab <- "Placebo RMSE by Method"
    sub_lab   <- sprintf("Horizon = %s periods, %s placebo runs",
                         attr(object, "horizon"), df$n_runs[1])
  }
  ggplot2::ggplot(df, ggplot2::aes(x = .data$method, y = .data$rmse,
                                   fill = .data$is_best)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.2,
                           na.rm = TRUE) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#2c7fb8", `FALSE` = "grey60"),
                               guide = "none") +
    ggplot2::labs(
      x = NULL, y = rmse_lab,
      title = title_lab,
      subtitle = sub_lab) +
    ggplot2::theme_minimal() +
    .center_titles()
}

#' @export
plot.cf_rmse_tbl <- function(x, ...) {
  print(autoplot.cf_rmse_tbl(x, ...))
  invisible(x)
}
