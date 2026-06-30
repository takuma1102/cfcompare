# Per-method estimation engines used by panel_compare().
#
# Native engines (DID, MC, TROP) are built on the shared .trop_engine() and have
# no external dependencies. SDID/SC wrap the synthdid package; an alternative MC
# can be routed through gsynth. Missing optional packages cause the method to be
# skipped with a message (the nonabsdid pattern), rather than erroring.

#' @keywords internal
#' @noRd
.engine_row <- function(method, engine, eng, pat, outcome, note = NA_character_) {
  data.frame(
    method = method,
    estimate = eng$estimate,
    std.error = eng$std.error %||% NA_real_,
    conf.low = eng$conf.low %||% NA_real_,
    conf.high = eng$conf.high %||% NA_real_,
    n_treated_cells = pat$n_treated_cells,
    n_treated_units = pat$n_treated_units,
    outcome = outcome,
    engine = engine,
    rank = eng$rank %||% NA_integer_,
    note = note,
    stringsAsFactors = FALSE
  )
}

#' Native DID / TWFE engine: lambda_nn = Inf, uniform weights.
#' @keywords internal
#' @noRd
.engine_did <- function(Y, W, pat, control, se, conf_level, X = NULL) {
  .trop_engine(Y, W, pat,
               lambda = list(time = 0, unit = 0, nn = Inf),
               grids = NULL, control = control,
               anchor = "pooled", se = se, conf_level = conf_level, X = X)
}

#' Native matrix-completion engine: uniform weights, CV over lambda_nn only.
#' @keywords internal
#' @noRd
.engine_mc <- function(Y, W, pat, control, se, conf_level, verbose = FALSE, X = NULL) {
  grids <- .trop_default_grids(Y, W, X = X)
  grids$time <- 0
  grids$unit <- 0
  .trop_engine(Y, W, pat, lambda = NULL, grids = grids, control = control,
               anchor = "pooled", se = se, conf_level = conf_level,
               verbose = verbose, X = X)
}

#' Native TROP engine (full).
#' @keywords internal
#' @noRd
.engine_trop <- function(Y, W, pat, control, anchor, se, conf_level,
                         verbose = FALSE, X = NULL) {
  if (anchor == "auto") {
    anchor <- if (pat$n_treated_cells <= control$max_cells) "per_cell" else "pooled"
  }
  grids <- .trop_default_grids(Y, W, X = X)
  .trop_engine(Y, W, pat, lambda = NULL, grids = grids, control = control,
               anchor = anchor, se = se, conf_level = conf_level,
               verbose = verbose, X = X)
}

#' synthdid engine (SDID or SC), block designs only.
#' @keywords internal
#' @noRd
.engine_synthdid <- function(Y, W, pat, which = c("sdid", "sc"), conf_level,
                             se = "jackknife") {
  which <- match.arg(which)
  if (!requireNamespace("synthdid", quietly = TRUE)) {
    return(structure(list(skip = TRUE,
      note = "install 'synthdid' to enable this method"), class = "skip"))
  }
  if (pat$type != "block") {
    return(structure(list(skip = TRUE,
      note = "synthdid requires a block design"), class = "skip"))
  }
  # synthdid wants units ordered control-first, treated-last; periods pre-first.
  N0 <- nrow(Y) - pat$n_treated_units
  T0 <- pat$block_t0 - 1L
  ord_u <- c(setdiff(seq_len(nrow(Y)), pat$treated_units), pat$treated_units)
  Ys <- Y[ord_u, , drop = FALSE]
  est_fun <- if (which == "sdid") synthdid::synthdid_estimate else synthdid::sc_estimate
  fit <- est_fun(Ys, N0 = N0, T0 = T0)
  att <- as.numeric(fit)
  # jackknife vcov refits N times; skip it unless an SE was actually requested
  # (and it needs >= 2 treated units to be defined).
  want_se <- !identical(se, "none") && pat$n_treated_units >= 2
  se_val <- NA_real_
  if (want_se) {
    v <- tryCatch(stats::vcov(fit, method = "jackknife"),
                  error = function(e) NA_real_)
    se_val <- if (is.matrix(v)) sqrt(v[1, 1]) else as.numeric(v)
  }
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  # Reconstruct the treated-cell counterfactual from the fitted weights. For SC
  # this is the omega-weighted control path; for SDID it adds the lambda-weighted
  # DID level correction. Guarded: any failure yields NULL (no counterfactual).
  # Not exercised in-package (synthdid is an optional dependency).
  cf <- tryCatch({
    wts <- attr(fit, "weights")
    omega <- wts$omega; lambda <- wts$lambda
    base <- as.numeric(crossprod(Ys[seq_len(N0), , drop = FALSE], omega))  # length T
    tr_avg <- colMeans(Ys[(N0 + 1L):nrow(Ys), , drop = FALSE])
    icpt <- if (which == "sdid")
      sum(lambda * (tr_avg[seq_len(T0)] - base[seq_len(T0)])) else 0
    cf_path <- base + icpt
    M <- matrix(NA_real_, nrow(Y), ncol(Y), dimnames = dimnames(Y))
    for (i in pat$treated_units) M[i, ] <- cf_path
    M
  }, error = function(e) NULL)
  list(estimate = att, std.error = se_val,
       conf.low = att - z * se_val, conf.high = att + z * se_val,
       se.method = if (want_se) "jackknife" else "none", rank = NA_integer_,
       counterfactual = cf, fit = fit)
}

#' gsynth engine (alternative MC / IFE), optional.
#' @keywords internal
#' @noRd
.engine_gsynth <- function(data, outcome, treatment, unit, time,
                           estimator = "mc", conf_level) {
  if (!requireNamespace("gsynth", quietly = TRUE)) {
    return(structure(list(skip = TRUE,
      note = "install 'gsynth' to enable this method"), class = "skip"))
  }
  fml <- stats::as.formula(paste(outcome, "~", treatment))
  fit <- gsynth::gsynth(fml, data = data, index = c(unit, time),
                        estimator = estimator, se = TRUE, inference = "parametric",
                        nboots = 200, parallel = FALSE)
  att <- fit$att.avg
  se <- tryCatch(fit$est.avg[1, "S.E."], error = function(e) NA_real_)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  # gsynth stores the estimated treated counterfactual in fit$Y.ct (periods x
  # treated units). Map it back onto an N x T matrix; guarded (optional dep).
  cf <- tryCatch({
    Yct <- fit$Y.ct
    trn <- fit$id.tr %||% colnames(Yct)
    M <- matrix(NA_real_, length(fit$id), length(fit$time),
                dimnames = list(as.character(fit$id), as.character(fit$time)))
    ri <- match(as.character(trn), rownames(M))
    M[ri, ] <- t(Yct)
    M
  }, error = function(e) NULL)
  list(estimate = as.numeric(att), std.error = as.numeric(se),
       conf.low = as.numeric(att) - z * se,
       conf.high = as.numeric(att) + z * se,
       se.method = "parametric", rank = NA_integer_,
       counterfactual = cf, fit = fit)
}

#' augsynth engine (Augmented Synthetic Control), block designs only, optional.
#' @keywords internal
#' @noRd
.engine_augsynth <- function(data, outcome, treatment, unit, time, pat,
                             conf_level) {
  if (!requireNamespace("augsynth", quietly = TRUE)) {
    return(structure(list(skip = TRUE,
      note = "install 'augsynth' (ebenmichael/augsynth) to enable this method"),
      class = "skip"))
  }
  if (pat$type != "block") {
    return(structure(list(skip = TRUE,
      note = "augsynth (single-period) requires a block design"), class = "skip"))
  }
  # augsynth: augsynth(outcome ~ trt, unit, time, t_int, data); t_int = first
  # treated period. Build the formula and call with bare column names.
  t_int <- pat$block_t0
  fml <- stats::as.formula(paste(outcome, "~", treatment))
  fit <- tryCatch(
    augsynth::augsynth(fml, unit = as.name(unit), time = as.name(time),
                       t_int = t_int, data = data,
                       progfunc = "None", scm = TRUE),
    error = function(e) e)
  if (inherits(fit, "error")) {
    return(structure(list(skip = TRUE,
      note = paste0("augsynth failed: ", conditionMessage(fit))), class = "skip"))
  }
  s <- summary(fit, inf_type = "jackknife")
  avg <- s$average_att
  att <- as.numeric(avg$Estimate[1])
  se <- tryCatch(as.numeric(avg$Std.Error[1]), error = function(e) NA_real_)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  list(estimate = att, std.error = se,
       conf.low = att - z * se, conf.high = att + z * se,
       se.method = "jackknife", rank = NA_integer_, fit = fit)
}

#' Callaway & Sant'Anna engine (did package), staggered/absorbing designs.
#' @keywords internal
#' @noRd
.engine_did_cs <- function(data, outcome, treatment, unit, time, pat,
                           conf_level) {
  if (!requireNamespace("did", quietly = TRUE)) {
    return(structure(list(skip = TRUE,
      note = "install 'did' (Callaway & Sant'Anna) to enable this method"),
      class = "skip"))
  }
  d <- as.data.frame(data)
  # derive group (first treated period) per unit; require absorbing treatment.
  sp <- split(d, d[[unit]])
  gmap <- vapply(sp, function(g) {
    g <- g[order(g[[time]]), ]
    tr <- which(g[[treatment]] == 1)
    if (!length(tr)) return(0)               # never treated
    first <- min(tr)
    # absorbing check: treated from first period onward
    if (!all(g[[treatment]][first:nrow(g)] == 1)) return(NA_real_)
    g[[time]][first]
  }, numeric(1))
  if (anyNA(gmap)) {
    return(structure(list(skip = TRUE,
      note = "did (Callaway-Sant'Anna) requires absorbing treatment"),
      class = "skip"))
  }
  if (all(gmap == 0)) {
    return(structure(list(skip = TRUE,
      note = "no treated groups for did"), class = "skip"))
  }
  d$.g <- gmap[as.character(d[[unit]])]
  cg <- if (any(gmap == 0)) "nevertreated" else "notyettreated"
  fit <- tryCatch({
    out <- did::att_gt(yname = outcome, tname = time, idname = unit,
                       gname = ".g", data = d, control_group = cg,
                       bstrap = TRUE, cband = FALSE, alp = 1 - conf_level)
    did::aggte(out, type = "simple")
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    return(structure(list(skip = TRUE,
      note = paste0("did failed: ", conditionMessage(fit))), class = "skip"))
  }
  att <- as.numeric(fit$overall.att)
  se <- as.numeric(fit$overall.se)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  list(estimate = att, std.error = se,
       conf.low = att - z * se, conf.high = att + z * se,
       se.method = "did-multiplier-bootstrap", rank = NA_integer_, fit = fit)
}

#' Simplex-constrained least squares (synthetic-control weights) via Frank-Wolfe.
#' Minimises ||a - X w||^2 over w >= 0, sum(w) = 1. Written from first principles.
#' @keywords internal
#' @noRd
.sc_simplex <- function(X, a, max_iter = 2000L, tol = 1e-8) {
  J <- ncol(X)
  if (J == 0L) return(numeric(0))
  w <- rep(1 / J, J)
  Xw <- as.numeric(X %*% w)
  for (it in seq_len(max_iter)) {
    resid <- a - Xw                      # n_pre vector
    grad <- -as.numeric(crossprod(X, resid))   # d/dw of ||a - Xw||^2 (up to 2)
    k <- which.min(grad)                 # vertex of the simplex
    s <- numeric(J); s[k] <- 1
    d <- s - w
    Xd <- as.numeric(X %*% d)
    denom <- sum(Xd * Xd)
    if (denom <= .Machine$double.eps) break
    gamma <- sum(resid * Xd) / denom
    gamma <- min(max(gamma, 0), 1)
    if (gamma < tol) break
    w <- w + gamma * d
    Xw <- Xw + gamma * Xd
  }
  w
}

#' DIFP point estimate: demeaned synthetic control (Ferman-Pinto recentering with
#' a Doudchenko-Imbens intercept). Block designs only. Independent implementation.
#' @keywords internal
#' @noRd
.difp_att <- function(Y, W, pat) {
  if (pat$type != "block") return(NA_real_)
  t0 <- pat$block_t0
  Tt <- ncol(Y)
  pre <- seq_len(t0 - 1L); post <- t0:Tt
  if (length(pre) < 1L) return(NA_real_)
  tu <- pat$treated_units
  co <- setdiff(seq_len(nrow(Y)), tu)
  if (length(co) < 1L) return(NA_real_)
  tr_avg <- colMeans(Y[tu, , drop = FALSE])
  tr_pre_mean <- mean(tr_avg[pre])
  Cc <- Y[co, , drop = FALSE] - rowMeans(Y[co, pre, drop = FALSE])  # recentre each control
  a <- tr_avg[pre] - tr_pre_mean                 # recentred treated pre-path
  X <- t(Cc[, pre, drop = FALSE])                # n_pre x n_control
  w <- .sc_simplex(X, a)
  cf_post <- tr_pre_mean + as.numeric(crossprod(Cc[, post, drop = FALSE], w))
  mean(tr_avg[post] - cf_post)
}

#' DIFP counterfactual matrix: the demeaned-SC group path on treated cells.
#'
#' Returns the N x T untreated-potential-outcome matrix implied by DIFP, with the
#' synthetic-control group path (`.difp_att`'s internal construction) written on
#' every treated row and `NA` elsewhere. Block designs only.
#' @keywords internal
#' @noRd
.difp_counterfactual <- function(Y, W, pat) {
  M <- matrix(NA_real_, nrow(Y), ncol(Y), dimnames = dimnames(Y))
  if (pat$type != "block") return(M)
  t0 <- pat$block_t0; Tt <- ncol(Y)
  pre <- seq_len(t0 - 1L)
  if (length(pre) < 1L) return(M)
  tu <- pat$treated_units
  co <- setdiff(seq_len(nrow(Y)), tu)
  if (length(co) < 1L) return(M)
  tr_avg <- colMeans(Y[tu, , drop = FALSE])
  tr_pre_mean <- mean(tr_avg[pre])
  Cc <- Y[co, , drop = FALSE] - rowMeans(Y[co, pre, drop = FALSE])
  a <- tr_avg[pre] - tr_pre_mean
  Xpre <- t(Cc[, pre, drop = FALSE])
  w <- .sc_simplex(Xpre, a)
  cf_path <- tr_pre_mean + as.numeric(crossprod(Cc, w))   # length Tt
  for (i in tu) M[i, ] <- cf_path
  M
}

#' DIFP engine (demeaned synthetic control) with optional SE.
#' @keywords internal
#' @noRd
.engine_difp <- function(Y, W, pat, control, conf_level, se = "auto") {
  if (pat$type != "block") {
    return(structure(list(skip = TRUE,
      note = "DIFP (demeaned SC) requires a block design"), class = "skip"))
  }
  att <- .difp_att(Y, W, pat)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  if (identical(se, "auto")) se <- if (pat$n_treated_units >= 2) "jackknife" else "none"
  se_val <- NA_real_; meth <- se
  if (se == "jackknife" && pat$n_treated_units >= 2) {
    tu <- pat$treated_units; G <- length(tu); a <- numeric(G)
    for (g in seq_len(G)) {
      keep <- setdiff(seq_len(nrow(Y)), tu[g])
      Yk <- Y[keep, , drop = FALSE]; Wk <- W[keep, , drop = FALSE]
      a[g] <- .difp_att(Yk, Wk, .assignment_pattern(Wk))
    }
    se_val <- sqrt((G - 1) / G * sum((a - mean(a))^2))
  } else {
    meth <- "none"
  }
  # counterfactual matrix (group SC path written on each treated row)
  Mhat <- .difp_counterfactual(Y, W, pat)
  list(estimate = att, std.error = se_val,
       conf.low = att - z * se_val, conf.high = att + z * se_val,
       se.method = meth, rank = NA_integer_, counterfactual = Mhat)
}

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
