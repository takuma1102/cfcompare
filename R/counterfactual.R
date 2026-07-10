# A common accessor for the estimated untreated-potential-outcome ("Y(0)")
# matrix produced by any supported estimator, so the same N x T counterfactual
# can be compared across methods (e.g. the unit-period prediction RMSE of the
# paper's Table 1; see panel_rmse(metric = "prediction")).

#' Estimated counterfactual (untreated potential outcome) matrix
#'
#' Returns a method's estimate of the untreated potential outcome \eqn{Y(0)} as
#' an `N x T` matrix on the panel's unit/time grid. This gives every estimator a
#' common interface, so counterfactuals -- and the unit-period prediction RMSE
#' built on them ([panel_rmse()] with `metric = "prediction"`) -- can be compared
#' on the same footing.
#'
#' Coverage of the returned matrix differs by method. The native engines
#' (`trop()`, and the DID/MC/TROP fits inside [panel_compare()]) return a fitted
#' value for *every* cell. Methods whose model only defines a counterfactual for
#' the treated units -- DIFP, and (when the optional package is installed) SDID,
#' SC and gsynth -- fill the treated rows and leave control cells `NA` (a
#' control's untreated outcome is simply its observed value).
#'
#' @param object A fitted object: a `trop` fit, a `cf_comparison` from
#'   [panel_compare()], a \pkg{synthdid} estimate, or any object carrying a
#'   `counterfactual` element.
#' @param ... Unused; for method extensibility.
#' @return For a single fit, an `N x T` numeric matrix (with `NA` where the
#'   method does not define a counterfactual). For a `cf_comparison`, a named
#'   list of such matrices, one per method.
#' @seealso [panel_rmse()], [trop()], [panel_compare()]
#' @examples
#' df <- sim_panel(N = 14, T = 9, n_treated = 3, t0 = 6, seed = 1)
#' fit <- trop(df, "y", "w", "id", "t",
#'             lambda = list(time = 0.1, unit = 0.5, nn = 2), se = "none")
#' M <- counterfactual_matrix(fit)
#' dim(M)
#' @export
counterfactual_matrix <- function(object, ...) UseMethod("counterfactual_matrix")

#' @rdname counterfactual_matrix
#' @details `predict_counterfactual()` is an alias for `counterfactual_matrix()`.
#' @export
predict_counterfactual <- function(object, ...) counterfactual_matrix(object, ...)

#' @export
counterfactual_matrix.trop <- function(object, ...) {
  Y <- object$panel$Y; W <- object$panel$W
  pat <- object$pattern; lam <- object$lambda
  ctrl <- trop_control()
  # full pooled counterfactual (covariate-adjusted when covariates were used),
  # regardless of the fit's anchor -- matches autoplot.trop() / the event study.
  M <- .trop_pooled_M(Y, W, lam, ctrl, pat, X = object$panel$X)
  # panel$Y is on the fitting scale; map back to the raw outcome scale for
  # standardized fits (identity otherwise).
  sc <- object$scaling %||% list(center = 0, scale = 1)
  sc$center + sc$scale * M
}

#' @export
counterfactual_matrix.cf_comparison <- function(object, ...) {
  Y <- object$panel$Y
  methods <- as.character(object$att$method)
  na_mat <- function() matrix(NA_real_, nrow(Y), ncol(Y), dimnames = dimnames(Y))
  cfs <- object$counterfactual %||% list()
  out <- lapply(methods, function(m) {
    cm <- cfs[[m]]
    if (is.null(cm)) na_mat() else cm
  })
  stats::setNames(out, methods)
}

#' @export
counterfactual_matrix.synthdid_estimate <- function(object,
                                                    which = c("sdid", "sc"), ...) {
  which <- match.arg(which)
  setup <- attr(object, "setup")
  wts <- attr(object, "weights")
  Ys <- setup$Y; N0 <- setup$N0; T0 <- setup$T0
  N <- nrow(Ys); Tt <- ncol(Ys)
  omega <- wts$omega; lambda <- wts$lambda
  base <- as.numeric(crossprod(Ys[seq_len(N0), , drop = FALSE], omega))
  tr_avg <- colMeans(Ys[(N0 + 1L):N, , drop = FALSE])
  icpt <- if (which == "sdid")
    sum(lambda * (tr_avg[seq_len(T0)] - base[seq_len(T0)])) else 0
  cf_path <- base + icpt
  # synthdid reorders units control-first, treated-last; rows are returned in
  # that ordering (the treated units are the final N - N0 rows).
  M <- matrix(NA_real_, N, Tt, dimnames = dimnames(Ys))
  for (i in (N0 + 1L):N) M[i, ] <- cf_path
  M
}

#' @export
counterfactual_matrix.default <- function(object, ...) {
  cm <- tryCatch(object[["counterfactual"]], error = function(e) NULL)
  if (!is.null(cm)) return(cm)
  stop("No counterfactual matrix available for an object of class ",
       paste(class(object), collapse = "/"), ".", call. = FALSE)
}
