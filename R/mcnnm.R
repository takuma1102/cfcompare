# Internal numerical engine for the TROP working model.
# Depends only on base R + stats (svd). No external packages.

#' Singular-value soft-thresholding
#'
#' Computes \eqn{\mathrm{SVT}_\tau(X) = U \, \mathrm{diag}(\max(d - \tau, 0)) \, V'}.
#'
#' @param X A numeric matrix.
#' @param tau Non-negative threshold.
#' @return A list with the thresholded matrix `L`, its `rank`, and the
#'   shrunken singular values `d`.
#' @keywords internal
#' @noRd
.svt <- function(X, tau) {
  if (is.infinite(tau)) {
    return(list(L = matrix(0, nrow(X), ncol(X)), rank = 0L, d = numeric(0)))
  }
  sv <- svd(X)
  d <- pmax(sv$d - tau, 0)
  r <- sum(d > 0)
  if (r == 0L) {
    return(list(L = matrix(0, nrow(X), ncol(X)), rank = 0L, d = d))
  }
  idx <- seq_len(r)
  L <- sv$u[, idx, drop = FALSE] %*%
    (d[idx] * t(sv$v[, idx, drop = FALSE]))
  list(L = L, rank = r, d = d)
}

#' Weighted nuclear-norm matrix completion with two-way fixed effects
#'
#' Solves the TROP working-model program
#' \deqn{\min_{\alpha,\beta,L} \sum_{j,s} m_{js} w_{js}
#'   (Y_{js} - \alpha_j - \beta_s - L_{js})^2 + \lambda \lVert L \rVert_*}
#' via proximal-gradient / soft-impute iterations. Only cells with `mask == 1`
#' enter the loss; the fitted matrix `M` provides counterfactual predictions for
#' every cell, including the excluded (treated) ones.
#'
#' Setting `lambda = Inf` forces `L = 0` and returns the weighted two-way
#' fixed-effects (DID/TWFE) fit; finite `lambda` with uniform weights yields the
#' matrix-completion (MC) estimator. See Athey, Imbens, Qu & Viviano (2025),
#' eq. (2).
#'
#' @param Y N x T outcome matrix (may contain `NA`).
#' @param mask N x T 0/1 matrix; 1 = cell is used in the loss.
#' @param w N x T non-negative observation weights (only used where `mask == 1`).
#' @param lambda Nuclear-norm penalty (use `Inf` for no low-rank term).
#' @param max_iter Maximum number of outer iterations.
#' @param tol Relative Frobenius convergence tolerance on the fitted matrix.
#' @param L_init Optional warm start for `L`.
#' @return A list with the fitted matrix `M`, low-rank part `L`, fixed effects
#'   `alpha`/`beta`/`grand`, estimated `rank`, iterations `iter`, and `lambda`.
#' @keywords internal
#' @noRd
.mcnnm_fit <- function(Y, mask, w, lambda,
                       max_iter = 200L, tol = 1e-5, L_init = NULL) {
  N <- nrow(Y); Tt <- ncol(Y)
  ww <- mask * w
  Lip <- max(ww)
  if (!is.finite(Lip) || Lip <= 0) Lip <- 1
  Yf <- Y
  Yf[is.na(Yf)] <- 0
  ones_T <- rep(1, Tt)
  ones_N <- rep(1, N)
  a <- numeric(N); b <- numeric(Tt); g <- 0
  L <- if (is.null(L_init)) matrix(0, N, Tt) else L_init
  fe_mat <- function() g + outer(a, ones_T) + outer(ones_N, b)
  M <- fe_mat() + L
  thr <- lambda / Lip
  rnk <- 0L
  it <- 0L
  for (it in seq_len(max_iter)) {
    M_old <- M
    resid <- ww * (Yf - M)            # zero on excluded cells
    Mc <- M + resid / Lip
    R <- Mc - L
    g <- mean(R)
    a <- rowMeans(R) - g
    b <- colMeans(R) - g
    FEm <- g + outer(a, ones_T) + outer(ones_N, b)
    sv <- .svt(Mc - FEm, thr)
    L <- sv$L
    rnk <- sv$rank
    M <- FEm + L
    denom <- sqrt(sum(M_old^2)) + 1e-12
    if (sqrt(sum((M - M_old)^2)) / denom < tol) break
  }
  list(M = M, L = L, alpha = a, beta = b, grand = g,
       rank = rnk, iter = it, lambda = lambda)
}
