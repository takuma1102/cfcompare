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
.svt <- function(X, tau, method = c("truncated", "full")) {
  method <- match.arg(method)
  if (is.infinite(tau)) {
    return(list(L = matrix(0, nrow(X), ncol(X)), rank = 0L, d = numeric(0)))
  }
  # Fast path (the default, `method = "truncated"`): when the nuclear-norm rank
  # is small, a truncated SVD (RSpectra) computes only the leading singular
  # triplets and is much faster than a full svd() on a large matrix. It is used
  # only when RSpectra is installed and the matrix is large enough to benefit;
  # on tiny matrices, when RSpectra is absent, or when `method = "full"`, the
  # exact full svd() below is used instead. The two give the same SVT to
  # numerical tolerance, so the choice is a speed/exactness trade-off only.
  # `method = "full"` forces the full SVD and is used for the numerical-agreement
  # checks against the official Python package (see README).
  if (identical(method, "truncated") &&
      min(dim(X)) > 50L &&
      requireNamespace("RSpectra", quietly = TRUE)) {
    res <- .svt_truncated(X, tau)
    if (!is.null(res)) return(res)               # else fall through to full svd
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

#' Truncated singular-value soft-thresholding via RSpectra.
#' Computes only the leading singular triplets, doubling the count until the
#' smallest computed singular value falls at/below `tau` (so every value above
#' the threshold is captured -- giving the same SVT as a full SVD, to numerical
#' tolerance). Returns NULL on failure so the caller can fall back to svd().
#' @keywords internal
#' @noRd
.svt_truncated <- function(X, tau) {
  mn <- min(dim(X))
  k <- min(20L, mn - 1L)
  repeat {
    s <- tryCatch(RSpectra::svds(X, k = k), error = function(e) NULL)
    if (is.null(s) || length(s$d) == 0L) return(NULL)
    if (min(s$d) <= tau || k >= mn - 1L) break
    k <- min(k * 2L, mn - 1L)
  }
  d <- pmax(s$d - tau, 0)
  r <- sum(d > 0)
  if (r == 0L) return(list(L = matrix(0, nrow(X), ncol(X)), rank = 0L, d = d))
  idx <- seq_len(r)
  L <- s$u[, idx, drop = FALSE] %*% (d[idx] * t(s$v[, idx, drop = FALSE]))
  list(L = L, rank = r, d = d)
}

#' Two-way within transformation (double demeaning) of a complete matrix.
#'
#' Subtracts row means, column means, and adds back the grand mean. For a
#' balanced (complete) matrix this is the exact residual maker for additive
#' two-way fixed effects, used to partial the fixed effects out of the
#' covariates (Frisch-Waugh-Lovell) in the covariate-augmented solver.
#' @keywords internal
#' @noRd
.double_demean <- function(M) {
  M - rowMeans(M) -
    matrix(colMeans(M), nrow(M), ncol(M), byrow = TRUE) + mean(M)
}

#' Normalise the `covariates` argument to a list of N x T matrices.
#'
#' Accepts `NULL` (no covariates), a single N x T matrix, an N x T x K array, or
#' a list of N x T matrices; validates the dimensions and that every covariate is
#' fully observed.
#' @keywords internal
#' @noRd
.as_cov_list <- function(X, N, Tt) {
  if (is.null(X)) return(list())
  if (is.list(X) && !is.data.frame(X)) {
    Xl <- X
  } else if (is.array(X) && length(dim(X)) == 3L) {
    Xl <- lapply(seq_len(dim(X)[3L]), function(k) X[, , k])
  } else if (is.matrix(X)) {
    Xl <- list(X)
  } else {
    stop("`covariates` must be a matrix, an N x T x K array, or a list of ",
         "matrices.", call. = FALSE)
  }
  for (k in seq_along(Xl)) {
    Mk <- Xl[[k]]
    if (!is.matrix(Mk) || nrow(Mk) != N || ncol(Mk) != Tt) {
      stop("Each covariate must be an ", N, " x ", Tt, " matrix.", call. = FALSE)
    }
    if (anyNA(Mk)) {
      stop("Covariates must be fully observed (no NA), including treated cells.",
           call. = FALSE)
    }
  }
  Xl
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
#' When `X` is supplied (a list of N x T covariate matrices, or an N x T x K
#' array), the low-rank term is split additively as \eqn{L_{js} = X_{js}'\phi +
#' R_{js}} (eq. of Section 6.2): the covariate-linear part is unpenalised and the
#' nuclear norm is applied to the residual \eqn{R}. Within each soft-impute
#' iteration the fixed effects and \eqn{\phi} are obtained jointly by least
#' squares (the fixed effects are partialled out of the covariates via the
#' two-way within transformation, Frisch-Waugh-Lovell), and the residual is
#' soft-thresholded. With no covariates the computation is unchanged.
#'
#' @param Y N x T outcome matrix (may contain `NA`).
#' @param mask N x T 0/1 matrix; 1 = cell is used in the loss.
#' @param w N x T non-negative observation weights (only used where `mask == 1`).
#' @param lambda Nuclear-norm penalty (use `Inf` for no low-rank term).
#' @param max_iter Maximum number of outer iterations.
#' @param tol Relative Frobenius convergence tolerance on the fitted matrix.
#' @param L_init Optional warm start for the low-rank part (the residual `R` when
#'   covariates are present), e.g. from a previous solve on a similar problem;
#'   speeds up convergence without changing the solution.
#' @param svd_method Singular-value decomposition used inside the soft-impute
#'   step: `"truncated"` (default; leading triplets via RSpectra when available
#'   and worthwhile) or `"full"` (exact base R `svd()`).
#' @param X Optional covariates: `NULL`, a single N x T matrix, an N x T x K
#'   array, or a list of N x T matrices. Each must be fully observed.
#' @return A list with the fitted matrix `M`, low-rank part `L` (the residual
#'   `R` when covariates are present), fixed effects `alpha`/`beta`/`grand`,
#'   covariate coefficients `phi` (length-K, empty if no covariates), estimated
#'   `rank`, iterations `iter`, and `lambda`.
#' @keywords internal
#' @noRd
.mcnnm_fit <- function(Y, mask, w, lambda,
                       max_iter = 200L, tol = 1e-5, L_init = NULL,
                       svd_method = c("truncated", "full"), X = NULL) {
  svd_method <- match.arg(svd_method)
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
  # Soft-threshold for the nuclear-norm prox. The loss is the natural
  # sum-of-squares  sum_js w_js (Y_js - a_j - b_s - L_js)^2  (paper eq. (2); no
  # 1/2 factor), matching the official Python and Stata packages. The gradient
  # step below moves M by  w (Y - M) / Lip, i.e. a step of 1 / (2 * Lip) on that
  # loss (whose gradient has Lipschitz constant 2 * Lip), so the matching
  # threshold is lambda / (2 * Lip). Using lambda / Lip would instead solve the
  # 1/2-loss variant and make lambda_nn half of the paper's / Python's / Stata's
  # scale.
  thr <- lambda / (2 * Lip)
  rnk <- 0L
  it <- 0L

  # ---- covariate preprocessing (Section 6.2: L = X phi + R) ----------------
  Xlist <- .as_cov_list(X, N, Tt)
  K <- length(Xlist)
  phi <- numeric(K)
  if (K > 0L) {
    # raw covariate design (for forming X phi) and its two-way-within version
    # (for estimating phi after partialling out the fixed effects). The within
    # design is constant across iterations, so factor it once.
    Xmat <- matrix(unlist(lapply(Xlist, as.numeric), use.names = FALSE),
                   ncol = K)
    Xdd  <- matrix(unlist(lapply(Xlist, function(M) as.numeric(.double_demean(M))),
                          use.names = FALSE), ncol = K)
    qrXdd <- qr(Xdd)
  }
  xphi_of <- function(phi) if (K > 0L) matrix(Xmat %*% phi, N, Tt) else 0

  M <- g + outer(a, ones_T) + outer(ones_N, b) + xphi_of(phi) + L
  for (it in seq_len(max_iter)) {
    M_old <- M
    resid <- ww * (Yf - M)            # zero on excluded cells
    Mc <- M + resid / Lip
    P <- Mc - L                       # target for fixed effects + X phi
    if (K > 0L) {
      phi <- qr.coef(qrXdd, as.numeric(.double_demean(P)))
      phi[is.na(phi)] <- 0            # covariates collinear with the FE drop out
      Xphi <- matrix(Xmat %*% phi, N, Tt)
    } else {
      Xphi <- 0
    }
    R0 <- P - Xphi
    g <- mean(R0)
    a <- rowMeans(R0) - g
    b <- colMeans(R0) - g
    FEm <- g + outer(a, ones_T) + outer(ones_N, b)
    structured <- FEm + Xphi
    sv <- .svt(Mc - structured, thr, method = svd_method)
    L <- sv$L
    rnk <- sv$rank
    M <- structured + L
    denom <- sqrt(sum(M_old^2)) + 1e-12
    if (sqrt(sum((M - M_old)^2)) / denom < tol) break
  }
  list(M = M, L = L, alpha = a, beta = b, grand = g,
       phi = phi, rank = rnk, iter = it, lambda = lambda)
}
