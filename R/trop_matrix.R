# Matrix-in entry point for the TROP working model, written independently from
# the paper's equations. Its matrix-in form makes it convenient to numerically
# compare results against other matrix-based panel-estimator software on
# identical inputs. It is not a port of, and shares no code with, any official
# TROP implementation.

#' Block-design TROP weights
#'
#' Builds the block-design weights from the paper's equation (3): unit weights
#' use the RMS distance, over pre-periods, to the average treated trajectory;
#' time weights decay with distance to the centre of the treated block at the
#' end of the panel.
#' @keywords internal
#' @noRd
.trop_reference_weights <- function(Y, treated_units, lambda_unit, lambda_time,
                                    treated_periods) {
  N <- nrow(Y); Tt <- ncol(Y); tp <- treated_periods
  dist_time <- abs((seq_len(Tt) - 1) - (Tt - tp / 2))
  avg_tr <- colMeans(Y[treated_units, , drop = FALSE])
  mask <- matrix(1, N, Tt)
  mask[, (Tt - tp + 1):Tt] <- 0
  A <- rowSums(((matrix(avg_tr, N, Tt, byrow = TRUE) - Y)^2) * mask)
  B <- rowSums(mask)
  if (any(B == 0)) stop("Some unit has no pre-periods under treated_periods.")
  du <- sqrt(A / B)
  outer(exp(-lambda_unit * du), exp(-lambda_time * dist_time))
}

#' TROP estimate on matrix input
#'
#' A thin, matrix-in wrapper around the TROP working model, written independently
#' from the paper. The matrix-in form is convenient for numerically comparing
#' [trop()] against other matrix-based implementations on identical inputs. For
#' data-frame input, cross-validated penalty selection, inference and the
#' multi-estimator comparison, use [trop()] and [panel_compare()] instead.
#'
#' The untreated outcome model is fitted on the control cells with the supplied
#' weights, and the returned effect is the average over treated cells of
#' `Y - alpha - beta - L`, exactly as in the paper's eq. (2). With
#' `lambda_nn = Inf` the low-rank term is dropped and the fit is weighted two-way
#' fixed effects, matching the reference to numerical tolerance. With a finite
#' `lambda_nn` the nuclear-norm term is solved by this package's
#' proximal-gradient routine; because that uses a different parameterisation from
#' the reference's convex solver, finite-penalty results agree in behaviour but
#' not to the last digit.
#'
#' @param Y N x T outcome matrix.
#' @param W N x T 0/1 treatment matrix (1 = actively treated cell).
#' @param treated_units Integer row indices of the treated units, used to anchor
#'   the unit weights.
#' @param lambda_unit,lambda_time Non-negative decay parameters for the unit and
#'   time weights.
#' @param lambda_nn Nuclear-norm penalty; use `Inf` to drop the low-rank term.
#' @param treated_periods Number of final columns treated as the post block, used
#'   to build the pre-period mask and the time-distance centre.
#' @param control A list of solver controls from [trop_control()].
#' @return A single numeric ATT estimate.
#' @references Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2025).
#'   \emph{Triply Robust Panel Estimators.} arXiv:2508.21536.
#' @seealso [trop()], [panel_compare()]
#' @export
#' @examples
#' df <- sim_panel(N = 15, T = 12, n_treated = 3, t0 = 9, att = 2, seed = 1)
#' Y <- matrix(0, max(df$id), max(df$t)); W <- Y
#' for (k in seq_len(nrow(df))) { Y[df$id[k], df$t[k]] <- df$y[k]
#'   W[df$id[k], df$t[k]] <- df$w[k] }
#' tu <- which(rowSums(W) > 0)
#' trop_matrix(Y, W, tu, 0.1, 0.1, Inf, treated_periods = 4)
trop_matrix <- function(Y, W, treated_units,
                        lambda_unit, lambda_time, lambda_nn,
                        treated_periods, control = trop_control()) {
  Y <- as.matrix(Y); W <- as.matrix(W)
  if (!identical(dim(Y), dim(W))) stop("Y and W must have the same dimensions.")
  N <- nrow(Y); Tt <- ncol(Y)
  if (treated_periods <= 0 || treated_periods >= Tt)
    stop("treated_periods must be in 1..T-1.")
  if (length(treated_units) == 0L) stop("treated_units must be non-empty.")
  if (lambda_unit < 0 || lambda_time < 0)
    stop("lambda_unit and lambda_time must be non-negative.")
  if (!any(W == 1)) stop("W contains no treated cells.")

  wmat <- .trop_reference_weights(Y, treated_units, lambda_unit, lambda_time,
                                  treated_periods)
  fit <- .trop_solve(Y, W, wmat, list(nn = lambda_nn), control)
  mean((Y - fit$M)[W == 1])
}
