# Tests for the nuclear-norm penalty grid and the loss-scaling convention it is
# tuned to. The MC-NNM solver targets the paper's loss (eq. (2))
#   sum_js w_js (Y_js - alpha_j - beta_s - L_js)^2 + lambda_nn ||L||_*
# with NO 1/2 factor, matching the official Python and Stata packages. Its
# soft-threshold is therefore lambda_nn / (2 * max(weight)), which means that
# with uniform weights the low-rank term vanishes exactly at lambda_nn = 2 * s1,
# where s1 is the largest singular value of the two-way-demeaned control matrix.
# The default grid is scaled to that boundary.

# sigma_max of the two-way-demeaned, mean-imputed control matrix (uniform
# weights) -- the same quantity .trop_default_grids() uses to scale nn.
.demeaned_sigma <- function(Y, W) {
  Yc <- Y
  Yc[W == 1 | is.na(Yc)] <- NA
  Yc[is.na(Yc)] <- mean(Yc, na.rm = TRUE)
  Z <- Yc - rowMeans(Yc)
  Z <- t(t(Z) - colMeans(Z))
  max(svd(Z)$d)
}

test_that("MC-NNM uses the paper's (no-1/2) loss: L vanishes at lambda_nn = 2 * sigma", {
  df <- sim_panel(N = 20, T = 12, n_treated = 3, t0 = 9, seed = 5)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  mask <- (m$W == 0) * 1
  w <- matrix(1, nrow(m$Y), ncol(m$Y))
  sig <- .demeaned_sigma(m$Y, m$W)

  rk <- function(mult) {
    cfcompare:::.mcnnm_fit(m$Y, mask, w, mult * sig,
                           max_iter = 1000L, tol = 1e-9)$rank
  }

  # Paper / Python / Stata convention: threshold = lambda_nn / (2 * Lip), so the
  # low-rank term is fully shrunk once lambda_nn >= 2 * sigma, and is still
  # present anywhere below that.
  expect_equal(rk(2.2), 0L)   # above 2 * sigma -> no low-rank term
  expect_gt(rk(1.5), 0L)      # between sigma and 2 * sigma -> L present
  # A regression to the 1/2-loss threshold (lambda_nn / Lip) would make L vanish
  # already at sigma, so rk(1.5) would wrongly be 0. This guards that fix.
})

test_that("default nuclear-norm grid is well-formed and scaled to the 2*sigma boundary", {
  df <- sim_panel(N = 24, T = 14, n_treated = 4, t0 = 11, seed = 7)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  g <- cfcompare:::.trop_default_grids(m$Y, m$W)
  fin <- g$nn[is.finite(g$nn)]

  expect_true(any(is.infinite(g$nn)))   # Inf (exactly L = 0) is always available
  expect_true(all(fin > 0))
  expect_true(all(diff(fin) < 0))       # strictly decreasing

  # The largest finite penalty sits at ~ s1 (half the 2*s1 vanishing point), so
  # it drives L to a near-zero rank; the smallest keeps a rich low-rank term.
  # Together they span the meaningful regularisation range.
  sig <- .demeaned_sigma(m$Y, m$W)
  expect_equal(max(fin), sig, tolerance = 1e-6)

  mask <- (m$W == 0) * 1
  w <- matrix(1, nrow(m$Y), ncol(m$Y))
  rank_top <- cfcompare:::.mcnnm_fit(m$Y, mask, w, max(fin))$rank
  rank_bot <- cfcompare:::.mcnnm_fit(m$Y, mask, w, min(fin))$rank
  expect_lte(rank_top, 3L)              # top finite penalty nearly kills L
  expect_gt(rank_bot, rank_top)         # grid spans a meaningful range
})
