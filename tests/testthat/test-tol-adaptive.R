# The soft-impute solver converges slowly at small lambda_nn, where a fixed
# relative tolerance stops prematurely (under-converging the fit). .mcnnm_fit()
# therefore tightens the tolerance as lambda_nn weakens -- tol_eff =
# tol * min(1, max(lambda / (2 * max(w) * s1), 1e-3)) -- and the trop_control()
# defaults (larger max_iter, smaller base tol) give it room to run. These tests
# lock that behaviour.

test_that("trop_control() defaults enable convergence at small penalties", {
  ctrl <- trop_control()
  expect_equal(ctrl$max_iter, 2000L)
  expect_equal(ctrl$tol, 1e-6)
})

test_that("small lambda_nn reaches the tight-tolerance solution under defaults", {
  df <- sim_panel(N = 36, T = 28, n_treated = 4, t0 = 20, seed = 7)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  mask <- (m$W == 0) * 1
  w <- matrix(1, nrow(m$Y), ncol(m$Y))
  Yf <- m$Y
  Yf[is.na(Yf)] <- 0
  s1 <- max(svd(cfcompare:::.double_demean(Yf))$d)
  nn <- 0.02 * s1                       # weak penalty: the slow-converging regime
  att <- function(fit) mean((m$Y - fit$M)[m$W == 1])

  ref <- cfcompare:::.mcnnm_fit(m$Y, mask, w, nn,
                                tol = 1e-9, max_iter = 8000L, svd_method = "full")
  def <- cfcompare:::.mcnnm_fit(m$Y, mask, w, nn)   # package defaults (nn-adaptive tol)

  # The default fit must reach essentially the same ATT as a tight-tolerance
  # solve. A fixed relative tolerance with the old low iteration cap (the
  # pre-adaptive behaviour) stops far short here, so this also guards against
  # reverting the adaptive tolerance.
  expect_lt(abs(att(def) - att(ref)), 0.05)
})
