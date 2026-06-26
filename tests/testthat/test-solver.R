test_that("SVT thresholds singular values correctly", {
  set.seed(1)
  X <- matrix(rnorm(30), 6, 5)
  expect_equal(cfcompare:::.svt(X, 0)$L, X)
  expect_true(all(cfcompare:::.svt(X, Inf)$L == 0))
  expect_equal(cfcompare:::.svt(X, Inf)$rank, 0L)
  # monotone rank in tau
  r_small <- cfcompare:::.svt(X, 0.1)$rank
  r_big   <- cfcompare:::.svt(X, 5)$rank
  expect_true(r_big <= r_small)
})

test_that("weighted MC-NNM reduces to TWFE when lambda = Inf", {
  df <- sim_panel(N = 18, T = 10, n_treated = 3, t0 = 8, seed = 11)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  mask <- (m$W == 0) * 1
  w <- matrix(1, nrow(m$Y), ncol(m$Y))
  fit <- cfcompare:::.mcnnm_fit(m$Y, mask, w, Inf)
  expect_equal(fit$rank, 0L)
  expect_true(all(fit$L == 0))
})

test_that("finite lambda yields a low-rank component", {
  df <- sim_panel(N = 20, T = 12, n_treated = 3, t0 = 9, seed = 12)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  mask <- (m$W == 0) * 1
  w <- matrix(1, nrow(m$Y), ncol(m$Y))
  s1 <- max(svd(m$Y - rowMeans(m$Y))$d)
  fit <- cfcompare:::.mcnnm_fit(m$Y, mask, w, s1 * 0.1)
  expect_gt(fit$rank, 0L)
})
