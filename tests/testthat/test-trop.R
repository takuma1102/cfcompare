test_that("trop() runs and returns a sensible ATT", {
  df <- sim_panel(N = 22, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 7)
  fit <- trop(df, "y", "w", "id", "t", se = "none",
              control = trop_control(n_cv_cells = 25, seed = 1, cv_cycles = 1))
  expect_s3_class(fit, "trop")
  expect_true(is.finite(fit$estimate))
  expect_true(abs(fit$estimate - 2) < 1.5)
  expect_equal(nrow(fit$tau_cells), fit$pattern$n_treated_cells)
  expect_equal(dim(fit$counterfactual), dim(fit$panel$Y))
})

test_that("fixed lambda skips CV and DID special case matches TWFE engine", {
  df <- sim_panel(N = 18, T = 10, n_treated = 3, t0 = 8, att = 1, seed = 3)
  fit <- trop(df, "y", "w", "id", "t",
              lambda = list(time = 0, unit = 0, nn = Inf),
              anchor = "pooled", se = "none")
  expect_equal(fit$rank, 0L)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  pat <- cfcompare:::.assignment_pattern(m$W)
  did <- cfcompare:::.engine_did(m$Y, m$W, pat, trop_control(), "none", 0.95)
  expect_equal(fit$estimate, did$estimate, tolerance = 1e-6)
})

test_that("jackknife SE is produced with multiple treated units", {
  df <- sim_panel(N = 20, T = 11, n_treated = 4, t0 = 8, seed = 9)
  fit <- trop(df, "y", "w", "id", "t", se = "jackknife",
              control = trop_control(n_cv_cells = 20, seed = 1, cv_cycles = 1))
  expect_true(is.finite(fit$std.error))
  expect_lt(fit$conf.low, fit$estimate)
  expect_gt(fit$conf.high, fit$estimate)
})

test_that("assignment pattern detection distinguishes block vs general", {
  df <- sim_panel(N = 15, T = 10, n_treated = 3, t0 = 7, seed = 2)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  expect_equal(cfcompare:::.assignment_pattern(m$W)$type, "block")
  W2 <- m$W
  W2[1, 9] <- 0  # turn treatment off -> non-absorbing
  expect_equal(cfcompare:::.assignment_pattern(W2)$type, "general")
})

test_that("input validation catches malformed panels", {
  df <- sim_panel(N = 12, T = 8, n_treated = 2, t0 = 6, seed = 1)
  expect_error(trop(df, "nope", "w", "id", "t"), "not found")
  df$w[1] <- 2
  expect_error(trop(df, "y", "w", "id", "t"), "0/1")
})
