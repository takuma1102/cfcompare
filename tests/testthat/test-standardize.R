test_that("standardize = TRUE leaves the TWFE/DID estimate unchanged", {
  # With lambda = (0, 0, Inf) the fit is plain weighted TWFE with uniform
  # weights; the ATT is invariant to an affine transform of the outcome, so the
  # standardized and raw fits must agree.
  df <- sim_panel(N = 18, T = 12, n_treated = 3, t0 = 9, att = 2, seed = 21)
  lam <- list(time = 0, unit = 0, nn = Inf)
  raw <- trop(df, "y", "w", "id", "t", lambda = lam, se = "none")
  std <- trop(df, "y", "w", "id", "t", lambda = lam, se = "none",
              standardize = TRUE)
  # small solver-tolerance differences remain because the relative
  # stopping rule acts on different scales
  expect_equal(std$estimate, raw$estimate, tolerance = 1e-4)
  # per-cell effects and the counterfactual are mapped back to the raw scale
  expect_equal(std$tau_cells$tau, raw$tau_cells$tau, tolerance = 1e-5)
  expect_equal(std$tau_cells$y, raw$tau_cells$y, tolerance = 1e-8)
  fin <- is.finite(raw$counterfactual)   # per-cell anchor leaves NAs off treated cells
  expect_equal(std$counterfactual[fin], raw$counterfactual[fin], tolerance = 1e-3)
  # scaling metadata is stored and panel$Y is on the fitting scale
  expect_true(abs(std$scaling$scale - stats::sd(raw$panel$Y)) < 1e-8)
  expect_lt(abs(mean(std$panel$Y)), 1e-8)
  expect_lt(abs(stats::sd(std$panel$Y) - 1), 1e-8)
})

test_that("standardize maps counterfactual_matrix back to the raw scale", {
  df <- sim_panel(N = 16, T = 10, n_treated = 3, t0 = 8, att = 1, seed = 22)
  lam <- list(time = 0, unit = 0, nn = Inf)
  raw <- trop(df, "y", "w", "id", "t", lambda = lam, se = "none")
  std <- trop(df, "y", "w", "id", "t", lambda = lam, se = "none",
              standardize = TRUE)
  expect_equal(counterfactual_matrix(std), counterfactual_matrix(raw),
               tolerance = 1e-4)
})

test_that("standardized event-study effects are on the raw outcome scale", {
  df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 3, seed = 23)
  lam <- list(time = 0, unit = 0, nn = Inf)
  std <- trop(df, "y", "w", "id", "t", lambda = lam, se = "none",
              standardize = TRUE)
  es <- trop_event_study(std, se = "none")
  post <- es$estimates[es$estimates$period == "post", ]
  # post-period effects should average near the (raw-scale) ATT
  expect_equal(mean(post$estimate), std$estimate, tolerance = 0.25)
})

test_that("n_cv_cells = Inf / 0 uses every control cell", {
  df <- sim_panel(N = 8, T = 6, n_treated = 2, t0 = 5, seed = 24)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  all_cells <- which(m$W == 0, arr.ind = TRUE)
  expect_equal(nrow(cfcompare:::.sample_control_cells(m$W, Inf)),
               nrow(all_cells))
  expect_equal(nrow(cfcompare:::.sample_control_cells(m$W, 0L)),
               nrow(all_cells))
  expect_lte(nrow(cfcompare:::.sample_control_cells(m$W, 5L, seed = 1)), 5L)
})
