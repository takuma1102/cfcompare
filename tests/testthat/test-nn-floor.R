# The solver adds a small, scale-relative stabilising floor to every finite
# lambda_nn (trop_control(nn_floor = "auto", nn_floor_scale = 1e-3), resolved
# as nn_floor_scale * s1 with s1 the operator norm of the two-way-demeaned
# control outcome on the fitting scale). lambda_nn = Inf is never affected, and
# nn_floor = 0 disables the adjustment so penalties apply exactly as given.
# These tests lock that behaviour.

test_that("trop_control() validates nn_floor and nn_floor_scale", {
  expect_error(trop_control(nn_floor = -1), "nn_floor")
  expect_error(trop_control(nn_floor = c(1, 2)), "nn_floor")
  expect_error(trop_control(nn_floor = "big"), "nn_floor")
  expect_error(trop_control(nn_floor = NA_real_), "nn_floor")
  expect_error(trop_control(nn_floor_scale = -0.1), "nn_floor_scale")
  ctrl <- trop_control()
  expect_identical(ctrl$nn_floor, "auto")
  expect_equal(ctrl$nn_floor_scale, 1e-3)
  expect_equal(trop_control(nn_floor = 0)$nn_floor, 0)
})

test_that("auto floor resolves to nn_floor_scale * operator norm and is stored", {
  df <- sim_panel(N = 18, T = 10, n_treated = 3, t0 = 8, att = 2, seed = 31)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  s1 <- cfcompare:::.trop_nn_scale(m$Y, m$W)
  fit <- trop(df, "y", "w", "id", "t",
              lambda = list(time = 0, unit = 0, nn = 2), se = "none")
  expect_gt(fit$lambda$nn_floor, 0)
  expect_equal(fit$lambda$nn_floor, 1e-3 * s1, tolerance = 1e-10)
  expect_equal(fit$lambda$nn, 2)          # reported penalty stays nominal
  # nn_floor_scale rescales the auto value
  f2 <- trop(df, "y", "w", "id", "t",
             lambda = list(time = 0, unit = 0, nn = 2), se = "none",
             control = trop_control(nn_floor_scale = 1e-2))
  expect_equal(f2$lambda$nn_floor, 10 * fit$lambda$nn_floor, tolerance = 1e-10)
  # a numeric nn_floor is used as-is
  f3 <- trop(df, "y", "w", "id", "t",
             lambda = list(time = 0, unit = 0, nn = 2), se = "none",
             control = trop_control(nn_floor = 0.25))
  expect_equal(f3$lambda$nn_floor, 0.25)
})

test_that("the floor is additive: (nn = v, floor = F) == (nn = v + F, floor = 0)", {
  df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 32)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  s1 <- cfcompare:::.trop_nn_scale(m$Y, m$W)
  v <- 0.1 * s1; fl <- 0.05 * s1
  a <- trop(df, "y", "w", "id", "t",
            lambda = list(time = 0.1, unit = 0.2, nn = v), se = "none",
            control = trop_control(nn_floor = fl))
  b <- trop(df, "y", "w", "id", "t",
            lambda = list(time = 0.1, unit = 0.2, nn = v + fl), se = "none",
            control = trop_control(nn_floor = 0))
  expect_equal(a$estimate, b$estimate)
  expect_equal(a$counterfactual[is.finite(a$counterfactual)],
               b$counterfactual[is.finite(b$counterfactual)])
  expect_equal(a$lambda$nn_floor, fl)
  expect_equal(b$lambda$nn_floor, 0)
})

test_that("lambda_nn = Inf (the DID/TWFE special case) is never affected", {
  df <- sim_panel(N = 16, T = 10, n_treated = 3, t0 = 8, att = 2, seed = 33)
  lamI <- list(time = 0, unit = 0, nn = Inf)
  a <- trop(df, "y", "w", "id", "t", lambda = lamI, se = "none")   # auto floor
  b <- trop(df, "y", "w", "id", "t", lambda = lamI, se = "none",
            control = trop_control(nn_floor = 0))
  expect_equal(a$estimate, b$estimate)
  expect_equal(a$rank, 0L)
})

test_that("lambda_nn = 0 is well posed under the default floor", {
  df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 34)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  s1 <- cfcompare:::.trop_nn_scale(m$Y, m$W)
  z <- trop(df, "y", "w", "id", "t",
            lambda = list(time = 0, unit = 0, nn = 0), se = "none")
  expect_true(is.finite(z$estimate))
  # ... and is identical to solving at the resolved floor with the adjustment
  # disabled (the additive convention at the zero limit).
  ref <- trop(df, "y", "w", "id", "t",
              lambda = list(time = 0, unit = 0, nn = 1e-3 * s1), se = "none",
              control = trop_control(nn_floor = 0))
  expect_equal(z$estimate, ref$estimate)
})

test_that("the auto floor is scale-relative: standardized = raw / sd(Y)", {
  df <- sim_panel(N = 18, T = 10, n_treated = 3, t0 = 8, att = 2, seed = 35)
  lam <- list(time = 0, unit = 0, nn = 1)
  raw <- trop(df, "y", "w", "id", "t", lambda = lam, se = "none")
  std <- trop(df, "y", "w", "id", "t", lambda = lam, se = "none",
              standardize = TRUE)
  expect_equal(std$lambda$nn_floor,
               raw$lambda$nn_floor / std$scaling$scale, tolerance = 1e-8)
})

test_that("trop_matrix() shares the convention; nn_floor = 0 restores exactness", {
  df <- sim_panel(N = 15, T = 12, n_treated = 3, t0 = 9, att = 2, seed = 37)
  pm <- panel_matrices(df, "y", "w", "id", "t")
  s1 <- cfcompare:::.trop_nn_scale(pm$Y, pm$W)
  v <- 0.2 * s1
  withf <- trop_matrix(pm$Y, pm$W, lambda_unit = 0.1, lambda_time = 0.1,
                       lambda_nn = v)
  exact <- trop_matrix(pm$Y, pm$W, lambda_unit = 0.1, lambda_time = 0.1,
                       lambda_nn = v + 1e-3 * s1,
                       control = trop_control(nn_floor = 0))
  expect_equal(withf, exact)
  # the Inf special case (weighted TWFE, the exact-agreement check) is unchanged
  expect_equal(
    trop_matrix(pm$Y, pm$W, lambda_unit = 0.1, lambda_time = 0.1,
                lambda_nn = Inf),
    trop_matrix(pm$Y, pm$W, lambda_unit = 0.1, lambda_time = 0.1,
                lambda_nn = Inf, control = trop_control(nn_floor = 0)))
})

test_that("downstream consumers reuse the fit's resolved floor via $lambda", {
  df <- sim_panel(N = 18, T = 10, n_treated = 3, t0 = 8, att = 2, seed = 38)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  s1 <- cfcompare:::.trop_nn_scale(m$Y, m$W)
  v <- 0.1 * s1; fl <- 0.05 * s1
  a <- trop(df, "y", "w", "id", "t",
            lambda = list(time = 0, unit = 0, nn = v), se = "none",
            control = trop_control(nn_floor = fl))
  b <- trop(df, "y", "w", "id", "t",
            lambda = list(time = 0, unit = 0, nn = v + fl), se = "none",
            control = trop_control(nn_floor = 0))
  # counterfactual_matrix() re-solves with a fresh trop_control(): the floor
  # stored in a$lambda must make it match b's exact-penalty counterfactual.
  expect_equal(counterfactual_matrix(a), counterfactual_matrix(b))
  ea <- trop_event_study(a, se = "none")
  eb <- trop_event_study(b, se = "none")
  expect_equal(ea$estimates$estimate, eb$estimates$estimate)
})
