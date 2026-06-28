test_that("panel_matrices returns matrices plus treated rows/periods", {
  df <- sim_panel(N = 15, T = 12, n_treated = 3, t0 = 9, att = 2, seed = 1)
  pm <- panel_matrices(df, "y", "w", "id", "t")
  expect_named(pm, c("Y", "W", "units", "times",
                     "treated_units", "treated_periods"))
  expect_equal(dim(pm$Y), c(max(df$id), max(df$t)))
  expect_identical(dim(pm$Y), dim(pm$W))
  expect_equal(pm$treated_units, which(rowSums(pm$W) > 0))
  # t0 = 9, T = 12 -> 4 post periods
  expect_equal(pm$treated_periods, 4L)
})

test_that("trop_matrix infers treated_units / treated_periods from W", {
  df <- sim_panel(N = 15, T = 12, n_treated = 3, t0 = 9, att = 2, seed = 2)
  pm <- panel_matrices(df, "y", "w", "id", "t")
  a_explicit <- trop_matrix(pm$Y, pm$W, pm$treated_units, 0, 0, Inf,
                            pm$treated_periods)
  a_inferred <- trop_matrix(pm$Y, pm$W, lambda_unit = 0, lambda_time = 0,
                            lambda_nn = Inf)
  expect_equal(a_inferred, a_explicit)
})

test_that("trop_matrix still validates an empty treatment matrix", {
  df <- sim_panel(N = 12, T = 10, n_treated = 2, t0 = 8, seed = 3)
  pm <- panel_matrices(df, "y", "w", "id", "t")
  expect_error(
    trop_matrix(pm$Y, pm$W * 0, lambda_unit = 0, lambda_time = 0, lambda_nn = Inf),
    "no treated cells"
  )
})

test_that("panel_matrices errors on degenerate designs", {
  df <- sim_panel(N = 10, T = 8, n_treated = 0L, att = 0, seed = 4)
  expect_error(panel_matrices(df, "y", "w", "id", "t"), "treated")
})
