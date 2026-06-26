test_that("panel_rmse returns a ranked tidy table for native methods", {
  df <- sim_panel(N = 24, T = 12, n_treated = 3, t0 = 10, rank = 2,
                  att = 2, noise = 1, seed = 2)
  r <- panel_rmse(df, "y", "w", "id", "t",
                  methods = c("DID", "MC", "TROP"),
                  horizon = 2, n_pseudo = 4, n_runs = 3,
                  control = trop_control(n_cv_cells = 30L, cv_cycles = 1L),
                  seed = 1)
  expect_s3_class(r, "cf_rmse_tbl")
  expect_true(all(c("method", "rmse", "rmse_se", "n_runs", "engine") %in% names(r)))
  expect_equal(nrow(r), 3L)
  expect_true(all(is.finite(r$rmse)))
  expect_true(all(r$rmse > 0))
})

test_that("panel_rmse skips wrapped methods gracefully when unavailable", {
  df <- sim_panel(N = 20, T = 12, n_treated = 3, t0 = 10, seed = 3)
  r <- panel_rmse(df, "y", "w", "id", "t",
                  methods = c("DID", "SDID"),
                  horizon = 2, n_pseudo = 4, n_runs = 2,
                  control = trop_control(n_cv_cells = 30L, cv_cycles = 1L),
                  seed = 1)
  expect_true(is.na(r$rmse[r$method == "SDID"]))
  expect_false(is.na(r$note[r$method == "SDID"]))
  expect_true(is.finite(r$rmse[r$method == "DID"]))
})

test_that("panel_rmse validates horizon and pseudo count", {
  df <- sim_panel(N = 12, T = 10, n_treated = 2, t0 = 8, seed = 1)
  expect_error(panel_rmse(df, "y", "w", "id", "t", horizon = 10))
  expect_error(panel_rmse(df, "y", "w", "id", "t", n_pseudo = 50,
                          methods = "DID"))
})
