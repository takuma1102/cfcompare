test_that("CS method is recognised and skips gracefully without did", {
  skip_if(requireNamespace("did", quietly = TRUE))
  df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 10, seed = 1)
  cmp <- panel_compare(df, "y", "w", "id", "t", methods = c("DID", "CS"),
                       se = "none",
                       control = trop_control(n_cv_cells = 20L, cv_cycles = 1L))
  expect_true("CS" %in% cmp$att$method)
  expect_true(is.na(cmp$att$estimate[cmp$att$method == "CS"]))
  expect_false(is.na(cmp$att$note[cmp$att$method == "CS"]))
})

test_that("panel_rmse placebo metric scores native methods and tags engines", {
  df <- sim_panel(N = 24, T = 12, n_treated = 3, t0 = 10, rank = 2,
                  att = 2, noise = 1, seed = 2)
  r <- panel_rmse(df, "y", "w", "id", "t",
                  methods = c("DID", "MC", "TROP"), metric = "placebo",
                  horizon = 2, n_pseudo = 3, n_runs = 2,
                  control = trop_control(n_cv_cells = 20L, cv_cycles = 1L),
                  seed = 1)
  expect_s3_class(r, "cf_rmse_tbl")
  expect_equal(attr(r, "metric"), "placebo")
  expect_true(all(is.finite(r$rmse)))
  expect_true(all(r$engine == "cfcompare"))
})

test_that("unknown methods are rejected by panel_rmse", {
  df <- sim_panel(N = 12, T = 10, n_treated = 2, t0 = 8, seed = 1)
  expect_error(panel_rmse(df, "y", "w", "id", "t", methods = "NOPE"))
})
