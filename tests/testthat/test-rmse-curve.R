test_that("rmse_curve returns a tidy cf_rmse_curve over the swept dimension", {
  skip_on_cran()
  cc <- rmse_curve("n_control", values = c(20, 30), n_runs = 2,
                   methods = c("DID", "TROP"), n_pre = 8, n_post = 4,
                   control = trop_control(n_cv_cells = 10L, cv_cycles = 1L),
                   seed = 1)
  expect_s3_class(cc, "cf_rmse_curve")
  expect_true(all(c("method", "x", "rmse", "bias", "n_runs") %in% names(cc)))
  expect_setequal(unique(cc$x), c(20, 30))
  expect_setequal(unique(cc$method), c("DID", "TROP"))
  expect_identical(attr(cc, "vary"), "n_control")
})

test_that("rmse_curve can sweep pre-treatment periods", {
  skip_on_cran()
  cp <- rmse_curve("n_pre", values = c(6, 10), n_runs = 2,
                   methods = c("DID", "TROP"), n_control = 24,
                   control = trop_control(n_cv_cells = 10L, cv_cycles = 1L),
                   seed = 1)
  expect_s3_class(cp, "cf_rmse_curve")
  expect_setequal(unique(cp$x), c(6, 10))
  expect_identical(attr(cp, "vary"), "n_pre")
})

test_that("rmse_curves bundles both sweeps", {
  skip_on_cran()
  g <- rmse_curves(values_control = c(20, 28), values_pre = c(6, 10),
                   n_runs = 2, methods = c("DID", "TROP"),
                   n_control = 24, n_pre = 8, n_post = 4,
                   control = trop_control(n_cv_cells = 10L, cv_cycles = 1L),
                   seed = 1)
  expect_s3_class(g, "cf_rmse_curves")
  expect_s3_class(g$n_control, "cf_rmse_curve")
  expect_s3_class(g$n_pre, "cf_rmse_curve")
})

test_that("rmse_curve can sweep the number of treated units", {
  skip_on_cran()
  cc <- rmse_curve("n_treated", values = c(2, 4), n_runs = 2,
                   methods = c("DID", "TROP"), n_control = 24,
                   n_pre = 8, n_post = 4,
                   control = trop_control(n_cv_cells = 10L, cv_cycles = 1L),
                   seed = 1)
  expect_s3_class(cc, "cf_rmse_curve")
  expect_setequal(unique(cc$x), c(2, 4))
  expect_identical(attr(cc, "vary"), "n_treated")
  expect_match(attr(cc, "xlab"), "treated units")
})

test_that("rmse_curve can sweep the number of post-treatment periods", {
  skip_on_cran()
  cc <- rmse_curve("n_post", values = c(3, 5), n_runs = 2,
                   methods = c("DID", "TROP"), n_control = 24,
                   n_treated = 4, n_pre = 8,
                   control = trop_control(n_cv_cells = 10L, cv_cycles = 1L),
                   seed = 1)
  expect_s3_class(cc, "cf_rmse_curve")
  expect_setequal(unique(cc$x), c(3, 5))
  expect_identical(attr(cc, "vary"), "n_post")
  expect_match(attr(cc, "xlab"), "post-treatment periods")
})
