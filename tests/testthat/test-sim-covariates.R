test_that("sim_panel is unchanged when no covariates are requested", {
  skip_on_cran()
  a <- sim_panel(N = 20, T = 10, n_treated = 3, t0 = 8, att = 2, seed = 5)
  b <- sim_panel(N = 20, T = 10, n_treated = 3, t0 = 8, att = 2, seed = 5,
                 n_cov = 0L)
  expect_identical(a, b)
  expect_setequal(names(a), c("id", "t", "y", "w", "y0"))
  expect_null(attr(a, "phi"))
})

test_that("sim_panel adds covariate columns and records the coefficients", {
  skip_on_cran()
  df <- sim_panel(N = 25, T = 12, n_treated = 4, t0 = 9, att = 2,
                  n_cov = 2, phi = c(1.2, -0.7), seed = 1)
  expect_true(all(c("x1", "x2") %in% names(df)))
  ph <- attr(df, "phi")
  expect_equal(unname(ph), c(1.2, -0.7))
  expect_identical(names(ph), c("x1", "x2"))
  expect_false(anyNA(df$x1))
})

test_that("sim_panel draws phi when only n_cov is given", {
  skip_on_cran()
  df <- sim_panel(N = 18, T = 10, n_treated = 3, t0 = 8, n_cov = 3, seed = 2)
  ph <- attr(df, "phi")
  expect_length(ph, 3L)
  expect_true(all(is.finite(ph)))
  expect_true(all(paste0("x", 1:3) %in% names(df)))
})

test_that("sim_panel validates n_cov / phi", {
  skip_on_cran()
  expect_error(
    sim_panel(N = 12, T = 8, n_treated = 2, t0 = 6, n_cov = 2, phi = c(1, 2, 3)),
    "phi"
  )
  expect_error(
    sim_panel(N = 12, T = 8, n_treated = 2, t0 = 6, n_cov = -1),
    "non-negative"
  )
})

test_that("covariate signal sits in Y(0): true ATT equals att", {
  skip_on_cran()
  df <- sim_panel(N = 25, T = 12, n_treated = 4, t0 = 9, att = 2,
                  n_cov = 2, phi = c(1.2, -0.7), seed = 3)
  # y - y0 is exactly att on treated cells regardless of covariates
  expect_equal(true_att(df), 2)
})

test_that("covariate coefficients are recoverable by trop(covariates=)", {
  skip_on_cran()
  df <- sim_panel(N = 30, T = 14, n_treated = 5, t0 = 10, att = 2,
                  n_cov = 2, phi = c(1.2, -0.7), seed = 11)
  fit <- trop(df, "y", "w", "id", "t", covariates = c("x1", "x2"),
              anchor = "pooled", se = "none",
              control = trop_control(n_cv_cells = 15L, cv_cycles = 1L, seed = 1))
  expect_length(fit$phi, 2L)
  expect_true(all(is.finite(fit$phi)))
})

test_that("true_att reads tau for semi-synthetic panels", {
  skip_on_cran()
  real <- sim_panel(N = 30, T = 16, n_treated = 0L, att = 0, seed = 1)
  ss <- sim_semisynthetic(real, "y", "id", "t",
                          n_treated = 5, t0 = 13, att = 3, seed = 2)
  expect_equal(true_att(ss), 3)
})

test_that("true_att errors clearly on malformed input", {
  skip_on_cran()
  expect_error(true_att(data.frame(a = 1)), "Treatment column")
  expect_error(true_att(data.frame(w = c(0, 1))), "true ATT")
  no_treated <- data.frame(w = c(0, 0), y = c(1, 2), y0 = c(1, 2))
  expect_error(true_att(no_treated), "No treated cells")
})
