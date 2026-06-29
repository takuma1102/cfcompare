fit_for <- function(att, seed) {
  df <- sim_panel(N = 36, T = 12, n_treated = 5, t0 = 9, rank = 3, att = att,
                  noise = 1, seed = seed)
  trop(df, "y", "w", "id", "t", se = "none",
       control = trop_control(n_cv_cells = 10L, cv_cycles = 1L, seed = 1))
}

test_that("trop_placebo_test returns the expected structure", {
  pt <- trop_placebo_test(fit_for(3, 1), B = 60,
                          control = trop_control(n_cv_cells = 10L, cv_cycles = 1L,
                                                 seed = 1))
  expect_s3_class(pt, "trop_placebo_test")
  expect_named(pt, c("observed", "placebo", "p.value", "mean", "sd",
                     "null.low", "null.high", "alternative", "B",
                     "n_treated_units", "outcome"))
  expect_true(pt$p.value >= 0 && pt$p.value <= 1)
  expect_equal(pt$n_treated_units, 5L)
  expect_length(pt$placebo, pt$B)
})

test_that("a clear effect is far from the placebo null", {
  pt <- trop_placebo_test(fit_for(4, 2), B = 100,
                          control = trop_control(n_cv_cells = 10L, cv_cycles = 1L,
                                                 seed = 1))
  expect_lt(pt$p.value, 0.05)
  expect_gt(abs(pt$observed), abs(pt$mean) + 2 * pt$sd)
})

test_that("alternatives run and the null distribution is centred near zero", {
  pt <- trop_placebo_test(fit_for(3, 3), B = 80, alternative = "greater",
                          control = trop_control(n_cv_cells = 10L, cv_cycles = 1L,
                                                 seed = 1))
  expect_identical(pt$alternative, "greater")
  expect_lt(abs(pt$mean), 0.5)   # placebo ATT averages near 0 by construction
})

test_that("trop_placebo_test validates its input", {
  expect_error(trop_placebo_test(list(), B = 10), "trop")
  # too few controls for the number of treated units
  df <- sim_panel(N = 6, T = 10, n_treated = 5, t0 = 7, att = 2, seed = 9)
  fit <- trop(df, "y", "w", "id", "t", se = "none",
              control = trop_control(n_cv_cells = 6L, cv_cycles = 1L))
  expect_error(trop_placebo_test(fit, B = 10), "control")
})

test_that("autoplot.trop_placebo_test returns a ggplot", {
  skip_if_not_installed("ggplot2")
  pt <- trop_placebo_test(fit_for(3, 4), B = 40,
                          control = trop_control(n_cv_cells = 10L, cv_cycles = 1L,
                                                 seed = 1))
  expect_s3_class(ggplot2::autoplot(pt), "ggplot")
})
