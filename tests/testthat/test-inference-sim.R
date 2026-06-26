test_that("bootstrap SE returns finite SE and a percentile CI", {
  df <- sim_panel(N = 24, T = 12, n_treated = 4, t0 = 10, att = 2, seed = 1)
  fit <- trop(df, "y", "w", "id", "t", se = "bootstrap",
              control = trop_control(n_cv_cells = 30L, cv_cycles = 1L,
                                     n_boot = 30L, seed = 1))
  expect_equal(fit$se.method, "bootstrap")
  expect_true(is.finite(fit$std.error) && fit$std.error > 0)
  expect_true(fit$conf.low < fit$estimate && fit$estimate < fit$conf.high)
})

test_that("sim_semisynthetic imposes the requested effect exactly", {
  real <- sim_panel(N = 30, T = 14, n_treated = 1L, att = 0, seed = 1)
  ss <- sim_semisynthetic(real, "y", "id", "t", n_treated = 5, t0 = 11,
                          att = 4, seed = 2)
  expect_true(all(c("id", "t", "y", "w", "y0", "tau") %in% names(ss)))
  expect_equal(mean(ss$tau[ss$w == 1]), 4)
  expect_equal(ss$y[ss$w == 0], ss$y0[ss$w == 0])
})

test_that("sim_semisynthetic supports a dynamic effect path", {
  real <- sim_panel(N = 20, T = 12, n_treated = 1L, att = 0, seed = 3)
  ss <- sim_semisynthetic(real, "y", "id", "t", n_treated = 4, t0 = 10,
                          effect = c(1, 2, 3), seed = 1)
  by_t <- tapply(ss$tau[ss$w == 1], ss$t[ss$w == 1], mean)
  expect_equal(as.numeric(by_t), c(1, 2, 3))
})
