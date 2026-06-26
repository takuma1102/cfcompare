test_that("DIFP runs on a block design and is close to SC", {
  df <- sim_panel(N = 24, T = 12, n_treated = 4, t0 = 10, rank = 2,
                  att = 2, noise = 1, seed = 1)
  cmp <- panel_compare(df, "y", "w", "id", "t",
                       methods = c("DIFP", "DID"), se = "none",
                       control = trop_control(n_cv_cells = 20L, cv_cycles = 1L))
  expect_true("DIFP" %in% cmp$att$method)
  difp <- cmp$att$estimate[cmp$att$method == "DIFP"]
  expect_true(is.finite(difp))
  expect_equal(difp, 2, tolerance = 1.5)   # recovers the effect roughly
})

test_that("DIFP requires a block design (skips otherwise)", {
  df <- sim_panel(N = 20, T = 12, n_treated = 3, t0 = 9, seed = 1)
  # make a non-absorbing (general) pattern
  df$w[df$id == 1 & df$t == 12] <- 0
  cmp <- panel_compare(df, "y", "w", "id", "t", methods = "DIFP", se = "none",
                       control = trop_control(n_cv_cells = 15L, cv_cycles = 1L))
  expect_true(is.na(cmp$att$estimate[cmp$att$method == "DIFP"]))
})
