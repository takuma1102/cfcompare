schema_cols <- c("method", "estimate", "std.error", "conf.low", "conf.high",
                 "n_treated_cells", "n_treated_units", "outcome", "engine",
                 "rank", "note")

test_that("compare_se_modes shares one estimate across inference methods", {
  df <- sim_panel(N = 24, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 1)
  tbl <- compare_se_modes(df, "y", "w", "id", "t",
                          se = c("jackknife", "placebo"),
                          control = trop_control(n_cv_cells = 12L, cv_cycles = 1L,
                                                 seed = 1))
  expect_s3_class(tbl, "cf_att_tbl")
  expect_equal(nrow(tbl), 2L)
  expect_true(all(schema_cols %in% names(tbl)))
  expect_setequal(tbl$method, c("TROP_jackknife", "TROP_placebo"))
  # the point estimate is fixed (penalties chosen once), so all rows agree
  expect_equal(length(unique(round(tbl$estimate, 8))), 1L)
})

test_that("compare_se_modes defaults to all three SE methods", {
  df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 2)
  tbl <- compare_se_modes(df, "y", "w", "id", "t",
                          control = trop_control(n_cv_cells = 10L, cv_cycles = 1L,
                                                 n_boot = 8L, seed = 1))
  expect_setequal(tbl$method,
                  c("TROP_bootstrap", "TROP_jackknife", "TROP_placebo"))
})

test_that("compare_se_modes accepts custom labels", {
  df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 3)
  tbl <- compare_se_modes(df, "y", "w", "id", "t",
                          se = c("jackknife", "placebo"),
                          labels = c("JK", "PL"),
                          control = trop_control(n_cv_cells = 10L, cv_cycles = 1L,
                                                 seed = 1))
  expect_setequal(tbl$method, c("JK", "PL"))
})

test_that("compare_se_modes validates labels length and SE names", {
  df <- sim_panel(N = 16, T = 10, n_treated = 3, t0 = 8, att = 2, seed = 4)
  ctrl <- trop_control(n_cv_cells = 10L, cv_cycles = 1L, seed = 1)
  expect_error(
    compare_se_modes(df, "y", "w", "id", "t", se = "jackknife",
                     labels = c("a", "b"), control = ctrl),
    "labels"
  )
  expect_error(
    compare_se_modes(df, "y", "w", "id", "t", se = "nonsense", control = ctrl)
  )
})
