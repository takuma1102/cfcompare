schema_cols <- c("method", "estimate", "std.error", "conf.low", "conf.high",
                 "n_treated_cells", "n_treated_units", "outcome", "engine",
                 "rank", "note")

test_that("panel_compare returns the shared tidy schema", {
  skip_on_cran()
  df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 4)
  cmp <- panel_compare(df, "y", "w", "id", "t",
                       methods = c("DID", "MC", "TROP"), se = "none",
                       control = trop_control(n_cv_cells = 20, seed = 1,
                                              cv_cycles = 1))
  expect_s3_class(cmp, "cf_comparison")
  expect_s3_class(cmp$att, "cf_att_tbl")
  expect_setequal(names(cmp$att), schema_cols)
  expect_setequal(cmp$att$method, c("DID", "MC", "TROP"))
  expect_true(all(is.finite(cmp$att$estimate)))
})

test_that("SDID is skipped gracefully when synthdid is unavailable", {
  skip_on_cran()
  skip_if(requireNamespace("synthdid", quietly = TRUE),
          "synthdid is installed; skip-path not exercised")
  df <- sim_panel(N = 16, T = 10, n_treated = 3, t0 = 8, seed = 6)
  expect_message(
    cmp <- panel_compare(df, "y", "w", "id", "t",
                         methods = c("DID", "SDID"), se = "none",
                         control = trop_control(n_cv_cells = 15, seed = 1,
                                                cv_cycles = 1)),
    "skipping SDID"
  )
  sdid_row <- cmp$att[cmp$att$method == "SDID", ]
  expect_true(is.na(sdid_row$estimate))
  expect_match(sdid_row$note, "synthdid")
})

test_that("unknown methods raise an informative error", {
  skip_on_cran()
  df <- sim_panel(N = 12, T = 8, n_treated = 2, t0 = 6, seed = 1)
  expect_error(panel_compare(df, "y", "w", "id", "t", methods = "FOO"),
               "Unknown method")
})

test_that("as_att coerces trop fits and lists", {
  skip_on_cran()
  df <- sim_panel(N = 16, T = 10, n_treated = 3, t0 = 8, seed = 8)
  fit <- trop(df, "y", "w", "id", "t", se = "none",
              control = trop_control(n_cv_cells = 15, seed = 1, cv_cycles = 1))
  row <- as_att(fit, method = "TROP")
  expect_s3_class(row, "cf_att_tbl")
  expect_equal(row$estimate, fit$estimate)

  combined <- as_att(list(fit, fit))
  expect_equal(nrow(combined), 2L)
})
