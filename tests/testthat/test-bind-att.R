schema_cols <- c("method", "estimate", "std.error", "conf.low", "conf.high",
                 "n_treated_cells", "n_treated_units", "outcome", "engine",
                 "rank", "note")

make_fit <- function(seed) {
  df <- sim_panel(N = 18, T = 10, n_treated = 3, t0 = 8, att = 2, seed = seed)
  trop(df, "y", "w", "id", "t", se = "none",
       control = trop_control(n_cv_cells = 12L, cv_cycles = 1L, seed = 1))
}

test_that("bind_att stacks fits and relabels single-row results by name", {
  f1 <- make_fit(1); f2 <- make_fit(2)
  tbl <- bind_att(pooled = f1, per_cell = f2)
  expect_s3_class(tbl, "cf_att_tbl")
  expect_equal(nrow(tbl), 2L)
  expect_setequal(tbl$method, c("pooled", "per_cell"))
  expect_true(all(schema_cols %in% names(tbl)))
})

test_that("bind_att keeps every row of a multi-row input", {
  f1 <- make_fit(3)
  df <- sim_panel(N = 18, T = 10, n_treated = 3, t0 = 8, att = 2, seed = 4)
  cmp <- panel_compare(df, "y", "w", "id", "t",
                       methods = c("DID", "TROP"), se = "none",
                       control = trop_control(n_cv_cells = 12L, cv_cycles = 1L))
  tbl <- bind_att(cmp, extra = f1)
  expect_equal(nrow(tbl), 3L)
  expect_true(all(c("DID", "TROP", "extra") %in% tbl$method))
})

test_that("bind_att warns when a name is given to a multi-row input", {
  df <- sim_panel(N = 16, T = 10, n_treated = 3, t0 = 8, att = 2, seed = 5)
  cmp <- panel_compare(df, "y", "w", "id", "t",
                       methods = c("DID", "TROP"), se = "none",
                       control = trop_control(n_cv_cells = 12L, cv_cycles = 1L))
  expect_warning(bind_att(both = cmp), "multi-row")
})

test_that("bind_att errors with no arguments", {
  expect_error(bind_att(), "at least one")
})

test_that("as_att.list honours list names for single-row parts", {
  f1 <- make_fit(6); f2 <- make_fit(7)
  tbl <- as_att(list(a = f1, b = f2))
  expect_equal(nrow(tbl), 2L)
  expect_setequal(tbl$method, c("a", "b"))
})

test_that(".rbind_att aligns differing columns and drops empty parts", {
  f1 <- make_fit(8)
  row <- as_att(f1, method = "TROP")
  # a part missing a schema column should be filled with NA, not error
  thin <- data.frame(method = "ext", estimate = 0.5, stringsAsFactors = FALSE)
  tbl <- bind_att(row, ext = thin)
  expect_equal(nrow(tbl), 2L)
  expect_true(all(schema_cols %in% names(tbl)))
  expect_true(is.na(tbl$std.error[tbl$method == "ext"]))
})
