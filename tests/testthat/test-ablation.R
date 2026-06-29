make_ab <- function(se = "none") {
  df <- sim_panel(N = 24, T = 10, n_treated = 4, t0 = 6, att = 2, seed = 1)
  trop_ablation(df, "y", "w", "id", "t", anchor = "pooled", se = se,
                control = trop_control(n_cv_cells = 12L, cv_cycles = 1L, seed = 1))
}

test_that("trop_ablation returns the six-spec table", {
  ab <- make_ab()
  expect_s3_class(ab, "trop_ablation")
  expect_equal(nrow(ab), 6L)
  expect_true(all(c("spec", "lambda_time", "lambda_unit", "lambda_nn",
                    "estimate", "rank") %in% names(ab)))
  # DID-like and no-regression rows force the nuclear-norm penalty to Inf
  expect_true(any(is.infinite(ab$lambda_nn)))
})

test_that("format() yields markdown and latex", {
  ab <- make_ab()
  md <- format(ab, "markdown")
  expect_true(any(grepl("Specification", md)))
  tex <- format(ab, "latex")
  expect_true(any(grepl("\\\\toprule", tex)))
})

test_that("plot() writes a PNG figure", {
  ab <- make_ab()
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)
  res <- plot(ab, file = f)
  expect_s3_class(res, "trop_ablation")   # returns the object invisibly
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 1000)
})

test_that("plot() writes a PDF and appends .png when no extension is given", {
  ab <- make_ab("jackknife")             # exercises the SE/CI columns too
  fp <- tempfile(fileext = ".pdf")
  fn <- tempfile()                        # no extension
  on.exit(unlink(c(fp, paste0(fn, ".png"))), add = TRUE)
  plot(ab, file = fp)
  plot(ab, file = fn)
  expect_true(file.exists(fp))
  expect_true(file.exists(paste0(fn, ".png")))
})

test_that("plot() renders to the active device and rejects other formats", {
  ab <- make_ab()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s3_class(plot(ab), "trop_ablation")
  expect_error(plot(ab, file = tempfile(fileext = ".svg")), "png or .pdf")
})
