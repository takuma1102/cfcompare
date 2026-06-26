test_that("trop_sensitivity returns a labelled grid with a CV-selected cell", {
  df <- sim_panel(N = 22, T = 12, n_treated = 3, t0 = 10, att = 2, seed = 1)
  g <- trop_sensitivity(df, "y", "w", "id", "t",
                        lambda_time = c(0, 0.25, 1), lambda_nn = NULL,
                        control = trop_control(n_cv_cells = 20L, cv_cycles = 1L),
                        seed = 1)
  expect_s3_class(g, "cf_trop_grid")
  expect_true(all(c("lambda_time", "lambda_nn", "lambda_unit", "att", "cv_loss")
                  %in% names(g)))
  sel <- attr(g, "selected")
  expect_equal(nrow(sel), 1L)
  expect_equal(sel$cv_loss, min(g$cv_loss))
})
