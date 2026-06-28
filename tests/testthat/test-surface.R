test_that("surface_matrix lays the grid out over the two swept penalties", {
  df <- sim_panel(N = 22, T = 12, n_treated = 3, t0 = 10, att = 2, seed = 1)
  g <- trop_sensitivity(df, "y", "w", "id", "t",
                        axes = c("nn", "time"),
                        lambda_time = c(0, 0.25, 1), lambda_nn = NULL,
                        control = trop_control(n_cv_cells = 15L, cv_cycles = 1L),
                        seed = 1)
  m <- surface_matrix(g, "cv_loss")
  ax <- attr(g, "axes")
  xcol <- paste0("lambda_", unname(ax["x"]))
  ycol <- paste0("lambda_", unname(ax["y"]))
  expect_true(is.matrix(m))
  expect_equal(nrow(m), length(unique(g[[ycol]])))
  expect_equal(ncol(m), length(unique(g[[xcol]])))
  expect_identical(names(dimnames(m)), c(ycol, xcol))

  # a sampled cell matches the corresponding long-grid row
  yv <- sort(unique(g[[ycol]])); xv <- sort(unique(g[[xcol]]))
  hit <- g[[ycol]] == yv[1] & g[[xcol]] == xv[1]
  expect_equal(unname(m[1, 1]), g$cv_loss[hit])
})

test_that("surface_matrix can also lay out the ATT", {
  df <- sim_panel(N = 20, T = 12, n_treated = 4, t0 = 9, att = 2, seed = 2)
  g <- trop_sensitivity(df, "y", "w", "id", "t",
                        lambda_time = c(0, 0.5), lambda_nn = NULL,
                        control = trop_control(n_cv_cells = 12L, cv_cycles = 1L),
                        seed = 1)
  m <- surface_matrix(g, "att")
  expect_true(is.matrix(m))
  expect_equal(length(m), nrow(g))
})

test_that("surface_matrix rejects non-grid input", {
  expect_error(surface_matrix(data.frame(a = 1)), "cf_trop_grid")
})

test_that("selected_lambda returns a trop()-ready penalty list", {
  df <- sim_panel(N = 22, T = 12, n_treated = 3, t0 = 10, att = 2, seed = 3)
  g <- trop_sensitivity(df, "y", "w", "id", "t",
                        lambda_time = c(0, 0.25, 1), lambda_nn = NULL,
                        control = trop_control(n_cv_cells = 15L, cv_cycles = 1L),
                        seed = 1)
  lam <- selected_lambda(g)
  expect_named(lam, c("time", "unit", "nn"))
  sel <- attr(g, "selected")
  expect_equal(lam$time, sel$lambda_time)
  expect_equal(lam$unit, sel$lambda_unit)
  expect_equal(lam$nn,   sel$lambda_nn)

  # the selected penalties feed straight back into a fit
  fit <- trop(df, "y", "w", "id", "t", lambda = lam, se = "none",
              control = trop_control(n_cv_cells = 10L, cv_cycles = 1L))
  expect_true(is.finite(fit$estimate))
})
