test_that("warm-started nuclear path matches cold-start CV to solver tolerance", {
  df <- sim_panel(N = 16, T = 12, n_treated = 3, t0 = 9, att = 2, seed = 7)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  ctrl <- trop_control(n_cv_cells = 6L, seed = 1)
  cells <- cfcompare:::.sample_control_cells(m$W, 6L, seed = 1)
  du_l <- lapply(seq_len(nrow(cells)), function(k) {
    Wk <- m$W; Wk[cells[k, 1], cells[k, 2]] <- 1
    cfcompare:::.unit_distance_to(m$Y, Wk, cells[k, 1], cells[k, 2])
  })
  nn_grid <- sort(cfcompare:::.trop_default_grids(m$Y, m$W)$nn,
                  decreasing = TRUE)
  lam_base <- list(time = 0.1, unit = 0.3)
  warm <- NULL
  qs_w <- qs_c <- numeric(length(nn_grid))
  for (gi in seq_along(nn_grid)) {
    l2 <- c(lam_base, list(nn = nn_grid[gi]))
    ev <- cfcompare:::.trop_cv_Q(m$Y, m$W, l2, ctrl, cells, du_list = du_l,
                                 warm = warm, keep_state = TRUE)
    qs_w[gi] <- ev$Q
    warm <- ev$warm
    qs_c[gi] <- cfcompare:::.trop_cv_Q(m$Y, m$W, l2, ctrl, cells,
                                       du_list = du_l)$Q
  }
  # convex program: warm starts change iteration counts, not the criterion
  # (beyond solver tolerance) nor the selected penalty
  expect_lt(max(abs(qs_w - qs_c)), 1e-4)
  expect_equal(which.min(qs_w), which.min(qs_c))
})

test_that("full-state warm start reproduces the cold-start solution", {
  df <- sim_panel(N = 14, T = 10, n_treated = 3, t0 = 8, seed = 9)
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  mask <- (m$W == 0) * 1
  w <- matrix(1, nrow(m$Y), ncol(m$Y))
  s1 <- max(svd(cfcompare:::.double_demean(m$Y))$d)
  f1 <- cfcompare:::.mcnnm_fit(m$Y, mask, w, 0.5 * s1)
  st <- cfcompare:::.solver_state(f1)
  f2 <- cfcompare:::.mcnnm_fit(m$Y, mask, w, 0.2 * s1, state_init = st)
  f2c <- cfcompare:::.mcnnm_fit(m$Y, mask, w, 0.2 * s1)
  expect_equal(f2$M, f2c$M, tolerance = 1e-4)
})
