test_that("trop_matrix uniform weights equals plain TWFE residual ATT", {
  df <- sim_panel(N = 15, T = 12, n_treated = 3, t0 = 9, att = 2, seed = 3)
  N <- max(df$id); T <- max(df$t)
  Y <- matrix(0, N, T); W <- matrix(0, N, T)
  for (k in seq_len(nrow(df))) { Y[df$id[k], df$t[k]] <- df$y[k]
    W[df$id[k], df$t[k]] <- df$w[k] }
  tu <- which(rowSums(W) > 0); tp <- as.integer(sum(W[tu[1], ]))
  tm <- trop_matrix(Y, W, tu, 0, 0, Inf, tp)
  did <- trop(df, "y", "w", "id", "t",
              lambda = list(time = 0, unit = 0, nn = Inf), se = "none")$estimate
  expect_equal(tm, did, tolerance = 1e-3)
})

test_that("trop_matrix validates inputs", {
  df <- sim_panel(N = 10, T = 10, n_treated = 2, t0 = 8, seed = 1)
  N <- max(df$id); T <- max(df$t)
  Y <- matrix(0, N, T); W <- matrix(0, N, T)
  for (k in seq_len(nrow(df))) { Y[df$id[k], df$t[k]] <- df$y[k]
    W[df$id[k], df$t[k]] <- df$w[k] }
  tu <- which(rowSums(W) > 0)
  expect_error(trop_matrix(Y[, -1], W, tu, 0, 0, Inf, 3))
  expect_error(trop_matrix(Y, W * 0, tu, 0, 0, Inf, 3))
  expect_error(trop_matrix(Y, W, tu, -1, 0, Inf, 3))
})

test_that("trop(anchor = 'pooled') matches trop_matrix on a trailing block", {
  # Both now use the same pooled convention (unit distances to the average
  # treated trajectory; time distances to the centre of the treated block), so
  # for fixed penalties on a trailing-block design they must agree to solver
  # tolerance -- and hence match the official Python TROP_TWFE_average
  # convention through the same code path.
  df <- sim_panel(N = 16, T = 12, n_treated = 3, t0 = 9, att = 2, seed = 7)
  pm <- panel_matrices(df, "y", "w", "id", "t")
  for (lam in list(c(0.4, 0.3, Inf), c(0.6, 0.2, 2))) {
    tm <- trop_matrix(pm$Y, pm$W, pm$treated_units,
                      lambda_unit = lam[1], lambda_time = lam[2],
                      lambda_nn = lam[3],
                      treated_periods = pm$treated_periods)
    # trop_matrix() keeps the Python time-centre convention (T - tp/2); ask
    # trop() for the same convention for a like-for-like comparison.
    fit <- trop(df, "y", "w", "id", "t",
                lambda = list(unit = lam[1], time = lam[2], nn = lam[3]),
                anchor = "pooled", se = "none",
                control = trop_control(time_center = "python"))
    expect_equal(fit$estimate, tm, tolerance = 1e-4)
  }
})

test_that("pooled time anchor tracks a non-trailing treated block", {
  # treated block in the middle of the panel: the centre must be the centre of
  # the actual block, not of a hypothetical trailing block. The default is the
  # true midpoint (Stata >= 0.2.4); "python" sits half a period later.
  expect_equal(cfcompare:::.treated_block_center(5:8), 6.5)
  expect_equal(cfcompare:::.treated_block_center(9:12), 10.5)  # trailing, T = 12
  expect_equal(cfcompare:::.treated_block_center(7), 7)        # single period
  expect_equal(cfcompare:::.treated_block_center(5:8, "python"), 7)
  # "python" trailing block equivalence with the (T - tp/2) + 1 formula
  Tt <- 12; tp <- 4
  expect_equal(cfcompare:::.treated_block_center((Tt - tp + 1):Tt, "python"),
               (Tt - tp / 2) + 1)
})
