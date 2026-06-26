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
