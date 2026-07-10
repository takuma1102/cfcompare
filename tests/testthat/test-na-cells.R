# NA outcome cells (unbalanced panels) must never enter the solver loss, the
# LOOCV cell sample, or the placebo-CV scoring. Regression tests for the bug
# where a control cell with a missing outcome was included in the loss with the
# zero-filled target Yf = 0, dragging the fit (and the unit/time fixed effects)
# towards zero, and where a held-out / placebo cell landing on an NA hole made
# the whole cross-validation criterion NA.

.na_panel <- function(N = 12, Tt = 10, level = 50, seed = 2) {
  set.seed(seed)
  level + outer(rnorm(N), rep(1, Tt)) + outer(rep(1, N), rnorm(Tt)) +
    matrix(rnorm(N * Tt, sd = .2), N, Tt)
}

test_that("non-finite cells are excluded from the solver loss", {
  Y <- .na_panel()
  N <- nrow(Y); Tt <- ncol(Y)
  W <- matrix(0, N, Tt); W[1, 8:10] <- 1
  Yna <- Y; Yna[5, 3] <- NA
  w <- matrix(1, N, Tt)

  # an NA cell must be numerically equivalent to excluding that cell from the
  # loss explicitly via the mask ...
  fa <- cfcompare:::.mcnnm_fit(Yna, (W == 0) * 1, w, 5)
  mb <- (W == 0) * 1; mb[5, 3] <- 0
  fb <- cfcompare:::.mcnnm_fit(Yna, mb, w, 5)
  expect_equal(fa$M, fb$M, tolerance = 1e-8)

  # ... so the completion at the NA cell sits at the data level (~ 50), not
  # dragged towards the zero-filled target, and the unit's fixed effect is not
  # contaminated.
  expect_lt(abs(fa$M[5, 3] - Y[5, 3]), 5)
  fc <- cfcompare:::.mcnnm_fit(Y, (W == 0) * 1, w, 5)
  expect_lt(abs(fa$alpha[5] - fc$alpha[5]), 1)

  # complete panels are untouched (the exclusion is a no-op there)
  expect_equal(cfcompare:::.mcnnm_fit(Y, (W == 0) * 1, w, 5)$M, fc$M)
})

test_that("LOOCV never holds out an unobserved cell and stays finite", {
  Y <- .na_panel(seed = 3)
  N <- nrow(Y); Tt <- ncol(Y)
  W <- matrix(0, N, Tt); W[1, 8:10] <- 1; W[2, 9:10] <- 1
  Yna <- Y; Yna[5, 3] <- NA; Yna[7, 6] <- NA

  cv <- cfcompare:::.sample_control_cells(W, Inf, Y = Yna)
  expect_false(any(is.na(Yna[cv])))
  # backward compatible: without Y the previous behaviour is unchanged
  expect_equal(nrow(cfcompare:::.sample_control_cells(W, Inf)), sum(W == 0))

  ctrl <- cfcompare::trop_control(seed = 3)
  q <- cfcompare:::.trop_cv_Q(Yna, W, list(time = 0, unit = 0, nn = Inf),
                              ctrl, cv)
  expect_true(is.finite(q$Q))
})

test_that("placebo CV scores only observed cells and stays finite", {
  skip_on_cran()
  Y <- .na_panel(seed = 4)
  N <- nrow(Y); Tt <- ncol(Y)
  W <- matrix(0, N, Tt); W[1, 8:10] <- 1; W[2, 9:10] <- 1
  # NA holes in post periods of control units: placebo cells can land on them
  Yna <- Y; Yna[6, 9] <- NA; Yna[8, 10] <- NA
  pat <- cfcompare:::.assignment_pattern(W)
  lam <- list(time = .1, unit = .1, nn = Inf)
  for (sd in 1:10) {
    ctrl <- cfcompare::trop_control(n_placebo = 5L, seed = sd)
    pb <- cfcompare:::.trop_placebo_setup(Yna, W, pat, ctrl)
    q <- cfcompare:::.trop_cv_Q_placebo(pb, lam, ctrl)
    expect_true(is.finite(q$Q) && q$Q >= 0)
  }
  # if NO placebo cell is ever observed, the criterion must fail loudly rather
  # than return NA (every post-period control cell missing)
  Ybad <- Y; Ybad[3:N, 8:10] <- NA
  ctrl <- cfcompare::trop_control(n_placebo = 3L, seed = 1)
  pb <- cfcompare:::.trop_placebo_setup(Ybad, W, pat, ctrl)
  expect_error(cfcompare:::.trop_cv_Q_placebo(pb, lam, ctrl), "scoreable")
})

test_that("trop() runs end-to-end on an unbalanced panel", {
  skip_on_cran()
  Y <- .na_panel(N = 16, Tt = 12, seed = 5)
  W <- matrix(0, nrow(Y), ncol(Y)); W[1, 9:12] <- 1; W[2, 11:12] <- 1
  df <- data.frame(id = rep(seq_len(nrow(Y)), each = ncol(Y)),
                   t = rep(seq_len(ncol(Y)), nrow(Y)),
                   y = as.numeric(t(Y)), w = as.numeric(t(W)))
  df <- df[-c(15, 40, 77), ]            # drop rows -> NA holes in the matrix
  m <- cfcompare:::.panel_to_matrices(df, "y", "w", "id", "t")
  expect_true(anyNA(m$Y))

  # fixed penalties
  f0 <- trop(df, "y", "w", "id", "t", se = "none",
             lambda = list(time = .1, unit = .1, nn = 2))
  expect_true(is.finite(f0$estimate))

  # LOOCV penalty selection (previously Q = NA broke the search)
  f1 <- suppressWarnings(trop(df, "y", "w", "id", "t", se = "none",
             control = cfcompare::trop_control(n_cv_cells = 6L,
                                               cv_cycles = 1L, seed = 2),
             grids = list(time = c(0, .2), unit = c(0, .5), nn = c(Inf, 1))))
  expect_true(is.finite(f1$estimate))

  # placebo-CV penalty selection with the staggered pattern bank
  f2 <- trop(df, "y", "w", "id", "t", se = "none",
             control = cfcompare::trop_control(cv_method = "placebo",
                                               n_placebo = 3L,
                                               cv_cycles = 1L, seed = 7),
             grids = list(time = c(0, .2), unit = c(0, .5), nn = c(Inf, 1)))
  expect_true(is.finite(f2$estimate))
})
