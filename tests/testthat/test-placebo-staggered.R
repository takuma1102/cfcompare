# Pattern-aware placebo cross-validation (staggered / block adoption).
#
# The placebo criterion must stamp the *actual* treated-unit adoption patterns
# onto resampled control units (as the official Stata command does), so that
# under staggered adoption each placebo draw reproduces the real cohort
# structure, and under block adoption it reduces to the previous
# common-post-block behaviour.

.make_panel <- function(N = 14, Tt = 12, seed = 1) {
  set.seed(seed)
  u <- rnorm(N); v <- rnorm(Tt)
  Y <- outer(u, rep(1, Tt)) + outer(rep(1, N), v) +
    0.5 * outer(rnorm(N), rnorm(Tt)) + matrix(rnorm(N * Tt, sd = .2), N, Tt)
  Y
}

test_that("placebo setup stamps staggered adoption patterns and cohorts", {
  skip_on_cran()
  Y <- .make_panel()
  W <- matrix(0, nrow(Y), ncol(Y))
  W[1, 9:12] <- 1                      # cohort A: adopts at t = 9
  W[2, 11:12] <- 1                     # cohort B: adopts at t = 11
  W[3, 11:12] <- 1                     # cohort B
  pat <- cfcompare:::.assignment_pattern(W)
  expect_identical(pat$type, "general")

  ctrl <- cfcompare::trop_control(n_placebo = 5L, seed = 42L)
  pb <- cfcompare:::.trop_placebo_setup(Y, W, pat, ctrl)

  # control panel drops the real treated rows
  expect_equal(nrow(pb$Yp), nrow(Y) - 3L)

  for (d in pb$draws) {
    # every draw has exactly the real number of treated cells,
    # with the real column profile (1 unit from t=9, 2 units from t=11)
    expect_equal(sum(d$Wp), sum(W))
    expect_equal(colSums(d$Wp), colSums(W[1:3, , drop = FALSE]),
                 ignore_attr = TRUE)
    # two adoption cohorts, sizes 1x4 and 2x2 cells
    expect_length(d$groups, 2L)
    sizes <- sort(vapply(d$groups, `[[`, numeric(1), "size"))
    expect_equal(sizes, c(4, 4))
    cols <- lapply(d$groups, `[[`, "cols")
    expect_setequal(vapply(cols, min, numeric(1)), c(9, 11))
    # per-cohort anchors and distances are cached
    for (g in d$groups) {
      expect_equal(g$t_center,
                   cfcompare:::.treated_block_center(g$cols))
      expect_length(g$du, nrow(pb$Yp))
    }
  }
})

test_that("block designs reduce to a single cohort (previous behaviour)", {
  skip_on_cran()
  Y <- .make_panel(seed = 2)
  W <- matrix(0, nrow(Y), ncol(Y))
  W[1:2, 10:12] <- 1
  pat <- cfcompare:::.assignment_pattern(W)
  expect_identical(pat$type, "block")

  ctrl <- cfcompare::trop_control(n_placebo = 4L, seed = 7L)
  pb <- cfcompare:::.trop_placebo_setup(Y, W, pat, ctrl)
  for (d in pb$draws) {
    expect_length(d$groups, 1L)
    expect_equal(d$groups[[1L]]$cols, 10:12)
  }

  # criterion runs and returns a finite scalar, with warm states per cohort
  lam <- list(time = 0.1, unit = 0.1, nn = Inf)
  q <- cfcompare:::.trop_cv_Q_placebo(pb, lam, ctrl, keep_state = TRUE)
  expect_true(is.finite(q$Q) && q$Q >= 0)
  expect_length(q$warm, length(pb$draws))
  expect_length(q$warm[[1L]], 1L)

  # warm-started re-evaluation at a weaker penalty agrees with a cold start.
  # The program is convex, so warm starts change only the iteration count, not
  # the solution -- but only up to the solver's stopping tolerance: at the
  # default tol the warm and cold iterates stop ~1e-4 (relative) apart, which a
  # 1e-6 comparison would flag. Solve the two evaluations to a tight tolerance
  # so the criterion gap (which scales with tol) sits well below the 1e-6
  # comparison level.
  ctrl_tight <- cfcompare::trop_control(n_placebo = 4L, seed = 7L, tol = 1e-10)
  lam2 <- list(time = 0.1, unit = 0.1, nn = 5)
  q_warm <- cfcompare:::.trop_cv_Q_placebo(pb, lam2, ctrl_tight, warm = q$warm)
  q_cold <- cfcompare:::.trop_cv_Q_placebo(pb, lam2, ctrl_tight)
  expect_equal(q_warm$Q, q_cold$Q, tolerance = 1e-6)
})

test_that("placebo criterion is finite on staggered designs and trop() runs", {
  skip_on_cran()
  Y <- .make_panel(N = 16, Tt = 12, seed = 3)
  W <- matrix(0, nrow(Y), ncol(Y))
  W[1, 9:12] <- 1
  W[2, 11:12] <- 1
  pat <- cfcompare:::.assignment_pattern(W)
  ctrl <- cfcompare::trop_control(n_placebo = 3L, seed = 11L,
                                  n_cv_cells = 6L, cv_cycles = 1L,
                                  cv_method = "placebo")
  pb <- cfcompare:::.trop_placebo_setup(Y, W, pat, ctrl)
  lam <- list(time = 0.1, unit = 0.1, nn = Inf)
  q <- cfcompare:::.trop_cv_Q_placebo(pb, lam, ctrl)
  expect_true(is.finite(q$Q) && q$Q >= 0)
})

test_that("loocv on a staggered design warns and suggests placebo CV", {
  skip_on_cran()
  Y <- .make_panel(N = 16, Tt = 12, seed = 4)
  W <- matrix(0, nrow(Y), ncol(Y))
  W[1, 9:12] <- 1
  W[2, 11:12] <- 1
  cv_cells <- cfcompare:::.sample_control_cells(W, 4L, seed = 5L)
  grids <- list(time = c(0, .5), unit = c(0, .5), nn = c(Inf, 20))
  ctrl <- cfcompare::trop_control(n_cv_cells = 4L, cv_cycles = 1L)
  expect_warning(
    cfcompare:::.trop_select_lambda(Y, W, grids, ctrl, cv_cells),
    regexp = "staggered"
  )
})
