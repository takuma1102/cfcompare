# Pooled anchoring over treatment spell groups (the Stata command's group(time)
# estimands) under staggered and non-absorbing (on-off) assignment.
#
# Conventions checked here, all following the official Stata command (v0.2.6+):
# * one pooled estimand per (start, end) spell, pooled across units; the ATT is
#   the treated-cell-weighted average of the group effects;
# * the group's unit distance is to the group's average path over the periods
#   in which no member is treated, with each unit's own treated cells masked;
# * the time weights are centred on the midpoint of the spell, so a single
#   treated cell gets the same estimate under the pooled and per-cell anchors;
# * a block design is a single group, so block results are unchanged.

.make_panel2 <- function(N = 20, Tt = 12, seed = 1) {
  set.seed(seed)
  outer(rnorm(N), rep(1, Tt)) + outer(rep(1, N), rnorm(Tt)) +
    0.4 * outer(rnorm(N), rnorm(Tt)) + matrix(rnorm(N * Tt, sd = .2), N, Tt)
}

.long <- function(Y, W) {
  N <- nrow(Y); Tt <- ncol(Y)
  data.frame(id = rep(seq_len(N), Tt), t = rep(seq_len(Tt), each = N),
             y = as.vector(Y), w = as.vector(W))
}

test_that(".treatment_spells() reproduces the Stata README example", {
  # The 5 x 6 example of the Stata command's README (group(time)): spells
  # {3,4} for units 1-2 (pooled), {4,5,6} for unit 3, {2,3} and {5,6} for
  # unit 4 -- four groups of sizes 4, 3, 2, 2.
  W <- rbind(c(0, 0, 1, 1, 0, 0),
             c(0, 0, 1, 1, 0, 0),
             c(0, 0, 0, 1, 1, 1),
             c(0, 1, 1, 0, 1, 1),
             c(0, 0, 0, 0, 0, 0))
  sp <- cfcompare:::.treatment_spells(W)
  expect_length(sp, 4L)
  expect_equal(vapply(sp, `[[`, character(1), "key"),
               c("2_3", "3_4", "4_6", "5_6"))
  expect_equal(vapply(sp, `[[`, numeric(1), "size"), c(2, 4, 3, 2))
  expect_equal(sp[[2L]]$us, c(1L, 2L))
  expect_equal(sp[[3L]]$cols, 4:6)
  # a unit with two spells belongs to two groups
  expect_true(4L %in% sp[[1L]]$us && 4L %in% sp[[4L]]$us)
  # every treated cell is covered exactly once
  cover <- matrix(0, nrow(W), ncol(W))
  for (g in sp) cover[g$us, g$cols] <- cover[g$us, g$cols] + 1
  expect_equal(cover, W)
  # the pattern description carries the spells
  pat <- cfcompare:::.assignment_pattern(W)
  expect_equal(pat$n_spells, 4L)
  expect_length(cfcompare:::.treatment_spells(W * 0), 0L)
})

test_that("group unit distance masks the members' other spells and the
          comparison units' treated cells", {
  Y <- .make_panel2()
  W <- matrix(0, nrow(Y), ncol(Y))
  W[1, c(3, 4, 9, 10)] <- 1          # unit 1: two spells
  W[2, 5:6] <- 1                      # unit 2: a spell that does not overlap
  g <- cfcompare:::.treatment_spells(W)[[1L]]      # spell {3,4} of unit 1
  d0 <- cfcompare:::.unit_distance_to_avg(Y, W, g$us, g$cols)
  # the member's own treated outcomes (both spells) do not enter the
  # reference path ...
  Y1 <- Y; Y1[1, c(3, 4, 9, 10)] <- 1e6
  expect_equal(cfcompare:::.unit_distance_to_avg(Y1, W, g$us, g$cols), d0)
  # ... and a comparison unit's own treated cells never enter its distance,
  # even though they do not overlap the target spell
  Y2 <- Y; Y2[2, 5:6] <- 1e6
  expect_equal(cfcompare:::.unit_distance_to_avg(Y2, W, g$us, g$cols), d0)
  # whereas an untreated period of a comparison unit does
  Y3 <- Y; Y3[2, 7] <- 1e6
  expect_gt(cfcompare:::.unit_distance_to_avg(Y3, W, g$us, g$cols)[2], d0[2])
  # the member is compared to the group average like any other unit: with a
  # single member its own distance is zero
  expect_equal(d0[1], 0)
})

test_that("block designs are a single group and the pooled fit is unchanged", {
  Y <- .make_panel2(seed = 2)
  W <- matrix(0, nrow(Y), ncol(Y)); W[1:3, 9:12] <- 1
  pat <- cfcompare:::.assignment_pattern(W)
  expect_identical(pat$type, "block")
  expect_equal(pat$n_spells, 1L)
  lam <- list(time = 0.3, unit = 0.5, nn = 2, nn_floor = 0)
  ctrl <- trop_control(nn_floor = 0)
  est <- cfcompare:::.trop_att(Y, W, lam, ctrl, "pooled", pat)
  # direct single-anchor solve (the previous pooled code path)
  du <- cfcompare:::.unit_distance_to_avg(Y, W, 1:3)
  wmat <- cfcompare:::.trop_weight_matrix(
    du, cfcompare:::.treated_block_center(9:12), ncol(Y), lam)
  M <- cfcompare:::.trop_solve(Y, W, wmat, lam, ctrl)$M
  expect_equal(est$Mhat, M)
  expect_equal(est$att, mean((Y - M)[W == 1]))
  expect_equal(nrow(est$tau_groups), 1L)
  expect_equal(est$tau_groups$estimate, est$att)
  expect_equal(est$tau_groups$n_cells, 12L)
})

test_that("staggered / on-off designs: one solve per spell, cell-weighted ATT", {
  Y <- .make_panel2(N = 24, seed = 3)
  W <- matrix(0, nrow(Y), ncol(Y))
  W[1:2, 9:12] <- 1                   # cohort A
  W[3, 11:12] <- 1                    # cohort B
  W[4, c(3, 4, 9, 10)] <- 1           # on-off unit: spells {3,4}, {9,10}
  pat <- cfcompare:::.assignment_pattern(W)
  expect_equal(pat$n_spells, 4L)
  lam <- list(time = 0.3, unit = 0.5, nn = 2, nn_floor = 0)
  # tight tolerance: groups warm-start from the previous group's solution,
  # which changes only the iteration count of the convex program -- up to tol.
  ctrl <- trop_control(nn_floor = 0, tol = 1e-9)
  est <- cfcompare:::.trop_att(Y, W, lam, ctrl, "pooled", pat)
  tg <- est$tau_groups
  expect_equal(nrow(tg), 4L)
  expect_equal(sum(tg$n_cells), sum(W))
  # ATT = cell-weighted average of the group effects (Stata's e(group_weight))
  expect_equal(est$att, sum(tg$estimate * tg$n_cells) / sum(tg$n_cells))
  expect_equal(est$att, mean(est$tau_cells$tau))
  # each group effect equals the mean gap of a single-group solve anchored to
  # that group alone, with all treated cells excluded from the loss
  for (g in pat$spells) {
    du <- cfcompare:::.unit_distance_to_avg(Y, W, g$us, g$cols)
    wmat <- cfcompare:::.trop_weight_matrix(
      du, cfcompare:::.treated_block_center(g$cols), ncol(Y), lam)
    Mg <- cfcompare:::.trop_solve(Y, W, wmat, lam, ctrl)$M
    expect_equal(est$Mhat[g$us, g$cols], Mg[g$us, g$cols], tolerance = 1e-5)
    tau_g <- mean((Y - Mg)[g$us, g$cols])
    expect_equal(tg$estimate[tg$start == g$start & tg$end == g$end], tau_g,
                 tolerance = 1e-5)
  }
  # the full counterfactual matrix has no holes
  expect_false(anyNA(est$Mhat))
  # the user-facing fit exposes the group table and prints it
  fit <- trop(.long(Y, W), "y", "w", "id", "t", lambda = lam[1:3],
              anchor = "pooled", se = "none", control = ctrl)
  expect_equal(fit$tau_groups$estimate, tg$estimate)
  expect_output(print(fit), "4 treatment spell group")
  # the pooled counterfactual used by the event study / plots is the same
  # (re-solved at the default solver tolerance, hence the loose comparison)
  expect_equal(counterfactual_matrix(fit), est$Mhat, ignore_attr = TRUE,
               tolerance = 1e-4)
})

test_that("a single treated cell gets the same estimate under both anchors", {
  # Pooled time weights are centred on the spell midpoint, which for one period
  # is the paper's |s - t|; the group distance with one member equals the
  # per-cell distance. (Stata command, v0.2.4 release note.)
  Y <- .make_panel2(seed = 4)
  W <- matrix(0, nrow(Y), ncol(Y)); W[5, 8] <- 1
  pat <- cfcompare:::.assignment_pattern(W)
  for (lam in list(list(time = 0.4, unit = 0.6, nn = Inf),
                   list(time = 0.4, unit = 0.6, nn = 3))) {
    lam$nn_floor <- 0
    ctrl <- trop_control(nn_floor = 0, tol = 1e-9)
    a_pool <- cfcompare:::.trop_att(Y, W, lam, ctrl, "pooled",   pat)$att
    a_cell <- cfcompare:::.trop_att(Y, W, lam, ctrl, "per_cell", pat)$att
    expect_equal(a_pool, a_cell, tolerance = 1e-5)
    # ... whereas the Python-style centre (half a period late) differs
    ctrl_py <- trop_control(nn_floor = 0, time_center = "python")
    a_py <- cfcompare:::.trop_att(Y, W, lam, ctrl_py, "pooled", pat)$att
    expect_false(isTRUE(all.equal(a_py, a_cell, tolerance = 1e-6)))
  }
})

test_that("uniform weights solve once (DID / MC) and match the grouped solve", {
  Y <- .make_panel2(seed = 5)
  W <- matrix(0, nrow(Y), ncol(Y))
  W[1:2, 9:12] <- 1; W[3, 11:12] <- 1
  pat <- cfcompare:::.assignment_pattern(W)
  lam <- list(time = 0, unit = 0, nn = Inf, nn_floor = 0)
  ctrl <- trop_control(nn_floor = 0)
  est <- cfcompare:::.trop_att(Y, W, lam, ctrl, "pooled", pat)
  # the single uniform-weight solve is the plain TWFE fit on the control cells
  wmat <- matrix(1, nrow(Y), ncol(Y))
  M <- cfcompare:::.trop_solve(Y, W, wmat, lam, ctrl)$M
  expect_equal(est$Mhat, M)
})

test_that("placebo CV groups the placebo draws by spell", {
  Y <- .make_panel2(N = 18, seed = 6)
  W <- matrix(0, nrow(Y), ncol(Y))
  W[1, c(3, 4, 9, 10)] <- 1          # on-off unit: two spells
  W[2, 9:10] <- 1                     # shares the second spell
  pat <- cfcompare:::.assignment_pattern(W)
  ctrl <- trop_control(n_placebo = 3L, seed = 9L)
  pb <- cfcompare:::.trop_placebo_setup(Y, W, pat, ctrl)
  for (d in pb$draws) {
    expect_equal(sum(d$Wp), sum(W))
    # two spell groups: {3,4} (one unit) and {9,10} (two units)
    expect_length(d$groups, 2L)
    expect_equal(vapply(d$groups, `[[`, numeric(1), "size"), c(2, 4))
    for (g in d$groups) {
      expect_equal(g$t_center, cfcompare:::.treated_block_center(g$cols))
      expect_equal(g$du,
                   cfcompare:::.unit_distance_to_avg(pb$Yp, d$Wp, g$us, g$cols))
    }
  }
  lam <- list(time = 0.1, unit = 0.1, nn = Inf)
  q <- cfcompare:::.trop_cv_Q_placebo(pb, lam, ctrl)
  expect_true(is.finite(q$Q) && q$Q >= 0)
})
