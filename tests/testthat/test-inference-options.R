# Own-group distance overrides (the Stata command's pooled_treat_distance()),
# the bootstrap resampling scheme (the Python package's resample = "controls")
# and the delete-one-unit jackknife over every unit (the Stata command's
# vce(jackknife)).

.make_panel3 <- function(N = 20, Tt = 12, seed = 1) {
  set.seed(seed)
  outer(rnorm(N), rep(1, Tt)) + outer(rep(1, N), rnorm(Tt)) +
    0.4 * outer(rnorm(N), rnorm(Tt)) + matrix(rnorm(N * Tt, sd = .2), N, Tt)
}

.long3 <- function(Y, W) {
  N <- nrow(Y); Tt <- ncol(Y)
  data.frame(id = rep(seq_len(N), Tt), t = rep(seq_len(Tt), each = N),
             y = as.vector(Y), w = as.vector(W))
}

test_that("pooled_treat_distance overrides the members' unit and in-spell time distances", {
  Y <- .make_panel3()
  W <- matrix(0, nrow(Y), ncol(Y)); W[1:3, 9:12] <- 1
  g <- cfcompare:::.treatment_spells(W)[[1L]]
  d0 <- cfcompare:::.unit_distance_to_avg(Y, W, g$us, g$cols)
  d1 <- cfcompare:::.unit_distance_to_avg(Y, W, g$us, g$cols, self_dist = 0)
  # members are generally not at zero distance from their own average ...
  expect_true(any(d0[g$us] > 0))
  # ... until overridden; non-members are untouched
  expect_equal(d1[g$us], c(0, 0, 0))
  expect_equal(d1[-g$us], d0[-g$us])
  lam <- list(time = 0.4, unit = 0.6, nn = Inf)
  wm0 <- cfcompare:::.trop_weight_matrix(d0, 10.5, ncol(Y), lam)
  wm1 <- cfcompare:::.trop_weight_matrix(d0, 10.5, ncol(Y), lam,
                                         self_cols = g$cols, self_time = 0)
  # in-spell periods get the full time weight, other periods are unchanged
  expect_equal(wm1[, g$cols] / exp(-lam$unit * d0), matrix(1, nrow(Y), 4),
               ignore_attr = TRUE)
  expect_equal(wm1[, 1:8], wm0[, 1:8])
  # a unit value alone leaves the time distances alone (and vice versa)
  wm2 <- cfcompare:::.trop_weight_matrix(d0, 10.5, ncol(Y), lam,
                                         self_cols = g$cols, self_time = NULL)
  expect_equal(wm2, wm0)
})

test_that("pooled_treat_distance is validated, stamped into the fit and reproduced", {
  expect_error(trop_control(pooled_treat_distance = list(foo = 1)),
               "pooled_treat_distance")
  expect_error(trop_control(pooled_treat_distance = list(unit = -1)),
               "non-negative")
  expect_error(trop_control(pooled_treat_distance = 0), "pooled_treat_distance")
  expect_null(trop_control()$pooled_treat_distance$unit)
  ctrl <- trop_control(pooled_treat_distance = list(unit = 0))
  expect_equal(ctrl$pooled_treat_distance, list(unit = 0, time = NULL))

  Y <- .make_panel3(seed = 2)
  W <- matrix(0, nrow(Y), ncol(Y)); W[1:3, 9:12] <- 1
  lam <- list(time = 0.4, unit = 0.6, nn = Inf)
  ctrl0 <- trop_control(nn_floor = 0)
  ctrl1 <- trop_control(nn_floor = 0, pooled_treat_distance = list(unit = 0, time = 0))
  f0 <- trop(.long3(Y, W), "y", "w", "id", "t", lambda = lam, anchor = "pooled",
             se = "none", control = ctrl0)
  f1 <- trop(.long3(Y, W), "y", "w", "id", "t", lambda = lam, anchor = "pooled",
             se = "none", control = ctrl1)
  expect_false(isTRUE(all.equal(f0$estimate, f1$estimate)))
  expect_equal(f1$lambda$pooled_treat_distance, list(unit = 0, time = 0))
  expect_output(print(f1), "own-group distances overridden")
  # the fit's setting travels with $lambda: a manual solve with the overridden
  # weights reproduces the estimate, and so does the pooled counterfactual
  # behind counterfactual_matrix() even under a fresh trop_control()
  g <- f1$pattern$spells[[1L]]
  du <- cfcompare:::.unit_distance_to_avg(Y, W, g$us, g$cols, 0)
  wm <- cfcompare:::.trop_weight_matrix(du, cfcompare:::.treated_block_center(g$cols),
                                        ncol(Y), f1$lambda, self_cols = g$cols,
                                        self_time = 0)
  M <- cfcompare:::.trop_solve(Y, W, wm, f1$lambda, ctrl1)$M
  expect_equal(f1$estimate, mean((Y - M)[W == 1]), tolerance = 1e-6)
  expect_equal(counterfactual_matrix(f1), M, ignore_attr = TRUE, tolerance = 1e-5)
  # the per-cell anchor ignores the option
  p0 <- trop(.long3(Y, W), "y", "w", "id", "t", lambda = lam, anchor = "per_cell",
             se = "none", control = ctrl0)
  p1 <- trop(.long3(Y, W), "y", "w", "id", "t", lambda = lam, anchor = "per_cell",
             se = "none", control = ctrl1)
  expect_equal(p0$estimate, p1$estimate)
  # and the placebo-CV setup carries it into the placebo draws
  pat <- cfcompare:::.assignment_pattern(W)
  pb <- cfcompare:::.trop_placebo_setup(Y, W, pat, trop_control(n_placebo = 2L, seed = 3L,
                                        pooled_treat_distance = list(unit = 0, time = 0)))
  for (d in pb$draws) for (gr in d$groups) {
    expect_equal(gr$du[gr$us], rep(0, length(gr$us)))
    expect_equal(gr$self_time, 0)
  }
})

test_that("bootstrap resampling: stratified (default) vs controls-only", {
  tu <- c(2L, 5L); co <- c(1L, 3L, 4L, 6L, 7L)
  set.seed(1)
  s <- cfcompare:::.boot_index_sets(tu, co, 50, "stratified")
  c1 <- cfcompare:::.boot_index_sets(tu, co, 50, "controls")
  expect_length(s, 50L); expect_length(c1, 50L)
  # both keep the unit counts fixed and draw controls with replacement
  expect_true(all(vapply(s, length, integer(1)) == 7L))
  expect_true(all(vapply(c1, length, integer(1)) == 7L))
  # controls-only keeps the treated rows exactly; stratified resamples them
  expect_true(all(vapply(c1, function(i) identical(i[1:2], tu), logical(1))))
  expect_true(any(vapply(s, function(i) !identical(i[1:2], tu), logical(1))))
  expect_true(all(unlist(lapply(c1, function(i) i[3:7] %in% co))))
  # a single treated unit is drawn from {tu}, never from 1:tu
  s1 <- cfcompare:::.boot_index_sets(5L, co, 20, "stratified")
  expect_true(all(vapply(s1, function(i) i[1] == 5L, logical(1))))

  Y <- .make_panel3(N = 16, seed = 4)
  W <- matrix(0, nrow(Y), ncol(Y)); W[1, 9:12] <- 1    # single treated unit
  lam <- list(time = 0.3, unit = 0.5, nn = Inf)
  f_c <- trop(.long3(Y, W), "y", "w", "id", "t", lambda = lam, anchor = "pooled",
              se = "bootstrap",
              control = trop_control(n_boot = 25L, seed = 2L, boot_resample = "controls"))
  f_s <- trop(.long3(Y, W), "y", "w", "id", "t", lambda = lam, anchor = "pooled",
              se = "bootstrap",
              control = trop_control(n_boot = 25L, seed = 2L, boot_resample = "stratified"))
  expect_true(is.finite(f_c$std.error) && is.finite(f_s$std.error))
  expect_identical(f_c$boot_resample, "controls")
  expect_identical(f_s$boot_resample, "stratified")
  expect_output(print(f_c), "controls resampled")
  expect_equal(f_c$estimate, f_s$estimate)    # point estimate is unaffected
  # with one treated unit the treated row cannot vary under either scheme, so
  # the two SEs differ only through the random control draws
  expect_equal(f_c$std.error, f_s$std.error, tolerance = 0.1)
  # the event study accepts the option as well
  es <- trop_event_study(f_c, se = "bootstrap", pre_periods = FALSE,
                         control = trop_control(n_boot = 10L, seed = 2L,
                                                boot_resample = "controls"))
  expect_true(all(is.finite(es$estimates$std.error)))
})

test_that("jackknife deletes every unit and matches a manual computation", {
  Y <- .make_panel3(N = 14, seed = 5)
  W <- matrix(0, nrow(Y), ncol(Y)); W[1:2, 9:12] <- 1
  lam <- list(time = 0.3, unit = 0.5, nn = Inf)
  ctrl <- trop_control(nn_floor = 0)
  fit <- trop(.long3(Y, W), "y", "w", "id", "t", lambda = lam, anchor = "pooled",
              se = "jackknife", control = ctrl)
  # manual delete-one-unit over all N units (Stata's vce(jackknife))
  N <- nrow(Y)
  a <- vapply(seq_len(N), function(i) {
    keep <- setdiff(seq_len(N), i)
    Yk <- Y[keep, ]; Wk <- W[keep, ]
    cfcompare:::.trop_att(Yk, Wk, fit$lambda, ctrl, "pooled",
                          cfcompare:::.assignment_pattern(Wk))$att
  }, numeric(1))
  se_manual <- sqrt((N - 1) / N * sum((a - mean(a))^2))
  expect_equal(fit$std.error, se_manual, tolerance = 1e-6)
  # replicates are over every unit, so the SE differs from the previous
  # treated-only jackknife
  a_tr <- a[1:2]
  expect_false(isTRUE(all.equal(fit$std.error,
                                sqrt(1 / 2 * sum((a_tr - mean(a_tr))^2)))))
  # a single treated unit: its own deletion is skipped, the remaining
  # (control-deletion) replicates still give a finite SE
  W1 <- W; W1[2, ] <- 0
  fit1 <- trop(.long3(Y, W1), "y", "w", "id", "t", lambda = lam, anchor = "pooled",
               se = "jackknife", control = ctrl)
  expect_true(is.finite(fit1$std.error))
  # the event study jackknife follows the same scheme
  es <- trop_event_study(fit, se = "jackknife", pre_periods = FALSE, control = ctrl)
  expect_true(all(is.finite(es$estimates$std.error)))
})
