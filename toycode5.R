# =============================================================================
# cfcompare end-to-end toy smoke test
#
# This script now exercises the package's OWN helpers instead of re-implementing
# them locally. Functions that used to be hand-rolled here are called directly:
#   - sim_panel(n_cov=, phi=)   : covariate panels with known coefficients (2b)
#   - surface_matrix(), selected_lambda()                                  (5)
#   - bind_att()                : stack ATT results into one table         (7b)
#   - trop_ablation()           : Table-5-style penalty ablation           (7d)
#   - compare_se_modes()        : bootstrap / jackknife / placebo on one fit(8)
#   - panel_matrices()          : long -> Y/W (+ treated rows/periods)     (9a)
#   - trop_matrix()             : treated_units/periods inferred from W     (9a)
#   - true_att()                : ground-truth ATT of a simulated panel    (9b)
#   - rmse_curve(vary=)         : now also "n_treated" / "n_post"          (12)
#
# Demonstrates, in order:
#   (1) panel_compare()          : static multi-method ATT comparison
#   (2) trop() + covariates      : single fit, native covariate support
#   (3) trop_event_study()       : per-period/event-study TROP effects
#   (4) panel_rmse()             : RMSE comparison across estimators
#   (5) trop_sensitivity()       : lambda-grid heatmap + CV/ATT surfaces
#   (6) rmse_curve()/curves()    : RMSE vs design dimensions
#   (7) tidiers / anchors / ablation
#   (8) SE modes                 : bootstrap / jackknife / placebo on one fit
#   (9) sim_semisynthetic() + trop_matrix()
#  (10) optional engines + non-block designs
#  (11) custom grids / full SVD / optional workers
#  (12) design-dimension RMSE curves over N_treated / N_post
#  (13) session info
#
# Notes:
# - Smoke-test/demo script, not a paper-quality simulation. For final numbers,
#   raise n_runs, n_boot, n_cv_cells, cv_cycles, max_iter.
# - Optional engines (synthdid, gsynth, augsynth, did) and plotting/parallel
#   packages (ggplot2, patchwork, future) are assumed installed; panel_compare()
#   skips any estimator whose package is missing on its own, with a message.
# =============================================================================


## ---- 0. Setup --------------------------------------------------------------

# Set TRUE on the first run, or after you update the GitHub package.

remove.packages("cfcompare")

INSTALL_OR_UPDATE_CFCOMPARE <- FALSE

if (!requireNamespace("cfcompare", quietly = TRUE) || INSTALL_OR_UPDATE_CFCOMPARE) {
  if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
  pak::pak("takuma1102/cfcompare")
}

# When developing from the source tree (after applying a patch), load in place:
devtools::load_all()

library(cfcompare)

cat("cfcompare version:", as.character(utils::packageVersion("cfcompare")), "\n")

# Optional packages used below (install once):
# install.packages(c("ggplot2", "patchwork", "gsynth", "augsynth", "did",
#                    "future", "future.apply"))
# remotes::install_github("synth-inference/synthdid")

# Control objects.
# ctrl: moderate toy settings.
ctrl <- trop_control(
  n_cv_cells = 60L,
  cv_cycles  = 2L,
  max_iter   = 200L,
  conf_level = 0.95,
  n_boot     = 50L,
  boot_ci    = "percentile",
  seed       = 1L
)

# ctrl_fast: faster smoke-test settings.
ctrl_fast <- trop_control(
  n_cv_cells = 15L,
  cv_cycles  = 1L,
  max_iter   = 90L,
  conf_level = 0.95,
  n_boot     = 18L,
  boot_ci    = "percentile",
  seed       = 1L
)

# Event-study controls.
ctrl_event_fit <- trop_control(
  n_cv_cells = 30L,
  cv_cycles  = 1L,
  max_iter   = 150L,
  conf_level = 0.95,
  seed       = 7L
)

ctrl_event_boot <- trop_control(
  n_cv_cells = 30L,
  cv_cycles  = 1L,
  max_iter   = 150L,
  conf_level = 0.95,
  n_boot     = 40L,
  boot_ci    = "percentile",
  seed       = 7L
)


## ---- 1. Example data -------------------------------------------------------

# Synthetic block-treatment panel with known ATT. Columns: id, t, y, w, y0.
df <- sim_panel(
  N = 40, T = 10, n_treated = 4, t0 = 4,
  rank = 3, att = 3, noise = 1, seed = 1
)

head(df)

outcome   <- "y"
treatment <- "w"
unit      <- "id"
time      <- "t"

# Full method list; panel_compare() skips any estimator whose package is
# missing (e.g. SDID/SC without 'synthdid'), so no need to pre-filter.
methods <- c("DID", "SDID", "SC", "MC", "DIFP", "TROP")


## ---- 2. panel_compare(): multi-method ATT comparison -----------------------

cmp <- panel_compare(
  df, outcome, treatment, unit, time,
  methods = methods,
  anchor  = "auto",
  se      = "auto",
  control = ctrl_fast,
  verbose = TRUE
)

print(cmp)

att_tbl <- as_att(cmp)
print(att_tbl)

p_compare <- autoplot(cmp)
print(p_compare)
# ggplot2::ggsave("compare_forest.png", p_compare, width = 8, height = 5, dpi = 150)

# Single TROP fit. trop() defaults to bootstrap, so use se = "none" for speed.
fit_nose <- trop(
  df, outcome, treatment, unit, time,
  anchor  = "auto",
  se      = "none",
  control = ctrl_fast
)

# se = c("bootstrap", "auto", "jackknife", "placebo", "none")
fit <- trop(
  df, outcome, treatment, unit, time,
  anchor  = "auto",
  se      = "bootstrap",
  control = ctrl_fast
)

print(fit_nose)
print(fit)

p_trop <- autoplot(fit)
print(p_trop)
# ggplot2::ggsave("trop_trajectory.png", p_trop, width = 8, height = 5, dpi = 150)


## ---- 2b. Covariates and the common counterfactual matrix ------------------

# Native TROP covariate support (Section 6.2): the low-rank term is split
# additively as L = X.phi + R -- the covariate-linear part is unpenalised and the
# nuclear norm acts on the residual R.
#
# sim_panel(n_cov=, phi=) builds the covariates directly: x1, x2 are part of
# Y(0) (both y and the counterfactual y0 carry the signal), with the true
# coefficients stored on attr(df_cov, "phi"). trop(covariates=) recovers them.
df_cov <- sim_panel(
  N = 40, T = 10, n_treated = 4, t0 = 4,
  rank = 3, att = 3, noise = 1,
  n_cov = 2, phi = c(1.2, -0.7), seed = 1
)

cat("true covariate coefficients (attr 'phi'):\n")
print(attr(df_cov, "phi"))

fit_cov <- trop(
  df_cov, outcome, treatment, unit, time,
  covariates = c("x1", "x2"),
  anchor  = "pooled",
  se      = "none",
  control = ctrl_fast
)

print(fit_cov)
# Estimated covariate coefficients (should be near c(1.2, -0.7)).
print(fit_cov$phi)

# Covariates flow through the event study and SE resampling automatically.
es_cov <- trop_event_study(fit_cov, se = "none", control = ctrl_fast)
print(es_cov)

# counterfactual_matrix(): a method's estimated untreated outcome Y(0) as an
# N x T matrix. predict_counterfactual() is an alias.
M_trop <- counterfactual_matrix(fit)
M_cov  <- predict_counterfactual(fit_cov)
cat("counterfactual_matrix() dims:",
    paste(dim(M_trop), collapse = " x "),
    "| any NA (native, full):", anyNA(M_trop), "\n")

# On a panel_compare() result it returns one matrix per method. Native engines
# fill every cell; treated-only methods fill treated rows, leaving controls NA.
cf_list <- counterfactual_matrix(cmp)
cat("counterfactual_matrix(cmp) methods:",
    paste(names(cf_list), collapse = ", "), "\n")


## ---- 3. trop_event_study(): event-study API -------------------------------

# Event-study toy data with true post-treatment ATT = 3.
df_es <- sim_panel(
  N = 36, T = 18, n_treated = 7, t0 = 13,
  rank = 3, att = 3, noise = 1, seed = 7
)

true_att_es <- 3

fit_es <- trop(
  df_es, outcome, treatment, unit, time,
  anchor  = "pooled",
  se      = "none",
  control = ctrl_event_fit,
  verbose = TRUE
)

print(fit_es)

# Main event-study run: bootstrap SE, pre-period placebo points included.
es_boot <- trop_event_study(
  fit_es,
  se = "bootstrap",
  pre_periods = TRUE,
  control = ctrl_event_boot
)

print(es_boot)
print(es_boot$estimates)

p_es_boot <- autoplot(es_boot)
print(p_es_boot)
# ggplot2::ggsave("trop_event_study_bootstrap.png", p_es_boot, width = 8, height = 5, dpi = 150)

# Same point estimates, different SE methods.
es_jackknife <- trop_event_study(
  fit_es, se = "jackknife", pre_periods = TRUE, control = ctrl_event_fit
)

es_placebo <- trop_event_study(
  fit_es, se = "placebo", pre_periods = TRUE, control = ctrl_event_fit
)

es_none <- trop_event_study(
  fit_es, se = "none", pre_periods = TRUE, control = ctrl_event_fit
)

print(es_jackknife)
print(es_placebo)
print(es_none)

# Check pre_periods = FALSE.
es_post_only <- trop_event_study(
  fit_es, se = "none", pre_periods = FALSE, control = ctrl_event_fit
)

print(es_post_only)

p_es_post_only <- autoplot(es_post_only)
print(p_es_post_only)

# Compact diagnostic summary.
summarise_event_study <- function(es, true_att) {
  e <- es$estimates
  e$truth <- ifelse(e$event_time < 0, 0, true_att)
  e$ci_covers_truth <- with(
    e,
    is.finite(conf.low) & is.finite(conf.high) &
      conf.low <= truth & conf.high >= truth
  )

  data.frame(
    se_method = es$se.method,
    n_event_times = nrow(e),
    n_pre = sum(e$period == "pre"),
    n_post = sum(e$period == "post"),
    mean_pre_estimate = if (any(e$period == "pre")) {
      mean(e$estimate[e$period == "pre"], na.rm = TRUE)
    } else NA_real_,
    mean_post_estimate = if (any(e$period == "post")) {
      mean(e$estimate[e$period == "post"], na.rm = TRUE)
    } else NA_real_,
    share_ci_covering_truth = if (any(is.finite(e$conf.low) & is.finite(e$conf.high))) {
      mean(e$ci_covers_truth, na.rm = TRUE)
    } else NA_real_,
    stringsAsFactors = FALSE
  )
}

event_summary <- do.call(
  rbind,
  lapply(
    list(bootstrap = es_boot, jackknife = es_jackknife,
         placebo = es_placebo, none = es_none),
    summarise_event_study,
    true_att = true_att_es
  )
)

print(event_summary)

# Explicitly test the trop() default se = "bootstrap".
df_small <- sim_panel(
  N = 18, T = 10, n_treated = 4, t0 = 7,
  rank = 2, att = 2, noise = 1, seed = 77
)

fit_default_boot <- trop(
  df_small, outcome, treatment, unit, time,
  anchor  = "pooled",
  control = trop_control(
    n_cv_cells = 10L, cv_cycles = 1L, max_iter = 80L, n_boot = 12L, seed = 77L
  )
)

print(fit_default_boot)
stopifnot(identical(fit_default_boot$se.method, "bootstrap"))


## ---- 4. panel_rmse(): RMSE comparison across estimators --------------------

# Placebo RMSE: randomly pseudo-treat controls, true effect = 0.
rmse_pl <- panel_rmse(
  df, outcome, treatment, unit, time,
  methods  = methods,
  metric   = "placebo",
  n_pseudo = 6,
  n_runs   = 20,
  horizon  = 4,
  control  = ctrl_fast,
  seed     = 1,
  verbose  = TRUE
)

print(rmse_pl)

p_rmse_pl <- autoplot(rmse_pl)
print(p_rmse_pl)
# ggplot2::ggsave("rmse_placebo.png", p_rmse_pl, width = 8, height = 5, dpi = 150)

# Prediction RMSE: hold out a block of control cells and predict Y(0) there.
pred_methods <- c("DID", "MC", "TROP", "DIFP", "SDID", "SC")

rmse_pr <- panel_rmse(
  df, outcome, treatment, unit, time,
  methods = pred_methods,
  metric  = "prediction",
  horizon = 3,
  n_runs  = 10,
  control = ctrl_fast,
  seed    = 1
)

print(rmse_pr)

p_rmse_pr <- autoplot(rmse_pr)
print(p_rmse_pr)
# ggplot2::ggsave("rmse_prediction.png", p_rmse_pr, width = 8, height = 5, dpi = 150)


## ---- 5. trop_sensitivity(): lambda-grid heatmap + surfaces -----------------

heatmap_axes <- c("nn", "unit")
# Alternatives: c("nn", "time"), c("nn", "unit"), c("time", "unit")

stopifnot(
  length(heatmap_axes) == 2L,
  all(heatmap_axes %in% c("time", "unit", "nn")),
  length(unique(heatmap_axes)) == 2L
)

grid <- trop_sensitivity(
  df, outcome, treatment, unit, time,
  axes = heatmap_axes,
  lambda_time = if ("time" %in% heatmap_axes) c(0, 0.1, 0.25, 0.5, 1, 1.5) else NULL,
  lambda_unit = if ("unit" %in% heatmap_axes) c(0, 0.1, 0.25, 0.5, 1, 1.5) else NULL,
  lambda_nn   = NULL,
  anchor  = "pooled",
  control = ctrl_fast,
  seed    = 1
)

print(grid)
print(head(as.data.frame(grid)))

# CV-selected penalties. selected_lambda() returns the list(time=, unit=, nn=)
# that trop(lambda=) accepts, so the data-driven choice feeds straight back in.
print(attr(grid, "selected"))
lam_sel <- selected_lambda(grid)
print(lam_sel)
fit_at_selected <- trop(
  df, outcome, treatment, unit, time,
  lambda = lam_sel, anchor = "pooled", se = "none", control = ctrl_fast
)
print(fit_at_selected)

# Metadata.
print(attr(grid, "axes"))
print(attr(grid, "fixed"))

# surface_matrix(): wide layout of one quantity over the two swept penalties.
cv_surface  <- surface_matrix(grid, "cv_loss")
att_surface <- surface_matrix(grid, "att")
print(cv_surface)
print(att_surface)

p_heat <- autoplot(grid)
print(p_heat)
# ggplot2::ggsave("trop_sensitivity.png", p_heat, width = 8, height = 5, dpi = 150)

surfaces <- plot_trop_surfaces(grid, which = "both", ask = FALSE)
print(surfaces$cv_loss)
print(surfaces$att)


## ---- 5b. trop_sensitivity(): fixed-unit-lambda version ----------------------

fixed_lambda_unit <- 0.25

grid_fix_unit <- trop_sensitivity(
  df, outcome, treatment, unit, time,
  axes = c("nn", "time"),
  lambda_nn   = NULL,
  lambda_time = c(0, 0.1, 0.25, 0.5, 1, 1.5),
  lambda_unit = fixed_lambda_unit,
  anchor  = "pooled",
  control = ctrl_fast,
  seed    = 1
)

print(grid_fix_unit)
print(selected_lambda(grid_fix_unit))
print(attr(grid_fix_unit, "axes"))
print(attr(grid_fix_unit, "fixed"))

print(autoplot(grid_fix_unit))

surfaces_fix_unit <- plot_trop_surfaces(grid_fix_unit, which = "both", ask = FALSE)
print(surfaces_fix_unit$cv_loss)
print(surfaces_fix_unit$att)


## ---- 6. rmse_curve()/rmse_curves(): RMSE vs design dimension ---------------

g_quick <- rmse_curves(
  values_control = seq(20, 45, by = 5),
  values_pre     = seq(6, 30,  by = 4),
  n_runs  = 3,
  methods = methods,
  control = trop_control(n_cv_cells = 15L, cv_cycles = 1L, max_iter = 80L, seed = 3L)
)

print(g_quick)

# Two separate figures, then the side-by-side paper-style view (needs patchwork).
plot(g_quick)
plot(g_quick, combined = TRUE)

# Single dimension directly.
cc <- rmse_curve(
  "n_control",
  values  = seq(20, 45, by = 5),
  n_runs  = 3,
  methods = methods,
  control = trop_control(n_cv_cells = 15L, cv_cycles = 1L, max_iter = 80L, seed = 4L)
)

print(cc)
print(autoplot(cc))
print(head(as.data.frame(g_quick$n_control)))

# Paper-quality run template, intentionally commented out.
# g <- rmse_curves(values_control = seq(20, 45, by = 5), n_runs = 500,
#                  methods = methods, control = ctrl)
# plot(g, combined = TRUE, file = "rmse_curves.png", width = 12, height = 5)

p_heat_fix_unit <- autoplot(grid_fix_unit) +
  ggplot2::labs(
    title = "TROP sensitivity: fixed unit penalty",
    subtitle = sprintf("lambda_unit fixed at %.3f", fixed_lambda_unit)
  )
print(p_heat_fix_unit)


## ---- 7. Tidiers, anchors, ablation -----------------------------------------

# 7a. Counterfactual path plot from panel_compare().
p_cf <- plot_counterfactual(
  cmp,
  methods = intersect(c("DID", "MC", "TROP"), names(cmp$counterfactual))
)
print(p_cf)
# ggplot2::ggsave("counterfactual_paths.png", p_cf, width = 8, height = 5, dpi = 150)

# 7b. as_att() tidiers + bind_att() to stack heterogeneous results.
# bind_att() replaces the old rbind(as.data.frame(...)) + class<- dance: it
# coerces each argument via as_att(), aligns columns and keeps every row, so a
# multi-row panel_compare() result and single external fits combine cleanly.
att_trop_direct <- as_att(fit, method = "TROP_direct")

external_fit <- structure(
  list(estimate = fit$estimate + 0.10),
  class = "external_fit"
)

att_external <- as_att(external_fit, method = "external_demo", outcome = outcome)

att_aug <- bind_att(cmp, att_trop_direct, att_external)
print(att_aug)

p_att_aug <- autoplot(att_aug)
print(p_att_aug)

# 7c. Anchor comparison: paper-faithful per-cell vs fast pooled.
# Named arguments to bind_att() become the method labels.
fit_per_cell <- trop(
  df, outcome, treatment, unit, time,
  anchor = "per_cell", se = "none", control = ctrl_fast
)

fit_pooled <- trop(
  df, outcome, treatment, unit, time,
  anchor = "pooled", se = "none", control = ctrl_fast
)

anchor_tbl <- bind_att(TROP_per_cell = fit_per_cell, TROP_pooled = fit_pooled)
print(anchor_tbl)

p_anchor <- autoplot(anchor_tbl)
print(p_anchor)

# 7d. Penalty ablation: the Table-5-style robustness check.
# trop_ablation() replaces the old hand-built variants list + trop() loop: it
# chooses the full-spec penalties once by CV, then refits each constrained
# specification (no regression / no unit weights / no time weights / MC-like /
# DID-like) reusing those penalties, and returns a print/format-able table.
ablate_tbl <- trop_ablation(
  df, outcome, treatment, unit, time,
  anchor  = "pooled",
  se      = "none",
  control = ctrl_fast
)

print(ablate_tbl)
# Paste-ready Markdown / LaTeX (booktabs):
cat(format(ablate_tbl, "markdown"), sep = "\n")
# cat(format(ablate_tbl, "latex"), sep = "\n")

# Finite nuclear-norm penalty reused by the custom-grid demo in section 11.
lam_cv <- fit$lambda
nn_finite <- if (is.finite(lam_cv$nn)) lam_cv$nn else stats::sd(df[[outcome]], na.rm = TRUE)


## ---- 8. Inference modes on one fit: bootstrap / jackknife / placebo --------

# compare_se_modes() replaces the old separate fits + rbind: the penalties are
# cross-validated once and reused, so the point estimate is identical across
# rows and only the standard error / CI differ between resampling schemes.
ctrl_boot <- trop_control(
  n_cv_cells = 20L, cv_cycles = 1L, max_iter = 120L,
  n_boot = 30L, boot_ci = "percentile", seed = 2L
)

inf_tbl <- compare_se_modes(
  df, outcome, treatment, unit, time,
  se      = c("bootstrap", "jackknife", "placebo"),
  anchor  = "pooled",
  control = ctrl_boot
)

print(inf_tbl)

p_inf <- autoplot(inf_tbl)
print(p_inf)


## ---- 9. sim_semisynthetic() and trop_matrix() ------------------------------

# 9a. Matrix-in TROP interface. panel_matrices() reshapes the long panel into
# the Y/W matrices (plus treated rows and post-block width); trop_matrix() then
# infers treated_units / treated_periods from W, so the two compose directly.
pm <- panel_matrices(df, outcome, treatment, unit, time)
cat(sprintf("panel_matrices(): %d units x %d periods; %d treated units, %d post periods\n",
            nrow(pm$Y), ncol(pm$Y), length(pm$treated_units), pm$treated_periods))

tm_att <- trop_matrix(
  pm$Y, pm$W,
  lambda_unit = 0.1, lambda_time = 0.1, lambda_nn = Inf,
  control = ctrl_fast
)

print(tm_att)

# 9b. Semi-synthetic panel from a real-like untreated base panel.
real_base <- sim_panel(
  N = 35, T = 18, n_treated = 0L, t0 = 14,
  rank = 3, att = 0, noise = 1, seed = 3
)

# Dynamic effect path; post length is 18 - 14 + 1 = 5.
ss_obs <- sim_semisynthetic(
  real_base, "y", "id", "t",
  n_treated = 5, t0 = 14,
  effect = c(0, 1, 2, 2, 3),
  baseline = "observed", seed = 4
)

ss_lr <- sim_semisynthetic(
  real_base, "y", "id", "t",
  n_treated = 5, t0 = 14, att = 2,
  baseline = "lowrank", lambda_nn = 0.5, noise = 0.1, seed = 5
)

# true_att() reads the per-cell truth (the tau column, here) and averages over
# treated cells -- the number every estimator is trying to recover.
print(c(
  true_ATT_observed_baseline = true_att(ss_obs),
  true_ATT_lowrank_baseline  = true_att(ss_lr)
))

cmp_ss <- panel_compare(
  ss_obs, "y", "w", "id", "t",
  methods = methods,
  anchor  = "pooled",
  se      = "none",
  control = ctrl_fast
)

print(cmp_ss)

rmse_ss <- panel_rmse(
  ss_obs, "y", "w", "id", "t",
  methods  = methods,
  metric   = "placebo",
  horizon  = 3,
  n_pseudo = 4,
  n_runs   = 5,
  control  = ctrl_fast,
  seed     = 6
)

print(rmse_ss)

p_rmse_ss <- autoplot(rmse_ss)
print(p_rmse_ss)


## ---- 10. Optional engines and non-block treatment designs ------------------

# 10a. Optional engines on the original block design. Missing packages are
# skipped automatically by panel_compare().
opt_methods <- c("DID", "MC", "TROP", "DIFP", "SDID", "SC",
                 "gsynth", "augsynth", "CS")

cmp_optional <- panel_compare(
  df, outcome, treatment, unit, time,
  methods = opt_methods,
  anchor  = "pooled",
  se      = "none",
  control = ctrl_fast,
  verbose = TRUE
)

print(cmp_optional)

rmse_optional <- tryCatch(
  panel_rmse(
    df, outcome, treatment, unit, time,
    methods  = opt_methods,
    metric   = "placebo",
    horizon  = 3,
    n_pseudo = 4,
    n_runs   = 5,
    control  = ctrl_fast,
    seed     = 11,
    verbose  = TRUE
  ),
  error = function(e) {
    message("Skipping optional-engine RMSE because of error: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(rmse_optional)) {
  print(rmse_optional)
  print(autoplot(rmse_optional))
}

# 10b. Staggered absorbing design.
base_general <- sim_panel(
  N = 36, T = 18, n_treated = 0L, t0 = 14,
  rank = 3, att = 0, noise = 1, seed = 10
)

df_stag <- base_general
df_stag$w <- 0
df_stag$w[df_stag$id %in% 1:3 & df_stag$t >= 9]  <- 1
df_stag$w[df_stag$id %in% 4:6 & df_stag$t >= 12] <- 1
df_stag$w[df_stag$id %in% 7:9 & df_stag$t >= 15] <- 1
df_stag$y <- df_stag$y0 + 2 * df_stag$w

cmp_stag <- panel_compare(
  df_stag, "y", "w", "id", "t",
  methods = c("DID", "MC", "TROP", "CS"),
  anchor  = "pooled",
  se      = "none",
  control = ctrl_fast,
  verbose = TRUE
)

print(cmp_stag)
print(cmp_stag$pattern)

# Event study on staggered design.
fit_stag_trop <- trop(
  df_stag, "y", "w", "id", "t",
  anchor = "pooled", se = "none", control = ctrl_fast
)

es_stag <- trop_event_study(
  fit_stag_trop, se = "none", pre_periods = TRUE, control = ctrl_fast
)

print(es_stag)

p_es_stag <- autoplot(es_stag)
print(p_es_stag)

# 10c. Non-absorbing / on-off treatment design.
df_onoff <- base_general
df_onoff$w <- as.numeric(
  (df_onoff$id %in% 1:3 & df_onoff$t %in% c(8, 9, 14)) |
    (df_onoff$id %in% 4:6 & df_onoff$t %in% c(11, 12, 13))
)
df_onoff$y <- df_onoff$y0 + 1.5 * df_onoff$w

# Include SDID/SC/DIFP intentionally to verify graceful skips on non-block designs.
general_methods <- c("DID", "MC", "TROP", "DIFP", "SDID", "SC")

cmp_onoff <- panel_compare(
  df_onoff, "y", "w", "id", "t",
  methods = general_methods,
  anchor  = "auto",
  se      = "none",
  control = ctrl_fast,
  verbose = TRUE
)

print(cmp_onoff)
print(cmp_onoff$att)
print(cmp_onoff$pattern)


## ---- 11. Custom grids, full SVD, optional parallel workers -----------------

# 11a. Custom penalty grid.
custom_grids <- list(
  time = c(0, 0.05, 0.2),
  unit = c(0, 0.25, 0.75),
  nn   = c(Inf, nn_finite)
)

fit_custom_grid <- trop(
  df, outcome, treatment, unit, time,
  grids   = custom_grids,
  anchor  = "pooled",
  se      = "none",
  control = ctrl_fast,
  verbose = TRUE
)

print(fit_custom_grid)
print(fit_custom_grid$lambda)

# 11b. Exact full SVD mode.
fit_full_svd <- trop(
  df, outcome, treatment, unit, time,
  anchor  = "pooled",
  se      = "none",
  control = trop_control(
    n_cv_cells = 15L, cv_cycles = 1L, max_iter = 80L, svd = "full", seed = 12L
  )
)

print(fit_full_svd)

# 11c. Optional workers smoke test. trop_control(workers=) falls back to
# sequential automatically if future/future.apply are not installed.
ctrl_parallel <- trop_control(
  n_cv_cells = 15L, cv_cycles = 1L, max_iter = 80L, workers = 2L, seed = 13L
)

rmse_parallel <- panel_rmse(
  df, outcome, treatment, unit, time,
  methods  = c("DID", "TROP"),
  metric   = "placebo",
  horizon  = 2,
  n_pseudo = 3,
  n_runs   = 4,
  control  = ctrl_parallel,
  seed     = 13,
  verbose  = TRUE
)

print(rmse_parallel)


## ---- 12. RMSE curves over N_treated / N_post (built-in) --------------------

# These design dimensions are now built into rmse_curve(): the old hand-rolled
# quick_design_curve() is gone. vary = "n_treated" / "n_post" sweep the number
# of treated units and post-treatment periods respectively.
fig3_ntr <- rmse_curve(
  "n_treated",
  values   = 1:4,
  n_runs   = 2,
  methods  = c("DID", "MC", "TROP"),
  n_control = 30L, n_pre = 12L, n_post = 5L,
  att      = 2,
  control  = ctrl_fast,
  seed     = 21
)

fig3_ttr <- rmse_curve(
  "n_post",
  values   = 1:4,
  n_runs   = 2,
  methods  = c("DID", "MC", "TROP"),
  n_control = 30L, n_treated = 5L, n_pre = 12L,
  att      = 2,
  control  = ctrl_fast,
  seed     = 22
)

print(fig3_ntr)
print(fig3_ttr)

print(autoplot(fig3_ntr))
print(autoplot(fig3_ttr))


## ---- 13. Session info ------------------------------------------------------

utils::sessionInfo()
