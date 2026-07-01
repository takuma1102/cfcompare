# Regenerate every image in man/figures/ from the CURRENT package.
#
# man/figures/ is the single source of truth for images: the README embeds them
# with <img src="man/figures/...">, the help pages via \figure{...}, and the
# vignette via knitr::include_graphics("../man/figures/..."). Change a plot or the
# API, rerun this script, and the README updates immediately (it just links the
# file); rebuild the vignette / pkgdown to refresh those (cheap -- no estimator
# code runs at build, only the image is re-embedded).
#
# Usage (from the package root):
#   devtools::load_all(".")          # or library(cfcompare) on the installed pkg
#   source("data-raw/make-figures.R")
#
# Requires ggplot2. The panel_compare()/panel_rmse() SDID & SC rows also need
# 'synthdid'; they are dropped automatically if it is not installed. This script
# is excluded from the build via .Rbuildignore (^data-raw$).

library(cfcompare)
library(ggplot2)

fig_dir <- "man/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
save_gg <- function(name, plot, width = 6, height = 3.5, dpi = 150) {
  ggsave(file.path(fig_dir, name), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  message("wrote ", file.path(fig_dir, name))
}

set.seed(1)
df  <- sim_panel(N = 16, T = 9, n_treated = 3, t0 = 7, rank = 2,
                 att = 2, noise = 1, seed = 1)
ctl <- trop_control(n_cv_cells = 12L, cv_cycles = 1L, seed = 1)

## a single TROP fit reused by several figures
fit <- trop(df, "y", "w", "id", "t", se = "jackknife", control = ctl)

## 1. method-comparison forest (compare_forest.png)
cmp <- panel_compare(df, "y", "w", "id", "t",
                     methods = c("DID", "MC", "TROP"),
                     se = "jackknife", control = ctl)
save_gg("compare_forest.png", autoplot(cmp))

## 2. SC-style trajectory for a single fit (trop_trajectory.png)
save_gg("trop_trajectory.png", autoplot(fit), height = 4.5)

## 3. event study (trop_event_study_bootstrap.png)
es <- trop_event_study(fit, se = "bootstrap")
save_gg("trop_event_study_bootstrap.png", autoplot(es), height = 4.5)

## 4. placebo (randomization) test (placebo_test.png)
set.seed(1)
pt <- trop_placebo_test(fit, B = 500L)
save_gg("placebo_test.png", autoplot(pt))

## 5. bootstrap vs jackknife SE comparison (se_comparison.png)
se_cmp <- compare_se_modes(df, "y", "w", "id", "t",
                           se = c("bootstrap", "jackknife"), control = ctl)
save_gg("se_comparison.png", autoplot(se_cmp))

## 6. prediction and placebo RMSE bars (rmse_prediction.png, rmse_placebo.png)
save_gg("rmse_prediction.png",
        autoplot(panel_rmse(df, "y", "w", "id", "t",
                            methods = c("DID", "MC", "TROP", "DIFP"),
                            metric = "prediction", horizon = 3,
                            n_pseudo = 5, n_runs = 5, seed = 1)))
save_gg("rmse_placebo.png",
        autoplot(panel_rmse(df, "y", "w", "id", "t",
                            methods = c("DID", "MC", "TROP", "DIFP"),
                            metric = "placebo", horizon = 3,
                            n_pseudo = 5, n_runs = 5, seed = 1)))

## 7. penalty ablation, paper Table-5 figure (trop_ablation.png)
##    plot.trop_ablation() writes straight to `file`.
ab <- trop_ablation(df, "y", "w", "id", "t",
                    rmse = TRUE, horizon = 3, n_pseudo = 5, n_runs = 100L,
                    control = ctl)
plot(ab, file = file.path(fig_dir, "trop_ablation.png"))
message("wrote ", file.path(fig_dir, "trop_ablation.png"))

## 8. sensitivity heatmap + CV / ATT surfaces (a lambda sweep -> cf_trop_grid)
grid <- trop_sensitivity(df, "y", "w", "id", "t",
                         axes = c("nn", "time"), control = ctl)
save_gg("trop_sensitivity.png", autoplot(grid), height = 4.5)

png(file.path(fig_dir, "cv_surface.png"), width = 900, height = 600, res = 150)
plot_trop_surfaces(grid, which = "cv_loss", ask = FALSE)
dev.off()
png(file.path(fig_dir, "att_surface.png"), width = 900, height = 600, res = 150)
plot_trop_surfaces(grid, which = "att", ask = FALSE)
dev.off()

message("man/figures/ regenerated.")
