# *cfcompare*: compare TROP and panel ATT estimators in R

<!-- badges: start -->
[![R-CMD-check](https://github.com/takuma1102/cfcompare/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/takuma1102/cfcompare/actions/workflows/R-CMD-check.yaml)
[![r-universe version](https://takuma1102.r-universe.dev/cfcompare/badges/version)](https://takuma1102.r-universe.dev/cfcompare)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

> **Note**: `cfcompare` is an independent R package. It is not endorsed or
> maintained by the authors of the TROP article (Athey, Imbens, Qu & Viviano,
> 2026). Official TROP software by the authors includes the Python package
> [`trop`](https://pypi.org/project/trop/)
> ([ostasovskyi/TROP-Estimator](https://github.com/ostasovskyi/TROP-Estimator))
> and the Stata command
> ([justinwaddy/TROP](https://github.com/justinwaddy/TROP)). Please cite the
> paper when using the TROP estimator; see [Citation](#citation).

`cfcompare` implements the triply robust panel estimator (TROP) and puts it in a
common R workflow with DID/TWFE, matrix completion, synthetic DID, synthetic
control, DIFP, `gsynth`, `augsynth`, and Callaway--Sant'Anna DID. The point is
not just to estimate TROP, but to compare estimators on the same panel, schema,
plots, counterfactual matrices, RMSE diagnostics, and penalty-sensitivity tools.

Use `cfcompare` when you want to:

- run several ATT estimators on one binary-treatment panel with `panel_compare()`;
- get a shared tidy ATT table and shared plotting methods;
- evaluate estimators by held-out RMSE with `panel_rmse()`;
- inspect TROP tuning through lambda-grid CV-loss and ATT surfaces.

## Installation

```r
# install.packages("pak")
pak::pak("takuma1102/cfcompare")
```

You can also install from R-universe.

```r
install.packages(
  "cfcompare",
  repos = c("https://takuma1102.r-universe.dev", "https://cloud.r-project.org")
)
```

The core native engines run without optional estimator packages. Install
`synthdid`, `gsynth`, `augsynth`, or `did` only for the wrapped methods you plan
to use.

## Quick start

```r
library(cfcompare)

# Compare multiple estimators.
cmp <- panel_compare(
  df,
  outcome = "y", treatment = "w", unit = "id", time = "t",
  methods = c("DID", "SDID", "MC", "TROP", "DIFP"),
  se = "bootstrap"
)

cmp$att                    # tidy ATT table, one row per method
autoplot(cmp)              # forest plot of ATT estimates and intervals
```

<img src="man/figures/compare_forest.png" alt="Forest plot of ATT estimators" />

```r
plot_counterfactual(cmp)   # observed vs fitted untreated paths
```

## Supported estimators

| Method | Engine | `panel_compare()` default? | Notes |
| --- | --- | --- | --- |
| `DID` | native R | yes | two-way fixed effects / difference-in-differences |
| `MC` | native R | yes | nuclear-norm matrix completion |
| `DIFP` | native R | yes | Doudchenko-Imbens-Ferman-Pinto demeaned SC with intercept |
| `TROP` | native R | yes | low-rank + two-way FE outcome model with unit/time weights |
| `SDID` | `synthdid` | yes | skipped if `synthdid` is unavailable or the design is unsupported |
| `SC` | `synthdid` | no | opt in with `methods =` |
| `gsynth` | `gsynth` | no | optional interactive-fixed-effects / MC-style engine |
| `augsynth` | `augsynth` | no | optional augmented synthetic control engine |
| `CS` | `did` | no | Callaway--Sant'Anna; requires an absorbing staggered/block treatment |

By default, `panel_compare()` runs `DID`, `SDID`, `MC`, `TROP`, and `DIFP`.
Use `methods =` to specify an explicit set, or `exclude =` to drop one method
from the default set. Optional engines whose package is missing, or whose design
requirements are not met, are skipped with a message while the remaining methods
still run.

## What cfcompare adds around TROP

The official Python and Stata projects are the reference software for estimating
TROP. `cfcompare` adds an R-native comparison and diagnostic layer around TROP:

| Goal | Entry point | Output |
| --- | --- | --- |
| Compare multiple ATT estimators | `panel_compare()` | `cf_comparison` with a common `cf_att_tbl`, fits, panel data, and counterfactuals |
| Compare held-out predictive performance | `panel_rmse()` | ranked RMSE table; this is separate from TROP tuning cross-validation |
| Track estimation RMSE vs design size | `rmse_curve()` / `rmse_curves()` | `sqrt(E[(tau_hat - tau)^2])` vs N, T, #treated, or #post, with a known true ATT |
| Inspect TROP penalty sensitivity | `trop_sensitivity()` + `autoplot()` | lambda grid with CV loss and ATT at each grid point |
| Plot CV-loss / ATT surfaces | `plot_trop_surfaces()` | separate full-width CV-loss and ATT surface plots; returns surface matrices invisibly |
| Audit penalty components | `trop_ablation()` | table moving from full TROP toward MC and DID by constraining penalties |
| Compare inference choices | `compare_se_modes()` | one ATT estimate with bootstrap, jackknife, and/or placebo uncertainty rows |
| Reuse fitted counterfactuals | `counterfactual_matrix()` | common `N x T` estimated untreated-outcome matrix interface |

## Out-of-sample RMSE comparison

`panel_rmse()` compares estimators by how well they predict held-out outcomes.
The default `metric = "placebo"` repeatedly assigns a placebo block to control
units and scores the resulting zero-effect placebo ATT. `metric = "prediction"`
scores per-cell counterfactual prediction error through `counterfactual_matrix()`.

```r
r <- panel_rmse(
  df, outcome = "y", treatment = "w", unit = "id", time = "t",
  methods = c("DID", "MC", "TROP", "DIFP"),
  horizon = 4, n_pseudo = 8, n_runs = 10, seed = 1
)

r            # ranked table: method, rmse, rmse_se, engine, note
autoplot(r)  # lower RMSE is better
```

For large panels, reduce `n_runs` and TROP CV work, for example through
`trop_control(n_cv_cells = , cv_cycles = )`.

> **Note**: This is a **predictive** error on held-out cells, not estimation error against a
> known effect. It is a different quantity from the *estimation* RMSE `sqrt(E[(tau_hat - tau)^2])` reported by `rmse_curve()` / `rmse_curves()` over
> Monte Carlo replications with a known true ATT. In short: `panel_rmse()` asks
> "how well does each method predict outcomes it never saw?", while `rmse_curve()`
> asks "how close is each method's ATT to the truth as the design changes?" — and the two can rank methods differently.

## TROP diagnostics

Run a single TROP fit when you need the selected penalties, per-treated-cell
effects, or the fitted untreated counterfactual matrix.

```r
fit <- trop(df, "y", "w", "id", "t")
fit$lambda                  # CV-selected (time, unit, nn) penalties
fit$tau_cells               # per-treated-cell effects
counterfactual_matrix(fit)  # fitted N x T untreated-outcome matrix
autoplot(fit)               # synthetic-control-style trajectory
```

Inspect the TROP penalty surface by sweeping any two of the three penalties
(`time`, `unit`, `nn`) and holding the third fixed.

```r
g <- trop_sensitivity(
  df, "y", "w", "id", "t",
  axes = c("nn", "time"),
  control = trop_control(n_cv_cells = 12L, cv_cycles = 1L)
)

autoplot(g)                    # compact ggplot2 heatmap
surfaces <- plot_trop_surfaces(g, which = "both", ask = FALSE)
surfaces$cv_loss               # matrix behind the CV-loss surface
surfaces$att                   # matrix behind the ATT surface
```

Other targeted diagnostics are available when needed:

```r
trop_ablation(df, "y", "w", "id", "t")
compare_se_modes(df, "y", "w", "id", "t", se = c("bootstrap", "jackknife"))
trop_event_study(fit, se = "bootstrap")
```

## Designs and inference

The treatment column is the active 0/1 treatment indicator for each unit-time
cell. The native engines handle block, staggered, and non-absorbing treatment
patterns. `SDID` and `SC` are block-design methods and are skipped otherwise.
`CS` is opt-in and requires absorbing staggered/block treatment.

Native estimators support `se = "bootstrap"`, `"jackknife"`, and
`"none"` where applicable. These are practical resampling choices for the R
workflow; use the TROP paper for formal inference conditions.

## Numerical agreement and deliberate differences

The TROP engine is written independently in base R from the paper's equations.
For correctness checks, it has been compared with the official Python `trop`
package on exactly comparable weighted-TWFE special cases, using
`trop_control(svd = "full")`; the results agree to numerical tolerance.
`trop_matrix()` provides a matrix-in interface for direct cross-checks.

Two differences are deliberate. First, `cfcompare`'s default `trop()` uses the
paper's general unit/time distances, so it applies beyond simple block designs.
Second, finite nuclear-norm penalties are solved with a proximal-gradient
soft-impute routine, so exact digits need not match convex-solver implementations
outside the comparable special cases.

## More detail

The README is intentionally short. For details, use the package help topics:

```r
?panel_compare
?panel_rmse
?trop
?trop_sensitivity
?trop_ablation
?compare_se_modes
```

The source vignette is in [`vignettes/cfcompare.Rmd`](vignettes/cfcompare.Rmd);
this can later be surfaced through pkgdown.

## Status

Experimental. The output schema is intended to be stable, but wrapped estimator
packages can change their internal objects. Pin package versions in production
code.

## Citation

`cfcompare` is an unofficial implementation. If you use the TROP estimator,
please cite the paper:

- Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026). *Triply Robust Panel
  Estimators.* Journal of Applied Econometrics, 1--16.
  [doi:10.1002/jae.70061](https://doi.org/10.1002/jae.70061). Working-paper
  version: [arXiv:2508.21536](https://arxiv.org/abs/2508.21536).

Also cite the underlying packages for any wrapped estimators you use
(`synthdid`, `gsynth`, `augsynth`, and `did`).
