# *cfcompare*: R package implementing TROP and comparing it with other ATT estimators

<!-- badges: start -->
[![R-CMD-check](https://github.com/takuma1102/cfcompare/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/takuma1102/cfcompare/actions/workflows/R-CMD-check.yaml)
[![r-universe version](https://takuma1102.r-universe.dev/cfcompare/badges/version)](https://takuma1102.r-universe.dev/cfcompare)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
<!-- badges: end -->

> **Note**: `cfcompare` is an independent R package. This is not officially endorsed or
> maintained by the authors of the TROP article (Athey, Imbens, Qu & Viviano,
> 2026). Official TROP software by the authors includes the Python package
> [`trop`](https://pypi.org/project/trop/)
> ([ostasovskyi/TROP-Estimator](https://github.com/ostasovskyi/TROP-Estimator))
> and the Stata command
> ([justinwaddy/TROP](https://github.com/justinwaddy/TROP)). Please cite the original
> TROP paper when using the TROP estimator. See [Citation](#citation) for more details.

`cfcompare` implements the triply robust panel estimator (TROP) and puts it in a
common R workflow with DID/TWFE, synthetic control, matrix completion, synthetic
DID, and DIFP. This package enables not just to estimate TROP, but also to compare various ATT estimators on the same schema.

Use `cfcompare` when you want to:

- run several ATT estimators on one binary-treatment panel with `panel_compare()`;
- get a shared tidy ATT table and shared plotting methods;
- evaluate estimators by held-out RMSE with `panel_rmse()`;
- inspect TROP tuning through lambda grids with `trop_sensitivity` and `plot_trop_surfaces` (for CV-loss and ATT surfaces).

> Note: all plots and tables shown in this README are based on a model sample dataset, and they are provided solely as
> examples of the visualizations that can be created with this package.

## Installation

```r
# install.packages("pak")
pak::pak("takuma1102/cfcompare")
```

You can also install from R-universe. This package will be submitted to CRAN as well down the line.

```r
install.packages(
  "cfcompare",
  repos = c("https://takuma1102.r-universe.dev", "https://cloud.r-project.org")
)
```

Most of the core native engines run without optional estimator packages. Install
`synthdid`, `gsynth`, `augsynth`, or `did` for the wrapped methods you plan
to use.

## Quick start

```r
library(cfcompare)

# Compare multiple estimators.
cmp <- panel_compare(
  df,
  outcome = "y", treatment = "w", unit = "id", time = "t",
  methods = c("DID", "SDID", "SC", "MC", "DIFP", "TROP"),
  se = "bootstrap",  # You can also choose jackknife or "none" for running this faster.
  # covariates = c("x1", "x2"),
  anchor  = "pooled",  # faster method
  # or "per-cell" (more accurate; in accordance with the original paper)
  # or "auto" (choosing "per-cell" or "pooled" depending on the number of treated cells)
)

cmp$att                    # tidy ATT table, one row per method
autoplot(cmp)              # forest plot of ATT estimates and intervals (see below for an example plot)
```

<img src="man/figures/compare_forest.png" alt="Forest plot of ATT estimators" />


## Supported estimators

| Method | Engine | `panel_compare()` default? | Notes |
| --- | --- | --- | --- |
| `TROP` | native R | yes | low-rank + two-way FE outcome model with unit/time weights |
| `DID` | native R | yes | two-way fixed effects |
| `SC` | `synthdid` | yes | skipped if `synthdid` is unavailable or the design is unsupported |
| `MC` | native R | yes | nuclear-norm matrix completion |
| `SDID` | `synthdid` | yes | skipped if `synthdid` is unavailable or the design is unsupported |
| `DIFP` | native R | yes | Doudchenko-Imbens-Ferman-Pinto estimator: demeaned SC with intercept |
| `gsynth` | `gsynth` | no | optional interactive-fixed-effects / MC-style engine |
| `augsynth` | `augsynth` | no | optional augmented synthetic control engine |
| `CS` | `did` | no | Callaway--Sant'Anna; requires an absorbing staggered treatment |

By default, `panel_compare()` runs `TROP`, `DID`, `SC`, `MC`, `SDID`, and `DIFP`.
Use `methods =` to specify an explicit set, or `exclude =` to drop one method
from the default set. Optional engines whose package is missing, or whose design
requirements are not met (e.g., in terms of covariates), are skipped with a message while the remaining methods
still run.

## Overview: what cfcompare can do

`cfcompare` has basic functions for estimating the TROP estimator (the official Python and Stata projects are referenced). It also adds an R-native comparison and diagnostic layer around TROP:

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

## Held-out RMSE comparison

`panel_rmse()` compares estimators by how well they predict held-out outcomes.
The default `metric = "placebo"` repeatedly assigns a placebo block to control
units and scores the resulting zero-effect placebo ATT.

```r
r <- panel_rmse(
  df, outcome = "y", treatment = "w", unit = "id", time = "t",
  methods = c("DID", "SDID", "SC", "MC", "DIFP", "TROP"),
  horizon = 4, n_pseudo = 6, n_runs = 100, seed = 1
)

r            # ranked table: method, rmse, rmse_se, engine, note
autoplot(r)  # see below for an example plot
```
<img src="man/figures/rmse_placebo.png" alt="Placebo RMSE" />

 `metric = "prediction"` scores per-cell counterfactual prediction error.
 
<img src="man/figures/rmse_prediction.png" alt="Prediction RMSE" />

For quick diagnostics on large panels, reduce the outer repetition counts, such as `n_runs`, and use lighter TROP controls via `trop_control(n_cv_cells = , cv_cycles = , max_iter = , workers = )`.

> **Note**: This is a **predictive** error on held-out cells, not estimation error against a
> known effect. It is a different quantity from the *estimation* RMSE `sqrt(E[(tau_hat - tau)^2])` reported
> by `rmse_curves()` and `rmse_curve()` (for a single-dimension plot) over
> Monte Carlo replications with a known true ATT. These two can rank methods differently.

## TROP diagnostics

Run a single TROP fit when you need the selected penalties, per-treated-cell
effects, or the fitted untreated counterfactual matrix.

```r
fit <- trop(df, "y", "w", "id", "t")
fit$lambda                  # CV-selected (time, unit, nn) penalties
fit$tau_cells               # per-treated-cell effects
counterfactual_matrix(fit)  # fitted N x T untreated-outcome matrix
autoplot(fit)               # synthetic-control-style trajectory; see below
```
<img src="man/figures/trop_trajectory.png" alt="SC-style trajectory" />

Inspect the TROP penalty surface by sweeping any two of the three penalties
(`time`, `unit`, `nn`) and holding the third fixed.

```r
g <- trop_sensitivity(
  df, "y", "w", "id", "t",
  axes = c("unit", "time"), # You can also choose "nn" for lambda_nn
  anchor  = "pooled"  # 
)

autoplot(g)                    # compact ggplot2 heatmap (see below)
```
<img src="man/figures/trop_sensitivity.png" alt="Heatmap" />

```r
surfaces <- plot_trop_surfaces(g, which = "both", ask = FALSE) (see below)
surfaces$cv_loss               # matrix behind the CV-loss surface
surfaces$att                   # matrix behind the ATT surface
```

<img src="man/figures/cv_surface.png" alt="Heatmap for cv" />
<img src="man/figures/att_surface.png" alt="Heatmap for att" />

Other targeted diagnostics are available when needed.

Use `trop_ablation` when you want to check how RMSEs change in accordance with the change in lambda penalties. (See also Table 5 in the original paper.)
```r
trop_ablation(df, "y", "w", "id", "t")  # see below
```

<img src="man/figures/trop_ablation.png" alt="Ablation Table" />

Use `compare_se_modes` when you want to compare bootstrap and jackknife standard errors.
```r
compare_se_modes(df, "y", "w", "id", "t", se = c("bootstrap", "jackknife"))  # see below
```
<img src="man/figures/se_comparison.png" alt="Bootstrap v. jackknife" />

For event-study plots, use `trop_event_study`.

```r
trop_event_study(fit, se = "bootstrap")
```
<img src="man/figures/trop_event_study_bootstrap.png" alt="Event Study Plot" />

For placebo test relabeling control units as pseudo-treated, use `trop_placebo_test()` (penalties are held fixed).

```r
pt <- trop_placebo_test(
  fit, B = 500L, alternative = "two.sided"
)

autoplot(pt)    # histogram of the null; observed ATT (red solid line) and 0 (dashed); see below.
```

<img src="man/figures/placebo_test.png" alt="Placebo test result" />

## Designs and inference

The treatment column is the active 0/1 treatment indicator for each unit-time
cell. The TROP diagnosis handles block, staggered, and non-absorbing treatment
patterns under the same schema. By contrast, `SDID` and `SC` are block-design methods and are skipped otherwise.
`CS` is opt-in and requires absorbing staggered treatment.

Native estimators support `se = "bootstrap"`, `"jackknife"`, and
`"none"` where applicable. These are practical resampling choices for the R
workflow; use the TROP paper for formal inference conditions.

## Numerical agreement with, and technical differences from the official packages

The TROP engine in this package is written in base R.
For correctness checks, it has been compared with the official Python `trop`
package on exactly comparable weighted-TWFE sample datasets, using
`trop_control(svd = "full")`, and the results agree to numerical tolerance.

The below is how the original TROP paper, the official Python/Stata packages, and `cfcompare` handle main implementation choices. The first table is the comparison as of July 2026 (kept for reference); the second is the comparison as of October 2026, after the official packages' summer releases (Python `trop` v0.1.11, Stata `trop` v0.2.9) and the corresponding changes in `cfcompare` (0.1.0.9007).

**As of July 2026**

|  | TROP Paper | Python | Stata | R (cfcompare)  |
|:--|:--|:--|:--|:--|
| **Solver** | Silent | CVXPY using SCS or OSQP | FISTA (accelerated proximal gradient) | Soft-impute-style proximal gradient (base R only)|
| **Weight anchoring** | Per-cell  | Pooled only | Pooled (default) and per-cell (option) | Auto (choosing per-cell or pooled by treated-cell count) by default; per-cell and pooled as options |
| **CV criterion** | LOOCV | Placebo resampling (default) and k-fold (option) | Placebo resampling (default) and k-fold (option) for pooled; LOOCV for per-cell | LOOCV (default) and placebo resampling (option) |
| **Placebo pattern for staggered treatment** | N/A (LOOCV) | Common block | Actual per-unit patterns | Actual per-unit patterns; warning under LOOCV + staggered |
| **Outcome pre-processing** | Raw Y | Raw Y | Standardized (default) and raw Y (option) | Raw Y (default) and standardized (option) |
| **Covariates** | Yes | No | No | Yes |
| **Inference** | Bootstrap | None | Bootstrap | Bootstrap or jackknife (when available) |

**As of October 2026** (Python `trop` v0.1.11, Stata `trop` v0.2.9, cfcompare 0.1.0.9007)

|  | TROP Paper | Python (v0.1.11) | Stata (v0.2.9) | R (cfcompare 0.1.0.9007) |
|:--|:--|:--|:--|:--|
| **Solver** | Silent | CVXPY (SCS for finite `lambda_nn`, OSQP for `lambda_nn = Inf`), `sum_squares` formulation | FISTA (accelerated proximal gradient with adaptive restart), warm-start nuclear path | Soft-impute-style proximal gradient (base R only; optional truncated SVD via `RSpectra`), warm-start nuclear path |
| **Estimands / weight anchoring** | Per-cell | Pooled only, over a single trailing block (`treated_periods`) | `group(time)` (default): one pooled estimand per *treatment spell* (same start and end, pooled across units), cell-weighted average; `group(cell)`: per-cell | Auto (per-cell when at most `max_cells` treated cells, else pooled) by default; `anchor = "per_cell"` or `"pooled"` (per treatment spell, as Stata's `group(time)`; per-spell effects in `fit$tau_groups`) |
| **Pooled unit distance** | N/A (per-cell eq. (3)) | RMS distance to the average treated path over the pre-periods | RMS distance to the group's average path over the periods in which no group member is treated; each comparison unit's own treated cells masked (v0.2.6) | Same as Stata (`.unit_distance_to_avg()`); identical to Python on a block design |
| **Pooled time-weight centre** | N/A (per-cell `\|t - s\|`) | `T - tp/2` (half a period after the block midpoint) | Spell midpoint `(first + last)/2` (v0.2.4) | Spell midpoint by default; `trop_control(time_center = "python")` reproduces the Python convention, which `trop_matrix()` always uses |
| **Own-group distance** | N/A | Treated units compared with their own average (non-zero in general) | Same; `pooled_treat_distance(unit(#) time(#))` overrides the members' unit distance and the in-spell time distance (v0.2.6) | Same; `trop_control(pooled_treat_distance = list(unit =, time =))` |
| **Units with no usable comparison period** | Undefined (0/0) | Not applicable (block design: every unit shares the pre-periods) | Unit weight set to 0 when `lambda_unit > 0` (included when `lambda_unit = 0`) | Largest finite distance (smallest positive weight) |
| **Other treated cells in the fit** | Excluded | Excluded (block design) | Excluded; tau is the mean residual over the target cells (v0.2.9) | Excluded; tau is the mean residual over the target cells |
| **CV criterion** | LOOCV | Placebo resampling (default), k-fold; cycle, joint or adaptive search | Placebo resampling (default) or k-fold for `group(time)`; LOOCV for `group(cell)`; cycle, joint or adaptive search | LOOCV (default) or placebo resampling (`cv_method = "placebo"`); coordinate-descent (cycle) search |
| **Placebo pattern for staggered / general treatment** | N/A (LOOCV) | Common trailing block | Actual per-unit patterns, grouped by treatment spell | Actual per-unit patterns, grouped by treatment spell; warning under LOOCV + staggered |
| **Outcome pre-processing** | Raw Y | Raw Y | Standardized (always); lambdas on the standardized scale | Raw Y (default) and standardized (`standardize = TRUE`) |
| **Covariates** | Yes | No | No | Yes |
| **Bootstrap** | Stratified block bootstrap | `trop.bootstrap`: unit-level bootstrap, controls-only resampling (default) or stratified; percentile CI | `vce(bootstrap)`: stratified block bootstrap; percentile CI (no inference by default) | `se = "bootstrap"` (default): stratified (default) or controls-only (`boot_resample = "controls"`); percentile or normal CI |
| **Jackknife** | -- | -- | `vce(jackknife)`: delete-one-unit over every unit | `se = "jackknife"`: delete-one-unit over every unit (replicates that leave no treated unit are skipped) |

See the notes below for technical differences for this package's estimation from official packages.
> Note 1 (solver): finite nuclear-norm penalties are solved with a proximal-gradient
> soft-impute routine, as opposed to FISTA or SCS, so that the implementation depends only on base R. Thus, exact digits need not perfectly match convex-solver implementations.
> To reconcile this difference, you can fix the penalties via `lambda = list(time=, unit=, nn=)` (bypassing CV), set `svd = "full"`, and match its `anchor`.

> Note 2 (cv): penalties are chosen by leave-one-control-cell-out prediction error by default; a placebo criterion is also available via `trop_control(cv_method = "placebo")`.

> Note 3 (per-cell or pooled): this package supports both anchoring modes of per-cell (solving a separate local weighting problem for each treated cell) and pooled methods (one set of weights per *treatment spell group* -- treated cells with the same start and end period, pooled across units, the `group(time)` estimands of the Stata command; faster than per-cell). Under block adoption there is one group; under staggered adoption one per cohort; a unit that switches in and out of treatment contributes one group per spell. The group's unit distance is to the average path of its units over the periods in which none of them is treated, each comparison unit's own treated cells are masked, and the time weights are centred on the spell midpoint. The per-group effects are returned in `fit$tau_groups`; the ATT is their cell-weighted average. `trop_control(pooled_treat_distance = list(unit = 0, time = 0))` gives the group's own units and periods the full weight (the Stata command's `pooled_treat_distance()`). `trop_matrix()` is always pooled over a single trailing block (the Python signature).

> Note 4 (covariates): this package allows to use covariates through `panel_compare`.

> Note 5 (lambda_nn's floor): for numerical stability, the solver adds a small scale-relative stabilising floor to every finite `lambda_nn` by default. (Reported `lambda_nn` as part of results do not include the samll floor.) `lambda_nn = Inf` (i.e., the DID/TWFE-like case) is never affected. `trop_control(nn_floor = 0)` disables this adjustment so penalties are applied exactly as supplied.

> Note 6 (standardization): estimation uses the raw, non-standardized outcome by default, so that `lambda` values are on the outcome's natural scale. `trop(standardize = TRUE)` optionally standardizes the outcome.

> Note 7 (inference): `se = "bootstrap"` resamples units with replacement keeping each unit's full time series -- treated and control units separately by default (the stratified block bootstrap of the paper and of the Stata command), or only the control units with `trop_control(boot_resample = "controls")` (the default of the Python package's `trop.bootstrap`). `se = "jackknife"` deletes one unit at a time over every unit, as the Stata command's `vce(jackknife)`; replicates that leave no treated unit are skipped.


## More details

For more details, use the package help topics:

```r
?panel_compare
?panel_rmse
?trop
?trop_sensitivity
?trop_ablation
?compare_se_modes
```

The source vignette is now being developed and the current version is stored in [`vignettes/cfcompare.Rmd`](vignettes/cfcompare.Rmd);
this will later be surfaced through pkgdown after being finalized.

## Status

Experimental. The output schema is intended to be stable, but wrapped estimator
packages can change their internal objects. Pin package versions in production
code.

## Citation

`cfcompare` is an independent package to implement TROP. If you use the TROP estimator,
please cite the original JAM paper:

- Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026). *Triply Robust Panel
  Estimators.* Journal of Applied Econometrics, 1--16.
  [doi:10.1002/jae.70061](https://doi.org/10.1002/jae.70061).

Also cite the underlying packages for any wrapped estimators you use
(`synthdid`, `gsynth`, `augsynth`, and `did`).
