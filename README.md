# cfcompare

<!-- badges: start -->
[![R-CMD-check](https://github.com/takuma1102/nonabsdid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/takuma1102/nonabsdid/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

> **Note: ** `cfcompare` is an independent R package
> and not endorsed or maintained by the authors of the TROP estimator (Athey, Imbens, Qu &
> Viviano, 2026). TROP software written by the authors
> includes the Python package
> [`trop`](https://pypi.org/project/trop/)
> ([ostasovskyi/TROP-Estimator](https://github.com/ostasovskyi/TROP-Estimator))
> and a Stata implementation
> ([justinwaddy/TROP](https://github.com/justinwaddy/TROP)). Please cite the paper (see
> [Citation](#citation)) when using this software.

`cfcompare` is a R package that implements TROP and
places it in a common comparison workflow with other ATT estimators, such as difference-in-differences (DID),
synthetic difference-in-differences (SDID), synthetic control (SC), matrix
completion (MC), `gsynth`, `augsynth`, and DIFP (Doudchenko–Imbens–Ferman–Pinto)
on the same tidy schema and shared plots, so that applied researchers can run them on the same data and compare. The TROP, DID,
MC and DIFP engines are written natively in base R from the equations in the
paper (Athey, Imbens, Qu & Viviano, 2026,
[doi:10.1002/jae.70061](https://doi.org/10.1002/jae.70061); working paper
[arXiv:2508.21536](https://arxiv.org/abs/2508.21536)); the others wrap existing R
packages. Existing official TROP software includes the Python package `trop` and
a Stata implementation ([justinwaddy/TROP](https://github.com/justinwaddy/TROP));
`cfcompare` is not affiliated with or
endorsed by those authors or maintainers.

*Estimation RMSE of the ATT on a toy factor-model panel where treatment is
selected on unobserved factor loadings. Two-way fixed effects (DID) is biased by
the confounding; the factor- and weight-aware estimators recover the effect more
accurately. (Illustrative Monte Carlo; relative performance varies by setting,
as emphasised in the paper.)*

Supported methods, on a common `ATT` schema:

- **DID** — two-way fixed effects / difference-in-differences (native).
- **MC** — matrix completion with nuclear-norm regularisation (native).
- **DIFP** — Doudchenko–Imbens–Ferman–Pinto: synthetic control after recentering
  the mean (a demeaned SC with an intercept), written natively. This is the DIFP
  comparison estimator used in the paper's simulations.
- **TROP** — the triply robust panel estimator: low-rank + two-way FE outcome
  model, exponential-decay **unit** weights and **time** weights, with penalties
  chosen by leave-one-out cross-validation (native).
- **SDID** / **SC** — synthetic difference-in-differences and synthetic control,
  via [`synthdid`](https://github.com/synth-inference/synthdid).
- **gsynth** — an alternative matrix-completion / interactive-fixed-effects
  estimator, via [`gsynth`](https://cran.r-project.org/package=gsynth).
- **augsynth** — the augmented synthetic control method, via
  [`augsynth`](https://github.com/ebenmichael/augsynth).
- **CS** — the Callaway & Sant'Anna staggered DID estimator, via
  [`did`](https://cran.r-project.org/package=did).

The native engines (DID, MC, TROP, DIFP) depend only on base R, so the core
comparison always runs. Optional engines whose package is missing — or that do
not apply to the design — are skipped with a message, and the remaining methods
still produce output.

By default `panel_compare()` runs `DID`, `SDID`, `MC` and `TROP`. **DIFP** and
**CS** are available but opt-in: add them via `methods = `. CS (Callaway &
Sant'Anna) in particular is *not* part of the TROP paper's comparison set and
requires an absorbing (staggered/block) treatment, so it is offered as an
option rather than a default.

## Why TROP

The TROP working model writes the untreated potential outcome as

```
Y_it(0) = alpha_i + beta_t + L_it + eps_it
```

and estimates the treatment effect on a treated cell `(i, t)` by

```
minimise  sum_{j,s} theta_s * omega_j * (1 - W_js) (Y_js - a_j - b_s - L_js)^2  +  lambda_nn ||L||_*
then      tau_it = Y_it - a_i - b_t - L_it
```

with `theta_s = exp(-lambda_time |t - s|)` (time weights) and
`omega_j = exp(-lambda_unit * dist(j, i))` (unit weights). This single program
**nests** DID (`lambda_nn = Inf`, uniform weights), matrix completion (uniform
weights, finite `lambda_nn`), and synthetic-control-type weighting, which is why
it makes a natural backbone for an estimator-comparison package.

## Installation

```r
# install.packages("pak")
pak::pak("takuma1102/cfcompare")
```

`synthdid`, `gsynth` and `fixest` are in `Suggests` — install the ones you plan
to use.

## Quick start

```r
library(cfcompare)

# simulate a factor-model panel with a block treatment
df <- sim_panel(N = 30, T = 16, n_treated = 5, t0 = 12, att = 2, seed = 1)

# compare estimators on the same data
cmp <- panel_compare(
  df,
  outcome = "y", treatment = "w", unit = "id", time = "t",
  methods = c("DID", "SDID", "MC", "TROP"),
  se = "jackknife"
)

cmp$att          # tidy ATT table (one row per method)
autoplot(cmp)    # forest plot of ATTs with CIs
plot_counterfactual(cmp)  # observed vs predicted control paths
```

## Out-of-sample RMSE comparison

Beyond comparing point estimates, you can compare estimators by how well each
*predicts held-out outcomes* — the cross-model RMSE validation from the
doubly/triply robust panel estimator paper. `panel_rmse()` runs the paper's
"random blocks" placebo procedure: it repeatedly treats a random set of control
units as a placebo cohort, holds out their final periods, predicts those cells
with each method (treating real treated cells and held-out cells as unobserved),
and records the root-mean-square error.

```r
r <- panel_rmse(
  df, outcome = "y", treatment = "w", unit = "id", time = "t",
  methods = c("DID", "MC", "TROP"),
  horizon = 4, n_pseudo = 8, n_runs = 10, seed = 1
)
r              # ranked tidy table: method, rmse, rmse_se, ...
autoplot(r)    # ranked bar chart, lowest RMSE first (best method highlighted)
```

Lower RMSE is better, and the comparison is the point of a triply/doubly robust
estimator: it should sit at or below the best of DID, SC, SDID and MC. Native
methods (DID, MC, TROP) are always computed; `SDID`/`SC` (via `synthdid`) and
`gsynth` are included when their package is installed and the design suits them,
and otherwise skipped with a note. The procedure refits per placebo run, so cap
`n_runs`/`n_cv_cells` (via `trop_control()`) on large panels.

## Running a single estimator

```r
fit <- trop(df, "y", "w", "id", "t", se = "jackknife")
fit
fit$tau_cells        # per treated-cell effects
fit$counterfactual   # estimated Y(0) matrix
fit$lambda           # cross-validated penalties
```

Fix the penalties to recover special cases directly:

```r
did <- trop(df, "y", "w", "id", "t", lambda = list(time = 0, unit = 0, nn = Inf))
mc  <- trop(df, "y", "w", "id", "t", lambda = list(time = 0, unit = 0, nn = 5))
```

## Working from existing results

Already ran an estimator elsewhere? Coerce it into the shared schema and drop it
into the same comparison and plots:

```r
library(synthdid)
sd <- synthdid_estimate(Y, N0, T0)
tidy <- as_att(sd, method = "SDID", outcome = "y")

# combine several
all <- as_att(list(fit, sd))
autoplot(all)
```

## Tidy schema

`panel_compare()` and `as_att()` return a `cf_att_tbl` with:

| column            | type | description                                  |
| ----------------- | ---- | -------------------------------------------- |
| `method`          | chr  | `"DID"`, `"SDID"`, `"SC"`, `"MC"`, `"TROP"`, … |
| `estimate`        | num  | ATT point estimate.                          |
| `std.error`       | num  | Standard error (may be `NA`).                |
| `conf.low`        | num  | Lower CI bound.                              |
| `conf.high`       | num  | Upper CI bound.                              |
| `n_treated_cells` | int  | Number of treated unit-time cells.           |
| `n_treated_units` | int  | Number of treated units.                     |
| `outcome`         | chr  | Outcome variable name.                       |
| `engine`          | chr  | `"cfcompare"`, `"synthdid"`, `"gsynth"`, …       |
| `rank`            | int  | Estimated rank of `L` (native engines).      |
| `note`            | chr  | Skips / warnings.                            |

## Designs supported

The native engines (DID, MC, TROP) handle **block**, **staggered**, and
**non-absorbing** binary treatments — `treatment` is the 0/1 indicator of
*active* treatment in each cell. `SDID`/`SC` apply to block designs only and are
skipped otherwise.

## Inference

Native engines support three variance estimators, selected with `se=`:

- **`"bootstrap"`** — a unit-level *stratified block* bootstrap (resample treated
  and control units separately, with replacement, keeping each unit's full time
  series), the approach used for TROP's standard errors and confidence intervals
  in the paper. Percentile intervals by default; tune with
  `trop_control(n_boot=, boot_ci=)`.
- **`"jackknife"`** — leave-one-treated-unit-out (needs ≥ 2 treated units).
- **`"placebo"`** — placebo standard errors (for a single treated unit).

```r
fit <- trop(df, "y", "w", "id", "t", se = "bootstrap",
            control = trop_control(n_boot = 200))
```

These are practical approximations; see the paper for the estimator's formal
inference results.

## Semi-synthetic simulation

To benchmark estimators against a known truth while keeping a realistic data
structure, `sim_semisynthetic()` takes a real panel, uses its outcomes (directly,
or via a low-rank fit) as the untreated potential outcomes `Y(0)`, and imposes a
known effect on a placebo block — the semi-synthetic experiment style of the
paper.

```r
ss <- sim_semisynthetic(real_df, "y", "id", "t",
                        n_treated = 6, t0 = 14, att = 3, seed = 1)
mean(ss$tau[ss$w == 1])                 # true ATT = 3
panel_compare(ss, "y", "w", "id", "t")  # which estimator recovers it?
```

## Visualising a single fit

`autoplot()` on a `trop` fit gives a synthetic-control-style trajectory plot: the
treated-unit average against the estimated `Y(0)`, the post-treatment gap (the
effect) shaded, and the TROP time weights drawn as a ribbon along the bottom.

```r
autoplot(trop(df, "y", "w", "id", "t"))
```

## Penalty sensitivity (heatmap)

`trop_sensitivity()` sweeps the time penalty against the nuclear-norm penalty
(holding the unit penalty fixed) and records, at each grid point, both the ATT
estimate and the cross-validation loss. `autoplot()` draws the diagnostic
heatmap: cells are coloured by CV loss, annotated with the ATT, and the
CV-selected penalty pair is outlined — a quick read on how sensitive the
estimate is to the penalties and where the data-driven choice lands.

```r
g <- trop_sensitivity(df, "y", "w", "id", "t")
autoplot(g)
g   # prints the CV-selected penalties and the ATT range over the grid
```

## Numerical agreement with the official software

The TROP engine in `cfcompare` is written independently in base R from the
equations in the paper. As a correctness check only — not because any code was
shared or ported — its output was compared against the authors' official Python
package (PyPI [`trop`](https://pypi.org/project/trop/),
`ostasovskyi/TROP-Estimator`) on identical panels, and the numbers line up where
the two are exactly comparable:

- On the comparable special case (uniform weights, `lambda_nn = Inf`, i.e.
  weighted two-way fixed effects), the two agree to numerical tolerance (e.g.
  `3.87754` here vs `3.87751`).
- Computing the block-design weights from the paper's equation (3) inside R and
  feeding them in gives a weighted-TWFE ATT that matches to the printed precision
  (e.g. `3.90339` vs `3.90339`).
- The paper's equation (2) defines the ATT as `tau = avg(Y - alpha - beta - L)`
  over treated cells; that is what `cfcompare` computes directly, and it agrees
  with the official package's joint weighted-TWFE coefficient to `< 1e-4`.

`trop_matrix()` is a convenience matrix-in form
(`Y, W, treated_units, lambda_unit, lambda_time, lambda_nn, treated_periods`)
that makes such numerical comparisons against other matrix-based implementations
straightforward.

Two deliberate differences remain, by design. `cfcompare`'s default `trop()` uses
the paper's general unit/time distances (pairwise / pooled unit distance and
per-cell `|t-s|` time distance), so it applies beyond block designs. And the
finite-`lambda_nn` low-rank term is solved here by a proximal-gradient
soft-impute, so matrix-completion penalties are parameterised differently from a
convex-solver implementation (qualitative behaviour matches; exact digits do
not). For `SDID`/`SC`, `cfcompare` wraps the R
[`synthdid`](https://github.com/synth-inference/synthdid) package; `augsynth`,
`gsynth` and `did` (Callaway–Sant'Anna) are wrapped as optional engines.

## Status

Experimental. The output schema is intended to be stable, but the upstream
estimator packages occasionally rearrange their internal structures, so please
pin versions in production code.

## Citation

`cfcompare` is an unofficial implementation. If you use the TROP estimator, please
cite the paper:

- Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026). *Triply Robust Panel
  Estimators.* Journal of Applied Econometrics, 1–16.
  [doi:10.1002/jae.70061](https://doi.org/10.1002/jae.70061). Working-paper
  version: [arXiv:2508.21536](https://arxiv.org/abs/2508.21536).

and the underlying packages for any wrapped estimators (`synthdid`, `gsynth`,
`augsynth`). The authors' official software is the Python package
[`trop`](https://pypi.org/project/trop/).
