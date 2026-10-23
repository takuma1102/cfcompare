# Triply RObust Panel (TROP) estimator

Fits the TROP estimator of Athey, Imbens, Qu & Viviano (2026) on a long
panel. TROP combines a low-rank-plus-two-way-fixed-effects outcome model
with exponential-decay unit weights (upweighting controls similar to the
treated) and time weights (upweighting periods near the treated
periods). Penalties are chosen by leave-one-out cross-validation on the
control cells. The estimator nests DID/TWFE, matrix completion and
synthetic-control-type weighting as special cases.

## Usage

``` r
trop(
  data,
  outcome,
  treatment,
  unit,
  time,
  covariates = NULL,
  lambda = NULL,
  standardize = FALSE,
  anchor = c("auto", "per_cell", "pooled"),
  se = c("bootstrap", "auto", "jackknife", "none"),
  grids = NULL,
  control = trop_control(),
  verbose = FALSE
)
```

## Arguments

- data:

  A long `data.frame` with one row per unit-time.

- outcome, treatment, unit, time:

  Column names (strings). `treatment` must be a 0/1 indicator of
  *active* treatment in that cell (works for block, staggered, and
  non-absorbing designs).

- covariates:

  Optional character vector of column names for time-varying covariates.
  When supplied, the low-rank term is augmented additively as \\L\_{js}
  = X\_{js}'\phi + R\_{js}\\ (Athey, Imbens, Qu & Viviano 2025, Section
  6.2): the covariate-linear part is unpenalised and the nuclear norm
  applies to the residual \\R\\. Penalty cross-validation is performed
  on the covariate-residualised model. Covariates must be fully observed
  (including in treated cells, where they enter the counterfactual). The
  fitted coefficients are returned in `$phi`.

- lambda:

  Optional named list `list(time=, unit=, nn=)` to fix the penalties and
  skip cross-validation. When `standardize = TRUE`, supplied penalties
  are interpreted on the standardized-outcome scale. Whether `nn` is
  supplied here or chosen by cross-validation, the solver adds the small
  scale-relative stabilising floor `nn_floor` (see
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md))
  to every *finite* `nn` (the reported `$lambda$nn` stays nominal and
  the resolved floor is stored in `$lambda$nn_floor`); set
  `control = trop_control(nn_floor = 0)` to apply `nn` exactly as given.

- standardize:

  Logical; if `TRUE`, the outcome is standardized internally as
  `(Y - mean(Y)) / sd(Y)` (over all observed cells) before fitting, and
  the ATT, standard error, confidence interval, per-cell effects and
  counterfactual are mapped back to the raw outcome scale. This matches
  the convention of the official Stata command (v0.2.1+); as there, the
  selected `lambda` values are reported on the *standardized* scale,
  which makes them comparable across outcomes and across
  implementations. Default `FALSE` (raw outcome, the previous
  behaviour). The returned `panel$Y` is on the fitting (standardized)
  scale, with the centre/scale stored in `$scaling`.

- anchor:

  How weights are anchored to treated cells: `"per_cell"` re-solves
  eq. (2) with cell-specific weights for every treated cell (faithful to
  the paper); `"pooled"` solves once per *treatment spell group* –
  treated cells with the same start and end period, pooled across units,
  the `group(time)` estimands of the official Stata command – with
  weights anchored to that group (unit distances to the average
  trajectory of the group's units over the periods in which none of them
  is treated, time distances to the midpoint of the spell; see
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md)'s
  `time_center`). Simultaneous permanent adoption is one group,
  staggered adoption gives one group per cohort, and a unit that
  switches in and out of treatment contributes one group per spell; the
  ATT is the treated-cell-weighted average of the group effects, which
  are returned in `$tau_groups`. `"auto"` (default) uses `per_cell` when
  there are at most `max_cells` treated cells and `pooled` otherwise.
  [`trop_matrix()`](https://takuma1102.github.io/cfcompare/reference/trop_matrix.md)
  has no `anchor` argument; it is always pooled over a single trailing
  block.

- se:

  Standard-error method: `"bootstrap"` (default; unit-level block
  bootstrap, stratified by treatment status unless
  `trop_control(boot_resample = "controls")`), `"jackknife"`
  (delete-one-unit over *every* unit, treated and control, as the Stata
  command's `vce(jackknife)`; replicates that leave no treated unit are
  skipped), `"auto"` (jackknife when there are at least 2 treated units,
  else bootstrap), or `"none"`.

- grids:

  Optional list of penalty grids; see
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md).

- control:

  A list of solver/CV settings from
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md).

- verbose:

  Logical; print CV progress.

## Value

An object of class `trop`: a list with the ATT `estimate`, `std.error`,
`conf.low`/`conf.high`, the selected penalties `lambda`, any covariate
coefficients `phi`, per-cell effects (`tau_cells`), the per-spell-group
effects `tau_groups` (pooled anchor only: one row per group with its
start/end period, unit and cell counts and `estimate`; the ATT is their
`n_cells`-weighted average), the estimated counterfactual matrix, the
reshaped panel, and `scaling` (the outcome centre/scale used internally;
the identity unless `standardize = TRUE`). All effect and outcome-level
components are on the raw outcome scale; `panel$Y` and `lambda` are on
the fitting scale.

## References

Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026). Triply Robust
Panel Estimators. *Journal of Applied Econometrics*.
[doi:10.1002/jae.70061](https://doi.org/10.1002/jae.70061) .

## Examples

``` r
df <- sim_panel(N = 14, T = 9, n_treated = 3, t0 = 6, seed = 1)
# fixed penalties keep this example fast and deterministic; omit `lambda`
# to choose (time, unit, nn) by cross-validation
fit <- trop(df, "y", "w", "id", "t",
            lambda = list(time = 0.1, unit = 0.5, nn = 2), se = "none")
fit
#> Triply RObust Panel (TROP) estimator
#>   ATT = 0.5904
#>   penalties: lambda_time=0.1, lambda_unit=0.5, lambda_nn=2
#>   lambda_nn stabilising floor: +0.00466 (added in the solver; trop_control(nn_floor = 0) disables)
#>   estimated rank(L)=2, anchor=per_cell, design=block
#>   treated cells=12 across 3 unit(s)
```
