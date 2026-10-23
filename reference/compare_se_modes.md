# Compare TROP standard-error methods on one fit

Refits the TROP estimator under several inference methods — the
unit-level block bootstrap and/or the delete-one-unit jackknife (over
every unit) — and returns the estimate, standard error and confidence
interval from each on the common
[`cf_att_tbl`](https://takuma1102.github.io/cfcompare/reference/as_att.md)
schema. The penalties are chosen once (a single `se = "none"` fit — by
cross-validation, or fixed up front with `lambda`) and reused for every
inference method, exactly as
[`trop_ablation()`](https://takuma1102.github.io/cfcompare/reference/trop_ablation.md)
reuses one fit's penalties, so the ATT is identical across rows and only
the uncertainty changes. This is the table to reach for when checking
how sensitive the confidence interval is to the choice of resampling
scheme — it is an *inference-mode* comparison, not an estimator
comparison, so all rows share one point estimate by construction.

## Usage

``` r
compare_se_modes(
  data,
  outcome,
  treatment,
  unit,
  time,
  covariates = NULL,
  lambda = NULL,
  se = c("bootstrap", "jackknife"),
  anchor = "pooled",
  control = trop_control(),
  labels = NULL
)
```

## Arguments

- data, outcome, treatment, unit, time, covariates:

  Passed to
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md).

- lambda:

  Optional named list `list(time=, unit=, nn=)` to fix the penalties and
  skip cross-validation, exactly as in
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md).
  When `NULL` (the default) the penalties are chosen once by
  cross-validation; either way the resolved penalties are shared by
  every inference method.

- se:

  Character vector of inference methods to compare; any of
  `"bootstrap"`, `"jackknife"`. Defaults to both. (To compare anchors or
  penalty constraints instead, see
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)'s
  `anchor` argument and
  [`trop_ablation()`](https://takuma1102.github.io/cfcompare/reference/trop_ablation.md).)

- anchor:

  Weight anchoring, as in
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)
  (default `"pooled"`).

- control:

  A
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md)
  list, shared by every fit. Set `n_boot` here to control the bootstrap
  replications.

- labels:

  Optional `method` labels, one per `se`; defaults to `"Bootstrap SE"` /
  `"Jackknife SE"`.

## Value

An object of class `cf_se_comparison` (extending `cf_att_tbl`) with one
row per inference method, sharing a common point estimate. If a method
does not apply to the design (for example the bootstrap with a single
control unit) it still contributes a row, carrying the shared estimate
with `NA` uncertainty and an explanatory `note`.

## See also

[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md),
[`trop_ablation()`](https://takuma1102.github.io/cfcompare/reference/trop_ablation.md),
[`bind_att()`](https://takuma1102.github.io/cfcompare/reference/bind_att.md)

## Examples

``` r
df <- sim_panel(N = 16, T = 9, n_treated = 3, t0 = 6, att = 2, seed = 1)
# \donttest{
compare_se_modes(df, "y", "w", "id", "t",
                 control = trop_control(n_cv_cells = 5L, cv_cycles = 1L,
                                        tol = 1e-3, n_boot = 20L,
                                        seed = 1L))
#>        method estimate std.error conf.low conf.high n_treated_cells
#>  Bootstrap SE   2.5253    0.2278   2.0600    2.8166              12
#>  Jackknife SE   2.5253    0.2901   1.9567    3.0940              12
#>  n_treated_units outcome    engine rank note
#>                3       y cfcompare    1 <NA>
#>                3       y cfcompare    1 <NA>

# Fix the penalties (skip cross-validation), as with trop()'s `lambda`:
compare_se_modes(df, "y", "w", "id", "t",
                 lambda = list(time = 0, unit = 0, nn = Inf),
                 control = trop_control(n_boot = 20L))
#>        method estimate std.error conf.low conf.high n_treated_cells
#>  Bootstrap SE   2.5252    0.2263   2.1471    2.9478              12
#>  Jackknife SE   2.5252    0.2705   1.9951    3.0553              12
#>  n_treated_units outcome    engine rank note
#>                3       y cfcompare    0 <NA>
#>                3       y cfcompare    0 <NA>
# }
```
