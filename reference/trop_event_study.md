# Per-period (event-study) effects from a TROP fit

Decomposes a fitted
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md) ATT
into per-period treatment effects indexed by event time (periods
relative to each unit's first treated period) and attaches pointwise
standard errors and confidence intervals. The same resampling used for
the overall SE is reused: each resample's per-cell effects are
aggregated by event time rather than into a single mean, so no extra
modelling assumptions are introduced. With `pre_periods = TRUE` the
pre-treatment gaps are returned as placebo / pre-trend points (a flat,
near-zero pre-period profile supports the design).

## Usage

``` r
trop_event_study(
  object,
  se = c("bootstrap", "jackknife", "none"),
  pre_periods = TRUE,
  control = NULL,
  ...
)
```

## Arguments

- object:

  A `trop` fit from
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md).

- se:

  Standard-error method for the per-period effects: `"bootstrap"`
  (default; unit-level block bootstrap, stratified unless
  `trop_control(boot_resample = "controls")`), `"jackknife"`
  (delete-one-unit over every unit, as for the overall ATT), or
  `"none"`.

- pre_periods:

  Logical; include pre-treatment event times as placebo / pre-trend
  points (default `TRUE`).

- control:

  A list of solver/CV/bootstrap settings from
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md);
  defaults to
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md).
  Use it to set `n_boot`, `boot_ci`, `seed`, and `workers` for the
  resampling.

- ...:

  Unused.

## Value

An object of class `trop_event_study`: a list whose `estimates` is a
`data.frame` with one row per event time (`event_time`, `estimate`,
`std.error`, `conf.low`, `conf.high`, `n_cells`, `period`), plus the
overall `att`, the `se.method`, and metadata.

## Details

Because the treated composition varies across bootstrap resamples and
each event time draws on relatively few cells, the per-period intervals
are wider than the overall ATT interval and cover each point
individually (they are not simultaneous bands).

## References

Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026). Triply Robust
Panel Estimators. *Journal of Applied Econometrics*.
[doi:10.1002/jae.70061](https://doi.org/10.1002/jae.70061) .

## See also

[`autoplot.trop_event_study()`](https://takuma1102.github.io/cfcompare/reference/autoplot.trop_event_study.md)
to plot the result.

## Examples

``` r
# \donttest{
df  <- sim_panel(N = 16, T = 9, n_treated = 4, t0 = 6, att = 3, seed = 1)
# fixed penalties keep the example fast; omit `lambda` to cross-validate
fit <- trop(df, "y", "w", "id", "t", se = "none",
            lambda = list(time = 0.1, unit = 0.5, nn = 2))
es  <- trop_event_study(fit, se = "bootstrap", pre_periods = TRUE,
                        control = trop_control(n_boot = 20L, seed = 1))
es
#> TROP event study (per-period ATT)
#>   overall ATT = 3.2006
#>   event times: 4 post-treatment, 5 pre-treatment (placebo)
#>   SE: bootstrap (20 reps); 95% pointwise CI
#>  event_time estimate std.error conf.low conf.high n_cells period
#>          -5 -0.38700    0.2014 -0.63272   0.05746       4    pre
#>          -4  0.23676    0.1404 -0.06144   0.39595       4    pre
#>          -3  0.10638    0.1967 -0.20472   0.44004       4    pre
#>          -2  0.04075    0.1970 -0.27970   0.37886       4    pre
#>          -1 -0.03997    0.2463 -0.46544   0.33532       4    pre
#>           0  3.19224    0.3285  2.73900   3.85245       4   post
#>           1  3.24039    0.4329  2.50574   3.99350       4   post
#>           2  3.74927    0.5909  2.94867   4.79362       4   post
#>           3  2.94782    0.6199  2.00870   3.98453       4   post
# }
```
