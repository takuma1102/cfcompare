# Build a semi-synthetic panel from real data

Takes a real long panel, uses its outcomes (optionally smoothed through
a low-rank-plus-two-way-fixed-effects fit) as the untreated potential
outcomes `Y(0)`, and imposes a *known* treatment effect on a chosen
block of units and periods. Because the baseline is real data but the
effect is known, the result is a ground-truth benchmark that "closely
matches" the real setting – the style of semi-synthetic experiment used
to evaluate panel estimators in Athey, Imbens, Qu & Viviano (2026). Pair
it with
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md)
or
[`panel_rmse()`](https://takuma1102.github.io/cfcompare/reference/panel_rmse.md)
to score estimators against the truth.

## Usage

``` r
sim_semisynthetic(
  data,
  outcome,
  unit,
  time,
  n_treated,
  t0 = NULL,
  att = 1,
  effect = NULL,
  baseline = c("observed", "lowrank"),
  lambda_nn = NULL,
  noise = 0,
  seed = NULL
)
```

## Arguments

- data:

  A real long `data.frame`, one row per unit-time.

- outcome, unit, time:

  Column names (strings).

- n_treated:

  Number of units to assign to the placebo treated group (sampled at
  random from all units).

- t0:

  First treated period (block design). Defaults to about three quarters
  of the way through the panel.

- att:

  Constant additive treatment effect imposed on treated cells. Ignored
  if `effect` is supplied.

- effect:

  Optional per-treated-period effect: a single number, or a numeric
  vector of length `T - t0 + 1` giving a dynamic effect path.

- baseline:

  `"observed"` uses the real outcomes directly as `Y(0)`; `"lowrank"`
  replaces them with a low-rank + two-way-FE fit (optionally plus
  resampled residuals scaled by `noise`) for a smoother synthetic
  baseline.

- lambda_nn:

  Nuclear-norm penalty for the `"lowrank"` baseline fit.

- noise:

  For `"lowrank"`, standard-deviation multiplier on resampled residuals
  added back to the fit (0 = noiseless baseline).

- seed:

  Optional RNG seed.

## Value

A long `data.frame` with columns `id`, `t`, `y`, `w`, `y0` (the imposed
untreated potential outcome) and `tau` (the true effect, 0 off
treatment).

## See also

[`sim_panel()`](https://takuma1102.github.io/cfcompare/reference/sim_panel.md),
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md),
[`panel_rmse()`](https://takuma1102.github.io/cfcompare/reference/panel_rmse.md)

## Examples

``` r
real <- sim_panel(N = 16, T = 10, n_treated = 0L, att = 0, seed = 1)
ss <- sim_semisynthetic(real, "y", "id", "t",
                        n_treated = 6, t0 = 8, att = 3, seed = 2)
mean(ss$tau[ss$w == 1])   # true ATT = 3
#> [1] 3
```
