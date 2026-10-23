# Estimation-RMSE curves over both design dimensions

Convenience wrapper that runs
[`rmse_curve()`](https://takuma1102.github.io/cfcompare/reference/rmse_curve.md)
for both the number of control units and the number of pre-treatment
periods and bundles the two curves, so they can be plotted individually
(default) or side by side (`combined = TRUE`), as in the paper's
two-panel figure.

## Usage

``` r
rmse_curves(values_control = NULL, values_pre = NULL, ...)
```

## Arguments

- values_control, values_pre:

  Optional grids for each sweep (see
  [`rmse_curve()`](https://takuma1102.github.io/cfcompare/reference/rmse_curve.md)
  defaults).

- ...:

  Passed to
  [`rmse_curve()`](https://takuma1102.github.io/cfcompare/reference/rmse_curve.md)
  (`n_runs`, `methods`, base design, `att`, `noise`, `control`, `seed`,
  `parallel`, `verbose`, ...). Set `parallel = TRUE` (plus a
  `future::plan`) to run both sweeps in parallel.

## Value

A `cf_rmse_curves` object: a list with `$n_control` and `$n_pre`, each a
`cf_rmse_curve`.

## See also

[`rmse_curve()`](https://takuma1102.github.io/cfcompare/reference/rmse_curve.md),
[`autoplot.cf_rmse_curves()`](https://takuma1102.github.io/cfcompare/reference/autoplot.cf_rmse_curves.md)

## Examples

``` r
# \donttest{
# defaults are a large simulation; use a small grid + few reps for a quick look
g <- rmse_curves(values_control = c(12, 16), values_pre = c(5, 8), n_post = 3,
                 n_treated = 2, n_runs = 1, methods = c("DID", "TROP"),
                 control = trop_control(n_cv_cells = 5L, cv_cycles = 1L,
                                        tol = 1e-3, seed = 1L))
plot(g)                  # two separate figures (default)


# combined = TRUE needs the optional 'patchwork' package:
if (requireNamespace("patchwork", quietly = TRUE))
  plot(g, combined = TRUE) # side-by-side, paper-style
# }
```
