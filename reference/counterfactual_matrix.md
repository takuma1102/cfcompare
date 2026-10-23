# Estimated counterfactual (untreated potential outcome) matrix

Returns a method's estimate of the untreated potential outcome \\Y(0)\\
as an `N x T` matrix on the panel's unit/time grid. This gives every
estimator a common interface, so counterfactuals – and the unit-period
prediction RMSE built on them
([`panel_rmse()`](https://takuma1102.github.io/cfcompare/reference/panel_rmse.md)
with `metric = "prediction"`) – can be compared on the same footing.

## Usage

``` r
counterfactual_matrix(object, ...)

predict_counterfactual(object, ...)
```

## Arguments

- object:

  A fitted object: a `trop` fit, a `cf_comparison` from
  [`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md),
  a synthdid estimate, or any object carrying a `counterfactual`
  element.

- ...:

  Unused; for method extensibility.

## Value

For a single fit, an `N x T` numeric matrix (with `NA` where the method
does not define a counterfactual). For a `cf_comparison`, a named list
of such matrices, one per method.

## Details

Coverage of the returned matrix differs by method. The native engines
([`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md),
and the DID/MC/TROP fits inside
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md))
return a fitted value for *every* cell. Methods whose model only defines
a counterfactual for the treated units – DIFP, and (when the optional
package is installed) SDID, SC and gsynth – fill the treated rows and
leave control cells `NA` (a control's untreated outcome is simply its
observed value).

`predict_counterfactual()` is an alias for `counterfactual_matrix()`.

## See also

[`panel_rmse()`](https://takuma1102.github.io/cfcompare/reference/panel_rmse.md),
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md),
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md)

## Examples

``` r
df <- sim_panel(N = 14, T = 9, n_treated = 3, t0 = 6, seed = 1)
fit <- trop(df, "y", "w", "id", "t",
            lambda = list(time = 0.1, unit = 0.5, nn = 2), se = "none")
M <- counterfactual_matrix(fit)
dim(M)
#> [1] 14  9
```
