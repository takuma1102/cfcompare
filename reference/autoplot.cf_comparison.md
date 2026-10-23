# Forest plot of ATT estimates across methods

Point estimate and (where available) confidence interval for each
method, on a single axis – the quick visual comparison this package is
built for.

## Usage

``` r
# S3 method for class 'cf_comparison'
autoplot(object, ...)

# S3 method for class 'cf_att_tbl'
autoplot(object, ...)
```

## Arguments

- object:

  A `cf_comparison` (from
  [`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md))
  or a `cf_att_tbl`.

- ...:

  Unused.

## Value

A ggplot2 object.

## Examples

``` r
# \donttest{
df <- sim_panel(N = 14, T = 9, n_treated = 3, t0 = 6, seed = 1)
cmp <- panel_compare(df, "y", "w", "id", "t",
                     methods = c("DID", "MC", "TROP"), se = "none",
                     control = trop_control(n_cv_cells = 5L, cv_cycles = 1L,
                                            tol = 1e-3, seed = 1L))
autoplot(cmp)

# }
```
