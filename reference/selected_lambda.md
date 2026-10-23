# CV-selected penalties from a TROP penalty-sensitivity grid

Returns the cross-validation-selected penalties of a
[`trop_sensitivity()`](https://takuma1102.github.io/cfcompare/reference/trop_sensitivity.md)
grid — the grid cell minimising the CV loss — as a
`list(time=, unit=, nn=)`, the same shape
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)
accepts for its `lambda` argument. So the data-driven choice on a grid
can be fed straight back into a fit, e.g.
`trop(..., lambda = selected_lambda(g))`, without digging the values out
of `attr(g, "selected")` by hand.

## Usage

``` r
selected_lambda(grid)
```

## Arguments

- grid:

  A `cf_trop_grid` from
  [`trop_sensitivity()`](https://takuma1102.github.io/cfcompare/reference/trop_sensitivity.md).

## Value

A named list with the selected `time`, `unit` and `nn` penalties. The
full selected row (including the ATT and CV loss there) remains
available in `attr(grid, "selected")`.

## See also

[`trop_sensitivity()`](https://takuma1102.github.io/cfcompare/reference/trop_sensitivity.md),
[`surface_matrix()`](https://takuma1102.github.io/cfcompare/reference/surface_matrix.md),
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)

## Examples

``` r
# \donttest{
df <- sim_panel(N = 14, T = 9, n_treated = 3, t0 = 6, att = 2, seed = 1)
g <- trop_sensitivity(df, "y", "w", "id", "t",
                      lambda_time = c(0, 0.25, 1),
                      control = trop_control(n_cv_cells = 5L, cv_cycles = 1L,
                                             tol = 1e-3, seed = 1L))
selected_lambda(g)
#> $time
#> [1] 1
#> 
#> $unit
#> [1] 0
#> 
#> $nn
#> [1] 4.658485
#> 
# }
```
