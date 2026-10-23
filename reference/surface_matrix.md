# Surface matrix from a TROP penalty-sensitivity grid

Reshapes a long
[`trop_sensitivity()`](https://takuma1102.github.io/cfcompare/reference/trop_sensitivity.md)
grid into a wide matrix of one quantity — the cross-validation loss or
the ATT — laid out over the two swept penalties (see the `axes` argument
of
[`trop_sensitivity()`](https://takuma1102.github.io/cfcompare/reference/trop_sensitivity.md)).
Rows are the y-axis penalty, columns the x-axis penalty; the matrix
`dimnames` are named after the two penalty columns so the layout is
self-describing. This is the same reshape
[`plot_trop_surfaces()`](https://takuma1102.github.io/cfcompare/reference/plot_trop_surfaces.md)
performs internally, exposed so the surfaces can be inspected or
post-processed without re-deriving them with
[`xtabs()`](https://rdrr.io/r/stats/xtabs.html).

## Usage

``` r
surface_matrix(grid, value = c("cv_loss", "att"))
```

## Arguments

- grid:

  A `cf_trop_grid` from
  [`trop_sensitivity()`](https://takuma1102.github.io/cfcompare/reference/trop_sensitivity.md).

- value:

  Which quantity to lay out: `"cv_loss"` (default) or `"att"`.

## Value

A numeric matrix (y-axis penalty in rows, x-axis penalty in columns).

## See also

[`trop_sensitivity()`](https://takuma1102.github.io/cfcompare/reference/trop_sensitivity.md),
[`selected_lambda()`](https://takuma1102.github.io/cfcompare/reference/selected_lambda.md),
[`plot_trop_surfaces()`](https://takuma1102.github.io/cfcompare/reference/plot_trop_surfaces.md)

## Examples

``` r
# \donttest{
df <- sim_panel(N = 14, T = 9, n_treated = 3, t0 = 6, att = 2, seed = 1)
g <- trop_sensitivity(df, "y", "w", "id", "t",
                      lambda_time = c(0, 0.25, 1),
                      control = trop_control(n_cv_cells = 5L, cv_cycles = 1L,
                                             tol = 1e-3, seed = 1L))
surface_matrix(g, "cv_loss")
#>            lambda_nn
#> lambda_time 0.186339404267937 0.465848510669843 0.931697021339686
#>        0            0.8943407         0.8943685         0.7885174
#>        0.25         0.8943148         0.8262141         0.6467792
#>        1            0.4546796         0.4106894         0.3766296
#>            lambda_nn
#> lambda_time 2.32924255334922 4.65848510669843
#>        0           0.5794293        0.4497187
#>        0.25        0.4355402        0.4332015
#>        1           0.3036015        0.2982285
surface_matrix(g, "att")
#>            lambda_nn
#> lambda_time 0.186339404267937 0.465848510669843 0.931697021339686
#>        0             1.660882          1.661209          1.661421
#>        0.25          1.660988          1.661288          1.658145
#>        1             1.886077          1.935183          1.926282
#>            lambda_nn
#> lambda_time 2.32924255334922 4.65848510669843
#>        0            1.659234         1.579937
#>        0.25         1.607711         1.522923
#>        1            1.891111         1.830940
# }
```
