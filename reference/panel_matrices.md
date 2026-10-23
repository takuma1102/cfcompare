# Reshape a long panel into TROP matrices

Converts a long panel `data.frame` into the `N x T` outcome (`Y`) and
treatment (`W`) matrices used by
[`trop_matrix()`](https://takuma1102.github.io/cfcompare/reference/trop_matrix.md),
and reports the treated rows and the width of the treated block. It is
the data-frame-to-matrix companion to
[`trop_matrix()`](https://takuma1102.github.io/cfcompare/reference/trop_matrix.md):
rather than building `Y`, `W`, `treated_units` and `treated_periods` by
hand, call `panel_matrices()` and pass the pieces straight through.

## Usage

``` r
panel_matrices(data, outcome, treatment, unit, time)
```

## Arguments

- data:

  A long `data.frame`/`tibble` with one row per unit-time cell.

- outcome, treatment, unit, time:

  Column names (strings). `treatment` must be a 0/1 indicator.

## Value

A list with components:

- `Y`:

  `N x T` outcome matrix (rows = units, columns = times).

- `W`:

  `N x T` 0/1 treatment matrix on the same grid.

- `units`,`times`:

  the sorted unit and time labels.

- `treated_units`:

  integer row indices of units with any treated cell.

- `treated_periods`:

  number of periods in which any treated unit is active, i.e. the width
  of the post block.

## See also

[`trop_matrix()`](https://takuma1102.github.io/cfcompare/reference/trop_matrix.md),
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)

## Examples

``` r
df <- sim_panel(N = 14, T = 9, n_treated = 3, t0 = 6, att = 2, seed = 1)
pm <- panel_matrices(df, "y", "w", "id", "t")
str(pm)
#> List of 6
#>  $ Y              : num [1:14, 1:9] 0.502 0.277 0.546 1.543 -0.964 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:14] "1" "2" "3" "4" ...
#>   .. ..$ : chr [1:9] "1" "2" "3" "4" ...
#>  $ W              : num [1:14, 1:9] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:14] "1" "2" "3" "4" ...
#>   .. ..$ : chr [1:9] "1" "2" "3" "4" ...
#>  $ units          : int [1:14] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ times          : int [1:9] 1 2 3 4 5 6 7 8 9
#>  $ treated_units  : Named int [1:3] 1 2 3
#>   ..- attr(*, "names")= chr [1:3] "1" "2" "3"
#>  $ treated_periods: int 4
# Feed straight into the matrix-in entry point (treated rows/periods inferred).
trop_matrix(pm$Y, pm$W, lambda_unit = 0.1, lambda_time = 0.1, lambda_nn = Inf)
#> [1] 1.419739
```
