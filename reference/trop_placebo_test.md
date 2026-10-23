# Placebo (randomisation) test for a TROP ATT

Assesses whether a fitted
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md) ATT
is distinguishable from zero by comparing it to a placebo null
distribution. Control units are repeatedly relabelled as pseudo-treated
— using the *same number of treated units and the same treatment timing*
as the real design, so block, staggered and non-absorbing patterns are
all reproduced faithfully — and the estimator is refit at the fit's
cross-validated penalties (the true effect is zero by construction). The
observed ATT is then placed against that null.

## Usage

``` r
trop_placebo_test(
  object,
  B = 500L,
  alternative = c("two.sided", "greater", "less"),
  control = NULL
)

# S3 method for class 'trop_placebo_test'
print(x, digits = 4, ...)
```

## Arguments

- object:

  A `trop` fit from
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md).

- B:

  Number of placebo assignments to draw.

- alternative:

  Direction of the test: `"two.sided"` (default), `"greater"` or
  `"less"`.

- control:

  A
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md)
  list for the placebo refits (penalties are held fixed at the fit's
  values; this mainly sets `seed` and `workers`). Defaults to
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md);
  set `seed` here for reproducible draws.

- x:

  A `trop_placebo_test` object.

- digits:

  Number of significant digits to print.

- ...:

  Unused.

## Value

An object of class `trop_placebo_test` (a list) with the `observed` ATT,
the vector of `placebo` ATTs, the `p.value`, the placebo-null `mean` and
`sd`, a central null interval `null.low`/`null.high`, `alternative`,
`B`, and `n_treated_units`. It has
[`print()`](https://rdrr.io/r/base/print.html) and
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
methods.

## Details

This is a randomisation/permutation test, not a standard error: it
answers "could a TROP design like mine produce an effect this large when
there is none?" For interval estimates use the bootstrap or jackknife
(see
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md) and
[`compare_se_modes()`](https://takuma1102.github.io/cfcompare/reference/compare_se_modes.md)).

## See also

[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md),
[`compare_se_modes()`](https://takuma1102.github.io/cfcompare/reference/compare_se_modes.md)

## Examples

``` r
# \donttest{
df  <- sim_panel(N = 16, T = 9, n_treated = 3, t0 = 6, att = 3, seed = 1)
# fixed penalties keep the example fast; omit `lambda` to cross-validate
fit <- trop(df, "y", "w", "id", "t", se = "none",
            lambda = list(time = 0.1, unit = 0.5, nn = 2))
pt  <- trop_placebo_test(fit, B = 30)
pt
#> TROP placebo test
#>   observed ATT : 3.464
#>   placebo null : mean -0.07527, sd 0.3276  (30 draws, 3 pseudo-treated units)
#>   null 95% range: [-0.6849, 0.5333]
#>   two.sided p-value : 0.0000
# }
```
