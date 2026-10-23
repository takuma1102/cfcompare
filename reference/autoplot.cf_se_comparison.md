# Plot a TROP standard-error comparison

Forest plot of a
[`compare_se_modes()`](https://takuma1102.github.io/cfcompare/reference/compare_se_modes.md)
result: every row carries the same TROP point estimate and only the
confidence interval differs between resampling schemes. Zero is always
kept in view and marked with a dashed reference line, since whether the
interval covers zero is the quantity of interest.

## Usage

``` r
# S3 method for class 'cf_se_comparison'
autoplot(object, ...)
```

## Arguments

- object:

  A `cf_se_comparison` from
  [`compare_se_modes()`](https://takuma1102.github.io/cfcompare/reference/compare_se_modes.md).

- ...:

  Unused.

## Value

A ggplot2 object.
