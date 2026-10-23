# Line plot of estimation RMSE versus a design dimension

One line per estimator; the y axis is on a log scale (as in the paper).
Lower is better.

## Usage

``` r
# S3 method for class 'cf_rmse_curve'
autoplot(object, log_y = TRUE, ...)
```

## Arguments

- object:

  A `cf_rmse_curve` from
  [`rmse_curve()`](https://takuma1102.github.io/cfcompare/reference/rmse_curve.md).

- log_y:

  Logical; log10 y axis (default `TRUE`).

- ...:

  Unused.

## Value

A `ggplot` object.
