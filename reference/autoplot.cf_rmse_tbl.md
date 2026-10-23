# Bar chart of out-of-sample RMSE across methods

Visualises a
[`panel_rmse()`](https://takuma1102.github.io/cfcompare/reference/panel_rmse.md)
result as a ranked bar chart (lowest RMSE first), with +/- 1
standard-error whiskers across placebo runs. This is the cross-model
RMSE comparison from the paper.

## Usage

``` r
# S3 method for class 'cf_rmse_tbl'
autoplot(object, ...)
```

## Arguments

- object:

  A `cf_rmse_tbl` from
  [`panel_rmse()`](https://takuma1102.github.io/cfcompare/reference/panel_rmse.md).

- ...:

  Unused.

## Value

A `ggplot` object.
