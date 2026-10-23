# Heatmap of the TROP penalty-sensitivity grid

Cells are coloured by cross-validation loss (darker = better
out-of-sample fit) and annotated with the ATT estimate at that penalty
pair; the CV-selected cell is outlined. The two swept penalties (see the
`axes` argument of
[`trop_sensitivity()`](https://takuma1102.github.io/cfcompare/reference/trop_sensitivity.md))
form the x and y axes, and the fixed penalty is reported in the
subtitle.

## Usage

``` r
# S3 method for class 'cf_trop_grid'
autoplot(object, ...)
```

## Arguments

- object:

  A `cf_trop_grid` from
  [`trop_sensitivity()`](https://takuma1102.github.io/cfcompare/reference/trop_sensitivity.md).

- ...:

  Unused.

## Value

A `ggplot` object.
