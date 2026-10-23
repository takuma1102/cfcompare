# Plot observed vs predicted counterfactual trajectories

For each native engine in a comparison (DID, MC, TROP), draws the
average outcome over the treated units against the predicted control
counterfactual, to show how the methods extrapolate through the
post-treatment period.

## Usage

``` r
plot_counterfactual(x, methods = NULL)
```

## Arguments

- x:

  A `cf_comparison` from
  [`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md).

- methods:

  Optional subset of methods to draw.

## Value

A ggplot2 object.
