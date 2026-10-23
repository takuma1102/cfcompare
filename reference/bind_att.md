# Stack ATT results into one comparison table

Row-binds several ATT results —
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)
fits,
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md)
results, synthdid estimates, existing `cf_att_tbl`s, or anything
[`as_att()`](https://takuma1102.github.io/cfcompare/reference/as_att.md)
understands — into a single `cf_att_tbl` ready for
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html) or
[`plot_counterfactual()`](https://takuma1102.github.io/cfcompare/reference/plot_counterfactual.md).
Each argument is passed through
[`as_att()`](https://takuma1102.github.io/cfcompare/reference/as_att.md)
first, so multi-row inputs (e.g. a
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md)
result) keep all their rows, and columns are aligned automatically. This
replaces the manual `rbind(as.data.frame(...), ...)` +
`class(...) <- "cf_att_tbl"` dance.

## Usage

``` r
bind_att(...)
```

## Arguments

- ...:

  Objects coercible by
  [`as_att()`](https://takuma1102.github.io/cfcompare/reference/as_att.md),
  optionally named to set the `method` label.

## Value

A `cf_att_tbl` (a `data.frame`) with the inputs stacked row-wise.

## Details

Names given to the arguments become `method` labels, overriding the
label already on a *single-row* result; a name on a multi-row input is
ignored with a warning (it cannot sensibly relabel several methods at
once). For per-object control of other fields such as `outcome`, coerce
that object with
[`as_att()`](https://takuma1102.github.io/cfcompare/reference/as_att.md)
first and pass the result in unnamed.

## See also

[`as_att()`](https://takuma1102.github.io/cfcompare/reference/as_att.md),
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md),
[`compare_se_modes()`](https://takuma1102.github.io/cfcompare/reference/compare_se_modes.md)

## Examples

``` r
df <- sim_panel(N = 14, T = 9, n_treated = 3, t0 = 6, att = 2, seed = 1)
lam <- list(time = 0.1, unit = 0.5, nn = 2)   # fixed penalties (fast example)
f_pooled   <- trop(df, "y", "w", "id", "t", anchor = "pooled",
                   lambda = lam, se = "none")
f_per_cell <- trop(df, "y", "w", "id", "t", anchor = "per_cell",
                   lambda = lam, se = "none")
bind_att(pooled = f_pooled, per_cell = f_per_cell)
#>    method estimate std.error conf.low conf.high n_treated_cells n_treated_units
#>    pooled   1.6018        NA       NA        NA              12               3
#>  per_cell   1.5904        NA       NA        NA              12               3
#>  outcome    engine rank note
#>        y cfcompare    3 <NA>
#>        y cfcompare    2 <NA>
```
