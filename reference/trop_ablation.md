# Penalty ablation for the TROP estimator

Refits TROP under a sequence of penalty constraints to show how the ATT
moves as the estimator is stripped back towards matrix completion and
difference-in-differences. This is the robustness exercise behind Table
5 of Athey, Imbens, Qu & Viviano (2026): constraining `lambda_unit`
and/or `lambda_time` to zero (dropping unit/time weights) and/or
`lambda_nn` to infinity (dropping the low-rank regression adjustment).

## Usage

``` r
trop_ablation(
  data,
  outcome,
  treatment,
  unit,
  time,
  covariates = NULL,
  control = trop_control(),
  anchor = "pooled",
  se = "none",
  lambda_full = NULL,
  rmse = FALSE,
  horizon = NULL,
  n_pseudo = NULL,
  n_runs = 200L
)

# S3 method for class 'trop_ablation'
print(x, digits = 3, ...)

# S3 method for class 'trop_ablation'
plot(x, file = NULL, width = NULL, height = NULL, res = 200, digits = 3, ...)
```

## Arguments

- data, outcome, treatment, unit, time, covariates:

  Passed to
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md).

- control:

  A
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md)
  list (used for every fit).

- anchor:

  Weight anchoring, as in
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)
  (default `"pooled"`).

- se:

  Standard-error method for every specification (default `"none"`; any
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)
  method is allowed, e.g. `"jackknife"`). When not `"none"`, the
  standard error and confidence interval appear in the table.

- lambda_full:

  Optional named list `list(time=, unit=, nn=)` to fix the full
  specification's penalties and skip cross-validation.

- rmse:

  Logical; if `TRUE`, also score each specification by a placebo RMSE. A
  pseudo block-treatment (true effect zero) is assigned to control
  units, every specification is refit at its penalties, and the RMSE is
  `sqrt(mean(estimate^2))` over `n_runs` draws. The same draws are
  scored under every specification, so the relative RMSE with TROP
  (full) = 1 is the quantity reported in Table 5 of the paper. Adds
  `rmse` and `rmse_rel` columns, and makes
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) show the
  (relative) RMSE instead of the ATT.

- horizon, n_pseudo:

  Placebo design used when `rmse = TRUE`: the number of final periods to
  pseudo-treat and the number of pseudo-treated control units. Default
  to the real design's post length and treated-unit count.

- n_runs:

  Number of placebo draws when `rmse = TRUE`.

- x:

  A `trop_ablation` object.

- digits:

  Number of decimal places for the penalties and estimates.

- ...:

  Unused.

- file:

  Optional output path. When supplied, the table is rendered to that
  file instead of the active graphics device; the format is taken from
  the extension (`.png`, the default, or `.pdf` for vector output). For
  LaTeX or Markdown source use
  [`format()`](https://rdrr.io/r/base/format.html) instead.

- width, height:

  Figure size in inches. Defaults adapt to the number of rows and
  columns.

- res:

  Resolution in PPI for the `.png` device (ignored for `.pdf`).

## Value

An object of class `trop_ablation` (a `data.frame`) with one row per
specification: `spec`, the three penalties `lambda_time`/`lambda_unit`/
`lambda_nn`, `estimate`, `std.error`, `conf.low`/`conf.high`, and
`rank`. When `rmse = TRUE` it also has `rmse` and `rmse_rel` (RMSE
relative to the full TROP).

## Details

The full specification's penalties are taken from a single
cross-validated
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)
fit; every constrained specification then reuses those values where they
are not constrained, so the rows are directly comparable. The result is
a table object: print it for a clean console table, or
[`format()`](https://rdrr.io/r/base/format.html) it to paste-ready LaTeX
(booktabs) or Markdown. It is also an ordinary `data.frame`, so it can
be passed to
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html),
`gt::gt()`, etc.

## References

Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026). Triply Robust
Panel Estimators. *Journal of Applied Econometrics*.
[doi:10.1002/jae.70061](https://doi.org/10.1002/jae.70061) .

## See also

[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md);
[`format()`](https://rdrr.io/r/base/format.html) for paste-ready
LaTeX/Markdown output

## Examples

``` r
df <- sim_panel(N = 14, T = 8, n_treated = 3, t0 = 5, att = 2, seed = 1)
# fixed base penalties keep the example fast; omit `lambda_full` to select
# them by cross-validation
ab <- trop_ablation(df, "y", "w", "id", "t",
                    lambda_full = list(time = 0.1, unit = 0.5, nn = 2))
ab                                   # clean console table
#> TROP penalty ablation  (outcome: y; anchor: pooled; SE: none)
#> 
#>  Specification             lambda_time lambda_unit lambda_nn ATT   rank(L)
#>  TROP (full)               0.100       0.500       2.000     1.766 3      
#>  No regression adjustment  0.100       0.500       Inf       1.702 0      
#>  No unit weights           0.100       0.000       2.000     1.833 4      
#>  No time weights           0.000       0.500       2.000     1.808 4      
#>  Matrix completion         0.000       0.000       2.000     1.867 5      
#>  Difference-in-differences 0.000       0.000       Inf       1.815 0      
# \donttest{
cat(format(ab, "latex"), sep = "\n") # paste-ready LaTeX (booktabs)
#> \begin{table}[!ht]
#> \centering
#> \caption{TROP penalty ablation}
#> \label{tab:trop-ablation}
#> \begin{tabular}{lrrrrr}
#> \toprule
#> Specification & $\lambda_{\mathrm{time}}$ & $\lambda_{\mathrm{unit}}$ & $\lambda_{\mathrm{nn}}$ & ATT & rank($L$) \\
#> \midrule
#> TROP (full) & 0.100 & 0.500 & 2.000 & 1.766 & 3 \\
#> No regression adjustment & 0.100 & 0.500 & $\infty$ & 1.702 & 0 \\
#> No unit weights & 0.100 & 0.000 & 2.000 & 1.833 & 4 \\
#> No time weights & 0.000 & 0.500 & 2.000 & 1.808 & 4 \\
#> Matrix completion & 0.000 & 0.000 & 2.000 & 1.867 & 5 \\
#> Difference-in-differences & 0.000 & 0.000 & $\infty$ & 1.815 & 0 \\
#> \bottomrule
#> \end{tabular}
#> \end{table}
# }
```
