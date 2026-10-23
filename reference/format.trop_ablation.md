# Render a TROP ablation table as paste-ready LaTeX or Markdown

Produces a publication-quality table from a
[`trop_ablation()`](https://takuma1102.github.io/cfcompare/reference/trop_ablation.md)
result. The `"latex"` output uses booktabs rules and math-mode penalty
headers, matching the style of the paper's tables (load the booktabs
LaTeX package).

## Usage

``` r
# S3 method for class 'trop_ablation'
format(
  x,
  output = c("latex", "markdown"),
  digits = 3,
  caption = "TROP penalty ablation",
  label = "tab:trop-ablation",
  ...
)
```

## Arguments

- x:

  A `trop_ablation` object from
  [`trop_ablation()`](https://takuma1102.github.io/cfcompare/reference/trop_ablation.md).

- output:

  `"latex"` (booktabs) or `"markdown"` (GitHub pipe table).

- digits:

  Number of decimal places.

- caption, label:

  Table caption and cross-reference label (LaTeX only).

- ...:

  Unused.

## Value

A character vector of lines (paste-ready). Combine with
`cat(..., sep = "\n")` or
[`writeLines()`](https://rdrr.io/r/base/writeLines.html).

## See also

[`trop_ablation()`](https://takuma1102.github.io/cfcompare/reference/trop_ablation.md)

## Examples

``` r
df <- sim_panel(N = 14, T = 8, n_treated = 3, t0 = 5, att = 2, seed = 1)
ab <- trop_ablation(df, "y", "w", "id", "t",
                    lambda_full = list(time = 0.1, unit = 0.5, nn = 2))
writeLines(format(ab, "markdown"))
#> | Specification | &lambda;_time | &lambda;_unit | &lambda;_nn | ATT | rank(L) |
#> | :-- | --: | --: | --: | --: | --: |
#> | TROP (full) | 0.100 | 0.500 | 2.000 | 1.766 | 3 |
#> | No regression adjustment | 0.100 | 0.500 | Inf | 1.702 | 0 |
#> | No unit weights | 0.100 | 0.000 | 2.000 | 1.833 | 4 |
#> | No time weights | 0.000 | 0.500 | 2.000 | 1.808 | 4 |
#> | Matrix completion | 0.000 | 0.000 | 2.000 | 1.867 | 5 |
#> | Difference-in-differences | 0.000 | 0.000 | Inf | 1.815 | 0 |
```
