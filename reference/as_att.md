# Coerce estimator output to the cfcompare ATT schema

Brings the output of
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md),
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md),
a synthdid estimate, or any object with an `estimate` into the common
`cf_att_tbl` schema, so results computed elsewhere can be slotted into
the same comparison and plots.

## Usage

``` r
as_att(x, ...)
```

## Arguments

- x:

  An object to tidy.

- ...:

  Passed to methods (e.g. `method`, `outcome`, `conf_level`).

## Value

A `cf_att_tbl` (a `data.frame`).
