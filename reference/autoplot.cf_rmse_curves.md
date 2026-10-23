# Autoplot for paired RMSE curves

Autoplot for paired RMSE curves

## Usage

``` r
# S3 method for class 'cf_rmse_curves'
autoplot(object, combined = FALSE, log_y = TRUE, ...)
```

## Arguments

- object:

  A `cf_rmse_curves`.

- combined:

  `FALSE` (default) returns a named list of two `ggplot`s; `TRUE`
  returns a single side-by-side `ggplot` (needs the `patchwork`
  package).

- log_y:

  Logical; log10 y axis.

- ...:

  Unused.

## Value

A list of two `ggplot`s, or one combined `ggplot`.
