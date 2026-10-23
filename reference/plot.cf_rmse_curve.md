# Plot one or both estimation-RMSE curves

Plot one or both estimation-RMSE curves

## Usage

``` r
# S3 method for class 'cf_rmse_curve'
plot(x, log_y = TRUE, file = NULL, width = 8, height = 5, ...)

# S3 method for class 'cf_rmse_curves'
plot(
  x,
  combined = FALSE,
  log_y = TRUE,
  file = NULL,
  width = NULL,
  height = NULL,
  ...
)
```

## Arguments

- x:

  A `cf_rmse_curve` or `cf_rmse_curves`.

- log_y:

  Logical; log10 y axis (default `TRUE`).

- file:

  Optional path. If given, the figure is written to a PNG at `width` x
  `height` inches; for separate panels two files are written with
  `-control`/`-pre` suffixes. If `NULL`, draws to the current device.

- width, height:

  Figure size in inches. Defaults are deliberately generous (single 8x5;
  combined 12x5).

- ...:

  Unused.

- combined:

  For `cf_rmse_curves`: `FALSE` (default) draws the two panels
  separately; `TRUE` draws them side by side (paper-style).

## Value

The input, invisibly.
