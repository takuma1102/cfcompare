# TROP estimate on matrix input

A thin, matrix-in wrapper around the TROP working model, written
independently from the paper. The matrix-in form is convenient for
numerically comparing
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)
against other matrix-based implementations on identical inputs. For
data-frame input, cross-validated penalty selection, inference and the
multi-estimator comparison, use
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md) and
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md)
instead.

## Usage

``` r
trop_matrix(
  Y,
  W,
  treated_units = NULL,
  lambda_unit,
  lambda_time,
  lambda_nn,
  treated_periods = NULL,
  control = trop_control()
)
```

## Arguments

- Y:

  N x T outcome matrix.

- W:

  N x T 0/1 treatment matrix (1 = actively treated cell).

- treated_units:

  Integer row indices of the treated units, used to anchor the unit
  weights. Defaults to `NULL`, in which case the treated rows are
  inferred from `W` as the rows with at least one treated cell.

- lambda_unit, lambda_time:

  Non-negative decay parameters for the unit and time weights.

- lambda_nn:

  Nuclear-norm penalty; use `Inf` to drop the low-rank term. As in
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md),
  the solver adds the small scale-relative stabilising floor `nn_floor`
  (see
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md))
  to a finite `lambda_nn`; `Inf` is never affected, so the weighted-TWFE
  special case used for the exact numerical-agreement checks is
  unchanged. For like-for-like comparisons against other implementations
  at identical finite penalty values, pass
  `control = trop_control(nn_floor = 0)` so `lambda_nn` is applied
  exactly as supplied.

- treated_periods:

  Number of final columns treated as the post block, used to build the
  pre-period mask and the time-distance centre. Defaults to `NULL`, in
  which case it is inferred from `W` as the number of periods in which
  any treated unit is active.

- control:

  A list of solver controls from
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md).

## Value

A single numeric ATT estimate.

## Details

The untreated outcome model is fitted on the control cells with the
supplied weights, and the returned effect is the average over treated
cells of `Y - alpha - beta - L`, exactly as in the paper's eq. (2). With
`lambda_nn = Inf` the low-rank term is dropped and the fit is weighted
two-way fixed effects, matching the reference to numerical tolerance.
With a finite `lambda_nn` the nuclear-norm term is solved by this
package's proximal-gradient / soft-impute routine, which targets the
same convex program as eq. (2) (loss
`sum_js w (Y - alpha - beta - L)^2`, with no 1/2 factor); `lambda_nn` is
therefore on the same scale as the official Python and Stata packages,
and finite-penalty results agree with them to solver tolerance.

`trop_matrix()` always uses the pooled weighting scheme anchored to the
centre of the treated block: it takes no `anchor` argument and has no
per-cell mode. It therefore computes the same pooled, block-centre
estimator as the official Python package's `TROP_TWFE_average`
(implemented independently), and exists as a like-for-like matrix-in
cross-check against that pooled reference.
`trop(..., anchor = "pooled")` uses the same convention (unit distances
to the *average* treated trajectory, time distances to the centre of the
treated block), so on a trailing-block design the two agree to solver
tolerance for the same penalties. For the per-cell estimator (the
paper's eq. (12) / Algorithm 2), use `trop(..., anchor = "per_cell")`,
which shares the data-frame core `.trop_att()`.

## References

Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026). *Triply Robust
Panel Estimators.* Journal of Applied Econometrics.
[doi:10.1002/jae.70061](https://doi.org/10.1002/jae.70061) .

## See also

[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md),
[`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md)

## Examples

``` r
df <- sim_panel(N = 14, T = 9, n_treated = 3, t0 = 6, att = 2, seed = 1)
Y <- matrix(0, max(df$id), max(df$t)); W <- Y
for (k in seq_len(nrow(df))) { Y[df$id[k], df$t[k]] <- df$y[k]
  W[df$id[k], df$t[k]] <- df$w[k] }
tu <- which(rowSums(W) > 0)
trop_matrix(Y, W, tu, 0.1, 0.1, Inf, treated_periods = 4)
#> [1] 1.419739

# treated_units and treated_periods are inferred from W when omitted, so the
# matrices produced by panel_matrices() can be fed in directly.
pm <- panel_matrices(df, "y", "w", "id", "t")
trop_matrix(pm$Y, pm$W, lambda_unit = 0.1, lambda_time = 0.1, lambda_nn = Inf)
#> [1] 1.419739
```
