# Simulate a panel from a low-rank factor model

Generates a long panel whose untreated potential outcomes follow an
interactive-fixed-effects (factor) model on top of two-way fixed
effects, in the spirit of the data-generating processes in Athey,
Imbens, Qu & Viviano (2026). A block treatment is applied to `n_treated`
units from period `t0` onward, with a constant additive effect `att`.

## Usage

``` r
sim_panel(
  N = 30,
  T = 20,
  n_treated = 5,
  t0 = NULL,
  rank = 3L,
  att = 1,
  noise = 0.5,
  seed = NULL,
  n_cov = 0L,
  phi = NULL,
  confounding = 0
)
```

## Arguments

- N:

  Number of units.

- T:

  Number of periods.

- n_treated:

  Number of treated units.

- t0:

  First treated period (block design).

- rank:

  Number of latent factors.

- att:

  True treatment effect added to treated cells.

- noise:

  Standard deviation of idiosyncratic noise.

- seed:

  Optional RNG seed.

- n_cov:

  Number of time-varying covariates to generate (default `0`, none).
  Generated covariates are cell-level standard-normal draws whose linear
  index is part of the untreated potential outcome `Y(0)`, matching the
  additive covariate model of
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)
  (`L = X.phi + R`); pass their names to `covariates` in
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md) /
  [`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md)
  to recover the coefficients. Defaults to `length(phi)` when `phi` is
  supplied.

- phi:

  Optional numeric coefficient vector for the covariates (length
  `n_cov`). When `NULL` and `n_cov > 0`, coefficients are drawn from a
  standard normal. The coefficients actually used are returned on the
  `"phi"` attribute of the result.

- confounding:

  Non-negative strength of selection-on-confounders in the treatment
  assignment. With `confounding = 0` (default) the first `n_treated`
  units are treated independently of the latent factors and covariates,
  so treatment is ignorable and DID/TWFE is unbiased. With
  `confounding > 0` treated units are selected on their latent factor
  loadings and, when covariates are present (`n_cov > 0`), on their
  post-period covariate index \\X'\phi\\ as well. The factor part biases
  estimators without a low-rank term (DID, SDID, SC); the covariate part
  biases estimators that ignore covariates (SDID, SC, uncontrolled DID).
  A fully robust estimator (e.g. TROP with `covariates=`) removes both.
  Larger values give stronger selection.

## Value

A long `data.frame` with columns `id`, `t`, `y`, `w`, and the noiseless
counterfactual `y0` (useful for evaluating estimators). When `n_cov > 0`
it additionally carries covariate columns `x1`, ..., `x{n_cov}`, and the
true coefficients are stored on `attr(., "phi")`.

## Examples

``` r
df <- sim_panel(N = 16, T = 9, n_treated = 3, t0 = 6, seed = 42)
head(df)
#>    id t          y w          y0
#> 1   1 1 1.99808336 0  1.99808336
#> 17  1 2 0.45625369 0  0.45625369
#> 33  1 3 0.08249937 0  0.08249937
#> 49  1 4 0.60979037 0  0.60979037
#> 65  1 5 0.39818588 0  0.39818588
#> 81  1 6 0.75972041 1 -0.24027959

# Two covariates with known coefficients, recoverable by trop(covariates=):
dc <- sim_panel(N = 16, T = 9, n_treated = 3, t0 = 6,
                n_cov = 2, phi = c(1.2, -0.7), seed = 1)
attr(dc, "phi")
#>   x1   x2 
#>  1.2 -0.7 
```
