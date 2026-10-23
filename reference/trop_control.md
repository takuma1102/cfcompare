# Control settings for [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)

Control settings for
[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)

## Usage

``` r
trop_control(
  max_iter = 2000L,
  tol = 1e-06,
  n_cv_cells = 120L,
  cv_cycles = 2L,
  max_cells = 60L,
  time_center = c("midpoint", "python"),
  pooled_treat_distance = NULL,
  conf_level = 0.95,
  n_boot = 200L,
  boot_resample = c("stratified", "controls"),
  boot_ci = c("percentile", "normal"),
  svd = c("truncated", "full"),
  cv_method = c("loocv", "placebo"),
  n_placebo = 30L,
  nn_floor = "auto",
  nn_floor_scale = 0.001,
  workers = 1L,
  seed = NULL
)
```

## Arguments

- max_iter:

  Maximum solver iterations. Higher by default because small `lambda_nn`
  converges slowly under the soft-impute iteration.

- tol:

  Reference solver convergence tolerance. The applied tolerance is
  tightened automatically as `lambda_nn` weakens (small penalties
  converge slowly), keeping accuracy roughly constant across the penalty
  grid; a fixed tolerance otherwise under-converges at small
  `lambda_nn`. See
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md).

- n_cv_cells:

  Number of control cells sampled for the CV criterion. Set to `Inf` (or
  `0`) to use *every* control cell, i.e. the paper's full eq. (5)
  criterion (the Stata command's `cells(0)`); slower but exact.

- cv_cycles:

  Number of coordinate-descent cycles in penalty selection.

- max_cells:

  Threshold for `anchor = "auto"` to switch to pooled weights.

- time_center:

  Where the pooled time weights `theta_s = exp(-lambda_time |s - c|)`
  are centred within a treated spell (`anchor = "pooled"`, the placebo
  CV criterion and the pooled counterfactual behind the plots / event
  study). `"midpoint"` (default) uses the true midpoint
  `c = (first + last) / 2` of the spell, the convention of the official
  Stata command (v0.2.4+); for a single treated period it reduces to the
  paper's per-cell distance `|s - t|`, so a single treated cell gets the
  same estimate under both anchors. `"python"` places the centre half a
  period later (`(first + last + 1) / 2`), the `T - tp/2` convention of
  the official Python package, which
  [`trop_matrix()`](https://takuma1102.github.io/cfcompare/reference/trop_matrix.md)
  always uses; choose it for like-for-like comparisons with that
  package.

- pooled_treat_distance:

  Optional overrides for the distances of a pooled group to *itself*, as
  a list with elements `unit` and/or `time` (single non-negative
  numbers; `NULL`, the default, keeps the usual distances). `unit`
  replaces the unit distance of every member of the group – which is
  otherwise its distance to the group average, hence generally non-zero
  with several members – by the given value, so `list(unit = 0)` gives
  the treated units the full unit weight; `time` replaces the time
  distance of the periods inside the spell by the given value, so
  `list(time = 0)` gives the treated periods the full time weight. This
  is the `pooled_treat_distance(unit(#) time(#))` option of the official
  Stata command (v0.2.6+); it applies to the pooled anchor, the
  placebo-CV criterion and the pooled counterfactual behind the plots /
  event study, and is ignored by the per-cell anchor.

- conf_level:

  Confidence level for intervals.

- n_boot:

  Number of replications for the bootstrap standard error
  (`se = "bootstrap"`).

- boot_resample:

  Which units the bootstrap resamples: `"stratified"` (default) draws
  the treated and the control units separately, each with replacement
  (the stratified block bootstrap of the paper and of the Stata
  command); `"controls"` keeps the treated units fixed and resamples
  only the control units (the default `resample = "controls"` of the
  official Python package's `trop.bootstrap`, v0.1.11+, suited to
  designs with a single or very few treated units).

- boot_ci:

  Bootstrap confidence-interval type: `"percentile"` or `"normal"`.

- svd:

  Singular-value decomposition used by the soft-impute solver:
  `"truncated"` (default) computes only the leading singular triplets
  with `RSpectra` when it is installed and the matrix is large enough to
  benefit, falling back to the full SVD otherwise; `"full"` always uses
  the exact base R [`svd()`](https://rdrr.io/r/base/svd.html). The two
  agree to numerical tolerance; `"truncated"` is faster on large panels,
  `"full"` is used for exact numerical-agreement checks.

- cv_method:

  Penalty cross-validation criterion. `"loocv"` (default) scores
  held-out control-cell prediction error (the paper's eq. (4)-(5));
  `"placebo"` instead stamps the *actual* treated-unit adoption patterns
  onto randomly resampled control units and minimises the mean squared
  placebo ATT (the placebo-RMSE / pattern-resampling criterion of the
  official Python and Stata packages). Under block adoption the two
  placebo variants coincide; under staggered adoption each placebo draw
  reproduces the real cohort structure, so `"placebo"` is the
  recommended choice there (a warning is emitted if `"loocv"` is used on
  a staggered design). Ignored when the penalties are supplied via
  `lambda` (CV is then skipped). The placebo criterion always uses the
  pooled (per-cohort block-centre, cell-weighted) ATT, even when the
  final fit uses `anchor = "per_cell"`.

- n_placebo:

  Number of placebo assignments drawn when `cv_method = "placebo"`
  (larger is more stable but slower).

- nn_floor:

  Small non-negative stabilising constant *added* to every finite
  `lambda_nn` inside the solver (`lambda_nn = Inf`, the DID/TWFE special
  case, is never affected). Rationale: the nuclear-norm proximal step
  soft-thresholds singular values by `lambda_nn / (2 max w)`, so at
  `lambda_nn = 0` the step degenerates to the identity, the completion
  of the treated cells is no longer pinned down by the program, and
  soft-impute converges slowly; a strictly positive floor keeps the
  weak-regularisation limit well posed. The default `"auto"` resolves to
  `nn_floor_scale * s1`, where `s1` is the operator norm (largest
  singular value) of the two-way-demeaned control outcome – the same
  data-scale reference the default `lambda_nn` grid is built from. The
  floor is therefore *relative to the outcome scale*: it is resolved on
  the fitting scale, so it adapts automatically when
  `standardize = TRUE` (both modes get the same floor in relative
  terms). Supply a single non-negative number to use that exact value
  instead (interpreted on the fitting scale, like `lambda`), or `0` to
  disable the adjustment entirely, in which case penalties are applied
  exactly as given – use `nn_floor = 0` for like-for-like numerical
  comparisons against other implementations at identical penalty values.
  The resolved value is stored in the fit as `$lambda$nn_floor` and
  printed when positive.

- nn_floor_scale:

  Relative constant used when `nn_floor = "auto"`: the resolved floor is
  `nn_floor_scale * s1` (default `1e-3`, i.e. 0.1\\ data's operator norm
  – an epsilon next to the default `lambda_nn` grid, whose smallest
  finite point is `0.04 * s1`). Ignored when `nn_floor` is numeric.

- workers:

  Number of parallel workers for the embarrassingly parallel loops
  (cross-validation cells, and the bootstrap / jackknife replicates).
  `1` (default) runs serially. Values `> 1` use `future.apply`/`future`
  when installed (a transient `multisession` plan is set up and restored
  automatically); if those packages are missing it falls back to serial
  with no error. Results are reproducible given `seed`.

- seed:

  Optional seed for CV cell sampling (reproducibility).

## Value

A list of control parameters.
