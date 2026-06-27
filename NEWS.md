# cfcompare 0.1.0

* Initial version.
* Native R implementations of DID/TWFE, matrix completion, and the Triply
  RObust Panel (TROP) estimator of Athey, Imbens, Qu & Viviano (2025).
* `panel_compare()` runs DID, SDID, SC, MC, TROP, gsynth and augsynth on one
  dataset and returns a common tidy ATT schema with shared plots.
* `trop_matrix()` provides a matrix-in form whose results can be numerically compared with the official Python `trop`
  matrix-in signature for direct cross-checking; verified to match the reference
  to numerical tolerance on the weighted-TWFE special case.
* Optional engines (`synthdid`, `gsynth`, `augsynth`) are skipped with a message
  when their package is unavailable or the design does not apply.

## Unreleased

* Performance: `trop_control()` gains `svd` and `workers`.
  - `svd = "truncated"` (new default) solves the soft-impute step with a
    truncated SVD (`RSpectra::svds`) when available and beneficial, falling back
    to the full base-R `svd()`; `svd = "full"` forces the exact full
    decomposition (used for the numerical-agreement checks). Replaces the old
    `options(cfcompare.truncated_svd=)` switch.
  - `workers > 1` runs the embarrassingly parallel loops -- cross-validation
    cells and the bootstrap / jackknife / placebo replicates, plus the
    `panel_rmse()` placebo runs -- in parallel via `future.apply`/`future`
    (serial fallback when absent). Resampling draws are generated up front so
    results stay reproducible given `seed` (exact bootstrap draws may differ from
    earlier serial versions).
  - Warm starts: in `anchor = "per_cell"`, each treated cell's soft-impute solve
    is initialised from the previous cell's low-rank fit.
  - The package still runs on base R alone; the optional paths activate only when
    `RSpectra` / `future.apply` are installed.
* `panel_rmse()` and `autoplot()` for it: cross-model out-of-sample RMSE
  comparison (the "random blocks" placebo validation from the doubly/triply
  robust panel estimator paper), visualised as a ranked bar chart.
* `se = "bootstrap"`: unit-level stratified block bootstrap for TROP standard
  errors and confidence intervals (the paper's inference approach), with
  `n_boot` / `boot_ci` controls.
* `sim_semisynthetic()`: build a ground-truth benchmark from a real panel by
  imposing a known (possibly dynamic) effect on a placebo block.
* `autoplot()` / `plot()` for a single `trop` fit: synthetic-control-style
  trajectory plot with the post-treatment gap and time-weight ribbon.
* README now states clearly that this is an unofficial package and cites the
  published Journal of Applied Econometrics (2026) version; added inst/CITATION.
* Added the Callaway & Sant'Anna staggered DID estimator as method `"CS"` (via
  the `did` package), and added `did` to Suggests.
* `panel_rmse()` gains `metric = "placebo"` (default): a placebo-ATT RMSE that
  scores **every** method -- including `SDID`/`SC` (synthdid), `gsynth`,
  `augsynth` and `CS` -- on a common footing. `metric = "prediction"` keeps the
  per-cell held-out RMSE for native methods.
* `.engine_synthdid()` now skips the expensive jackknife variance when no SE is
  requested (`se = "none"`), making large comparisons and placebo loops fast.
* README now shows an RMSE comparison figure under unobserved confounding.

# cfcompare 0.1.0 (continued)

* Added **DIFP** (Doudchenko–Imbens–Ferman–Pinto): a native demeaned synthetic
  control (SC after recentering the mean, with an intercept), the DIFP comparison
  estimator from the paper. Available in `panel_compare()` and `panel_rmse()` as
  method `"DIFP"`.
* Added `trop_sensitivity()` plus `autoplot()`/`plot()` methods: a penalty
  (lambda_time x lambda_nn) heatmap coloured by cross-validation loss and
  annotated with the ATT, with the CV-selected cell highlighted.
* Documentation reworded throughout to make clear that cfcompare is an
  independent, unofficial implementation written from the paper — not a port of,
  and sharing no code with, the official Python or Stata software.
* `CS` (Callaway & Sant'Anna) is documented as an opt-in method, outside the
  paper's comparison set and requiring an absorbing design.
