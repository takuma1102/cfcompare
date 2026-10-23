# cfcompare 0.2.0

## Alignment with Stata `trop` v0.2.9 and Python `trop` v0.1.11 (October 2026)

Changes since the July 2026 development release on GitHub / R-universe,
following the official Stata command's (v0.2.9) and Python package's
(v0.1.11) summer releases; see the two comparison tables in the README (as of
July 2026 and as of October 2026).

* `trop(anchor = "pooled")` now pools over *treatment spell groups* -- treated
  cells with the same start and end period, pooled across units -- solving
  eq. (2) once per group with weights anchored to that group, and returns the
  per-group effects in `$tau_groups` (start/end, unit and cell counts,
  estimate; the ATT is their cell-weighted average). This is the `group(time)`
  estimand of the official Stata command (v0.2.6+). Simultaneous permanent
  adoption is a single group, so block-design results are unchanged; staggered
  adoption now gets one anchor per cohort and a unit that switches in and out
  of treatment one anchor per spell, instead of a single anchor spanning the
  whole treated set. The group unit distance is to the average path of the
  group's units over the periods in which none of them is treated, with each
  comparison unit's own treated cells masked (the Stata masking); the
  placebo-CV criterion, the event study, `counterfactual_matrix()` and
  `autoplot()` use the same grouping.
* Pooled time weights are now centred on the true midpoint of the spell,
  `(first + last) / 2` (the Stata command's convention since v0.2.4), so a
  single treated cell gets the same estimate under the pooled and per-cell
  anchors. The previous centre, half a period later (`T - tp/2` in 0-based
  indexing, the official Python package's convention), remains available as
  `trop_control(time_center = "python")` and is what `trop_matrix()` always
  uses, for like-for-like checks against that package.
* New `trop_control(pooled_treat_distance = list(unit =, time =))`: overrides
  the distances of a pooled group to itself -- the members' unit distance
  (otherwise their distance to the group average, hence generally non-zero)
  and the in-spell time distance -- e.g. `list(unit = 0, time = 0)` gives the
  group's own units and periods the full weight. The equivalent of the Stata
  command's `pooled_treat_distance(unit(#) time(#))` option (v0.2.6+); it
  applies to the pooled anchor, the placebo CV and the pooled counterfactual,
  is stamped into `$lambda` like the nuclear-norm floor so downstream
  re-solves reproduce it, and is reported by `print()`.
* New `trop_control(boot_resample = c("stratified", "controls"))`:
  `"controls"` keeps the treated units fixed and resamples only the control
  units (the default `resample = "controls"` of the Python package's
  `trop.bootstrap`, v0.1.11+); `"stratified"` (default, unchanged) resamples
  treated and control units separately. Applies to the ATT and event-study
  bootstraps; the fit records the scheme in `$boot_resample`. The resampling
  now indexes through `sample.int()`, so a single treated unit is resampled
  from itself rather than from `1:unit`.
* `se = "jackknife"` now deletes one unit at a time over *every* unit
  (treated and control), as the Stata command's `vce(jackknife)`, instead of
  over the treated units only; replicates that leave no treated unit (or
  fail) are skipped and the variance uses the number of finite replicates. It
  therefore no longer requires two treated units. The event-study jackknife
  follows the same scheme. (`se = "auto"` still picks the jackknife when there
  are at least two treated units and the bootstrap otherwise.)

## First development release on GitHub / R-universe (July 2026)

* Native R implementations of the Triply
  RObust Panel (TROP) estimator of Athey, Imbens, Qu & Viviano (2026).
* `panel_compare()` runs DID, SDID, SC, MC, TROP, gsynth and augsynth on one
  dataset and returns a common tidy ATT schema with shared plots.
* `trop_matrix()` provides a matrix-in form whose results can be numerically compared with the official Python `trop`
  matrix-in signature for direct cross-checking; verified to match the reference
  to numerical tolerance on the weighted-TWFE special case.* `autoplot()` on a `trop` fit now draws the TROP time weights
  `theta_s = exp(-lambda_time * |t - s|)` as a filled band along the bottom
  (synthdid-style), replacing the earlier bar strip, and labels it so the band
  reads as time weights rather than observation counts. The band is shown for
  uniform weights too (`lambda_time = 0`, flat) and can be hidden with
  `show_weights = FALSE`.
* `compare_se_modes()` compares the TROP bootstrap, jackknife and placebo
  standard errors on a single fit: the cross-validated penalties are chosen once
  and reused, so the point estimate is identical across rows and only the
  uncertainty differs (the same fixed-penalty pattern as `trop_ablation()`).
