# cfcompare 0.1.0.9006

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
