# cfcompare 0.1.0

* Initial version.
* Native R implementations of the Triply
  RObust Panel (TROP) estimator of Athey, Imbens, Qu & Viviano (2026).
* `panel_compare()` runs DID, SDID, SC, MC, TROP, gsynth and augsynth on one
  dataset and returns a common tidy ATT schema with shared plots.
* `trop_matrix()` provides a matrix-in form whose results can be numerically compared with the official Python `trop`
  matrix-in signature for direct cross-checking; verified to match the reference
  to numerical tolerance on the weighted-TWFE special case.