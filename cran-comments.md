## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

  Possibly misspelled words in DESCRIPTION are author surnames (Athey, Imbens,
  Qu, Viviano), package names ('synthdid', 'gsynth', 'augsynth') and the
  estimator's acronym (TROP, "Triply RObust Panel"); they are spelled
  correctly.

## Suggested packages not on CRAN

The suggested packages 'synthdid' and 'augsynth' are not on CRAN. They are
used only conditionally (`requireNamespace()`); every feature that needs them
is skipped with an informative message when they are absent, and the tests
skip accordingly. Where to obtain them is stated in the Description field
(GitHub: synth-inference/synthdid, ebenmichael/augsynth).

## Test environments

* local: (fill in, e.g. macOS 15, R 4.5.x)
* win-builder: R-devel (devtools::check_win_devel())
* R-hub: linux, windows, macos (rhub::rhub_check()) -- optional
