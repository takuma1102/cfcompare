# Contributing to cfcompare

Thanks for taking the time to contribute!

## Filing issues

When filing a bug report, please include a minimal reproducible example
(a [reprex](https://reprex.tidyverse.org)) and the output of
`sessionInfo()`. For feature requests, describe the use case and, where it
helps, the estimator/paper you have in mind.

## Pull requests

* Discuss non-trivial changes in an issue first so we agree on the approach.
* Fork the repo and create a branch from `main`.
* Keep the package passing `R CMD check`:

  ```r
  devtools::check()
  ```

* Style: this package follows the tidyverse style guide and ships an
  `.lintr`. Run `lintr::lint_package()` and `styler::style_pkg()` before
  pushing.
* Documentation is generated with roxygen2 — edit the roxygen comments in
  `R/`, then run `devtools::document()`. Do not hand-edit files in `man/` or
  `NAMESPACE`.
* Add a test under `tests/testthat/` for any behaviour change, and a bullet
  to `NEWS.md`.
* New optional estimators belong in `Suggests` (or `Remotes`) and must
  degrade gracefully — skip with a clear message when the engine package is
  not installed, the way the existing engines do.

## Code of Conduct

By participating you agree to abide by the
[Code of Conduct](CODE_OF_CONDUCT.md).
