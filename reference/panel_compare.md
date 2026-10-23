# Compare DID, SDID, MC and TROP on the same panel

Runs a set of panel estimators on a single long panel and returns their
average-treatment-effect-on-the-treated (ATT) estimates on a common tidy
schema, so applied researchers can compare them at a glance. DID, MC and
TROP are computed natively (no external dependencies); SDID and SC are
routed through the synthdid package when available, and an alternative
MC/IFE can be routed through gsynth. Methods whose optional package is
missing, or that do not apply to the design, are skipped with a message.

## Usage

``` r
panel_compare(
  data,
  outcome,
  treatment,
  unit,
  time,
  covariates = NULL,
  methods = c("TROP", "DID", "SC", "MC", "SDID", "DIFP"),
  exclude = NULL,
  anchor = "auto",
  se = c("auto", "jackknife", "bootstrap", "none"),
  control = trop_control(),
  verbose = FALSE
)
```

## Arguments

- data:

  A long `data.frame` with one row per unit-time.

- outcome, treatment, unit, time:

  Column names (strings). `treatment` must be a 0/1 indicator of
  *active* treatment in that cell (works for block, staggered, and
  non-absorbing designs).

- covariates:

  Optional character vector of time-varying covariate column names.
  Passed to the native engines (`"DID"`, `"MC"`, `"TROP"`), which share
  the additive covariate model \\L\_{js} = X\_{js}'\phi + R\_{js}\\; the
  synthdid-based methods (`"SDID"`, `"SC"`), `"DIFP"` and the optional
  engines (`"gsynth"`, `"augsynth"`, `"CS"`) do not receive covariates
  here and are run without them. Covariates must be fully observed
  across all unit-time cells.

- methods:

  Character vector of methods to run. Any of `"DID"`, `"SDID"`, `"SC"`,
  `"MC"`, `"DIFP"`, `"TROP"`, `"gsynth"`, `"augsynth"`, `"CS"`. Defaults
  to the native engines plus SDID and SC:
  `c("TROP", "DID", "SC", "MC", "SDID", "DIFP")`.

- exclude:

  Optional character vector of methods to drop from `methods` (after
  defaults are applied). Convenient for running "everything except one",
  e.g. `exclude = "DIFP"`. Unknown names are ignored with a warning.

- anchor:

  Weight anchoring for TROP; see
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md).

- se:

  Standard-error method for the native engines; see
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md).

- control:

  A list of solver/CV settings from
  [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md).

- verbose:

  Logical; print CV progress.

## Value

An object of class `cf_comparison`: a list with the tidy table `att`
(class `cf_att_tbl`), per-method counterfactual matrices
`counterfactual`, the native fit objects `fits`, and the reshaped
`panel`.

## See also

[`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md),
[`autoplot.cf_comparison()`](https://takuma1102.github.io/cfcompare/reference/autoplot.cf_comparison.md),
[`plot_counterfactual()`](https://takuma1102.github.io/cfcompare/reference/plot_counterfactual.md)

## Examples

``` r
# \donttest{
# panel_compare() fits TROP with cross-validated penalties, so it is wrapped
# in \donttest; few CV cells and a loose solver tolerance keep it quick
# (use the trop_control() defaults for real analyses)
df <- sim_panel(N = 14, T = 9, n_treated = 3, t0 = 6, seed = 3)
cmp <- panel_compare(df, "y", "w", "id", "t",
                     methods = c("DID", "MC", "TROP"), se = "none",
                     control = trop_control(n_cv_cells = 5L, cv_cycles = 1L,
                                            tol = 1e-3, seed = 1L))
cmp$att
#>  method estimate std.error conf.low conf.high n_treated_cells n_treated_units
#>     DID   1.2562        NA       NA        NA              12               3
#>      MC   1.3339        NA       NA        NA              12               3
#>    TROP   1.1285        NA       NA        NA              12               3
#>  outcome    engine rank note
#>        y cfcompare    0 <NA>
#>        y cfcompare    8 <NA>
#>        y cfcompare    0 <NA>
# }
```
