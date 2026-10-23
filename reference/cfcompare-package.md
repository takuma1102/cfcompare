# cfcompare: Compare DID, SDID, MC and the Triply Robust Panel estimator

`cfcompare` lets applied researchers run several panel causal-inference
estimators on the **same** data and compare them on one tidy schema and
one plot. It provides a native R implementation of the Triply RObust
Panel (TROP) estimator of Athey, Imbens, Qu & Viviano (2026), and a
comparison harness that also wraps difference-in-differences (TWFE),
matrix completion, and – via the synthdid package – synthetic
difference-in-differences and synthetic control.

## Details

Key functions:

- [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md) –
  the native TROP estimator.

- [`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md)
  – run DID / SDID / MC / TROP together.

- [`as_att()`](https://takuma1102.github.io/cfcompare/reference/as_att.md)
  – coerce external results into the shared schema.

- [`autoplot.cf_comparison()`](https://takuma1102.github.io/cfcompare/reference/autoplot.cf_comparison.md),
  [`plot_counterfactual()`](https://takuma1102.github.io/cfcompare/reference/plot_counterfactual.md)
  – visuals.

- [`sim_panel()`](https://takuma1102.github.io/cfcompare/reference/sim_panel.md)
  – simulate factor-model panels for experiments.

## References

Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026). Triply Robust
Panel Estimators. *Journal of Applied Econometrics*.
[doi:10.1002/jae.70061](https://doi.org/10.1002/jae.70061) .

## See also

Useful links:

- <https://github.com/takuma1102/cfcompare>

- <https://takuma1102.github.io/cfcompare/>

- Report bugs at <https://github.com/takuma1102/cfcompare/issues>

## Author

**Maintainer**: Takuma Iwasaki <iwasakit@stanford.edu> \[copyright
holder\]

Authors:

- Takuma Iwasaki <iwasakit@stanford.edu> \[copyright holder\]
