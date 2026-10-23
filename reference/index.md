# Package index

## Compare estimators

Run DID/TWFE, synthetic control, SDID, matrix completion, DIFP and TROP
(plus optional gsynth, augsynth and Callaway–Sant’Anna) on one panel,
and put the results on a single tidy ATT schema.

- [`panel_compare()`](https://takuma1102.github.io/cfcompare/reference/panel_compare.md)
  : Compare DID, SDID, MC and TROP on the same panel
- [`as_att()`](https://takuma1102.github.io/cfcompare/reference/as_att.md)
  : Coerce estimator output to the cfcompare ATT schema
- [`bind_att()`](https://takuma1102.github.io/cfcompare/reference/bind_att.md)
  : Stack ATT results into one comparison table
- [`autoplot(`*`<cf_comparison>`*`)`](https://takuma1102.github.io/cfcompare/reference/autoplot.cf_comparison.md)
  [`autoplot(`*`<cf_att_tbl>`*`)`](https://takuma1102.github.io/cfcompare/reference/autoplot.cf_comparison.md)
  : Forest plot of ATT estimates across methods
- [`counterfactual_matrix()`](https://takuma1102.github.io/cfcompare/reference/counterfactual_matrix.md)
  [`predict_counterfactual()`](https://takuma1102.github.io/cfcompare/reference/counterfactual_matrix.md)
  : Estimated counterfactual (untreated potential outcome) matrix
- [`plot_counterfactual()`](https://takuma1102.github.io/cfcompare/reference/plot_counterfactual.md)
  : Plot observed vs predicted counterfactual trajectories

## TROP estimator

The Triply RObust Panel estimator of Athey, Imbens, Qu and Viviano
(2026): data-frame and matrix entry points, and solver /
cross-validation / inference settings.

- [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md) :
  Triply RObust Panel (TROP) estimator

- [`trop_control()`](https://takuma1102.github.io/cfcompare/reference/trop_control.md)
  :

  Control settings for
  [`trop()`](https://takuma1102.github.io/cfcompare/reference/trop.md)

- [`trop_matrix()`](https://takuma1102.github.io/cfcompare/reference/trop_matrix.md)
  : TROP estimate on matrix input

- [`panel_matrices()`](https://takuma1102.github.io/cfcompare/reference/panel_matrices.md)
  : Reshape a long panel into TROP matrices

## TROP diagnostics

Inspect a TROP fit: trajectories, event-study effects, placebo tests,
standard-error comparisons, penalty ablations and penalty-sensitivity
surfaces.

- [`autoplot(`*`<trop>`*`)`](https://takuma1102.github.io/cfcompare/reference/autoplot.trop.md)
  : Synthetic-control-style trajectory plot for a single TROP fit
- [`trop_event_study()`](https://takuma1102.github.io/cfcompare/reference/trop_event_study.md)
  : Per-period (event-study) effects from a TROP fit
- [`autoplot(`*`<trop_event_study>`*`)`](https://takuma1102.github.io/cfcompare/reference/autoplot.trop_event_study.md)
  : Event-study plot for a TROP fit
- [`trop_placebo_test()`](https://takuma1102.github.io/cfcompare/reference/trop_placebo_test.md)
  [`print(`*`<trop_placebo_test>`*`)`](https://takuma1102.github.io/cfcompare/reference/trop_placebo_test.md)
  : Placebo (randomisation) test for a TROP ATT
- [`autoplot(`*`<trop_placebo_test>`*`)`](https://takuma1102.github.io/cfcompare/reference/autoplot.trop_placebo_test.md)
  : Plot a TROP placebo test
- [`compare_se_modes()`](https://takuma1102.github.io/cfcompare/reference/compare_se_modes.md)
  : Compare TROP standard-error methods on one fit
- [`autoplot(`*`<cf_se_comparison>`*`)`](https://takuma1102.github.io/cfcompare/reference/autoplot.cf_se_comparison.md)
  : Plot a TROP standard-error comparison
- [`trop_ablation()`](https://takuma1102.github.io/cfcompare/reference/trop_ablation.md)
  [`print(`*`<trop_ablation>`*`)`](https://takuma1102.github.io/cfcompare/reference/trop_ablation.md)
  [`plot(`*`<trop_ablation>`*`)`](https://takuma1102.github.io/cfcompare/reference/trop_ablation.md)
  : Penalty ablation for the TROP estimator
- [`format(`*`<trop_ablation>`*`)`](https://takuma1102.github.io/cfcompare/reference/format.trop_ablation.md)
  : Render a TROP ablation table as paste-ready LaTeX or Markdown
- [`trop_sensitivity()`](https://takuma1102.github.io/cfcompare/reference/trop_sensitivity.md)
  : TROP penalty-sensitivity grid (heatmap data)
- [`autoplot(`*`<cf_trop_grid>`*`)`](https://takuma1102.github.io/cfcompare/reference/autoplot.cf_trop_grid.md)
  : Heatmap of the TROP penalty-sensitivity grid
- [`surface_matrix()`](https://takuma1102.github.io/cfcompare/reference/surface_matrix.md)
  : Surface matrix from a TROP penalty-sensitivity grid
- [`selected_lambda()`](https://takuma1102.github.io/cfcompare/reference/selected_lambda.md)
  : CV-selected penalties from a TROP penalty-sensitivity grid
- [`plot_trop_surfaces()`](https://takuma1102.github.io/cfcompare/reference/plot_trop_surfaces.md)
  : Base-R surface views of a TROP penalty-sensitivity grid

## Held-out RMSE and simulation

Score estimators on held-out control data, trace estimation RMSE over a
design dimension, and simulate panels with a known effect.

- [`panel_rmse()`](https://takuma1102.github.io/cfcompare/reference/panel_rmse.md)
  : Out-of-sample RMSE across panel estimators
- [`autoplot(`*`<cf_rmse_tbl>`*`)`](https://takuma1102.github.io/cfcompare/reference/autoplot.cf_rmse_tbl.md)
  : Bar chart of out-of-sample RMSE across methods
- [`rmse_curve()`](https://takuma1102.github.io/cfcompare/reference/rmse_curve.md)
  : Estimation-RMSE curve over one design dimension
- [`rmse_curves()`](https://takuma1102.github.io/cfcompare/reference/rmse_curves.md)
  : Estimation-RMSE curves over both design dimensions
- [`autoplot(`*`<cf_rmse_curve>`*`)`](https://takuma1102.github.io/cfcompare/reference/autoplot.cf_rmse_curve.md)
  : Line plot of estimation RMSE versus a design dimension
- [`autoplot(`*`<cf_rmse_curves>`*`)`](https://takuma1102.github.io/cfcompare/reference/autoplot.cf_rmse_curves.md)
  : Autoplot for paired RMSE curves
- [`plot(`*`<cf_rmse_curve>`*`)`](https://takuma1102.github.io/cfcompare/reference/plot.cf_rmse_curve.md)
  [`plot(`*`<cf_rmse_curves>`*`)`](https://takuma1102.github.io/cfcompare/reference/plot.cf_rmse_curve.md)
  : Plot one or both estimation-RMSE curves
- [`sim_panel()`](https://takuma1102.github.io/cfcompare/reference/sim_panel.md)
  : Simulate a panel from a low-rank factor model
- [`sim_semisynthetic()`](https://takuma1102.github.io/cfcompare/reference/sim_semisynthetic.md)
  : Build a semi-synthetic panel from real data
- [`true_att()`](https://takuma1102.github.io/cfcompare/reference/true_att.md)
  : True ATT of a (semi-)synthetic panel

## Package overview

- [`cfcompare`](https://takuma1102.github.io/cfcompare/reference/cfcompare-package.md)
  [`cfcompare-package`](https://takuma1102.github.io/cfcompare/reference/cfcompare-package.md)
  : cfcompare: Compare DID, SDID, MC and the Triply Robust Panel
  estimator
- [`reexports`](https://takuma1102.github.io/cfcompare/reference/reexports.md)
  [`autoplot`](https://takuma1102.github.io/cfcompare/reference/reexports.md)
  : Objects exported from other packages
