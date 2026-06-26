#' cfcompare: Compare DID, SDID, MC and the Triply Robust Panel estimator
#'
#' `cfcompare` lets applied researchers run several panel causal-inference estimators
#' on the **same** data and compare them on one tidy schema and one plot. It
#' provides a native R implementation of the Triply RObust Panel (TROP)
#' estimator of Athey, Imbens, Qu & Viviano (2026), and a comparison harness that
#' also wraps difference-in-differences (TWFE), matrix completion, and -- via the
#' \pkg{synthdid} package -- synthetic difference-in-differences and synthetic
#' control.
#'
#' Key functions:
#' \itemize{
#'   \item [trop()] -- the native TROP estimator.
#'   \item [panel_compare()] -- run DID / SDID / MC / TROP together.
#'   \item [as_att()] -- coerce external results into the shared schema.
#'   \item [autoplot.cf_comparison()], [plot_counterfactual()] -- visuals.
#'   \item [sim_panel()] -- simulate factor-model panels for experiments.
#' }
#'
#' @references Athey, S., Imbens, G. W., Qu, Z., & Viviano, D. (2026).
#'   Triply Robust Panel Estimators. arXiv:2508.21536.
#'
#' @keywords internal
#' @importFrom stats as.formula median na.omit qnorm rnorm sd vcov
#' @importFrom utils modifyList
#' @importFrom ggplot2 autoplot .data
"_PACKAGE"

# Re-export the ggplot2 generic so autoplot() works without attaching ggplot2.
#' @export
ggplot2::autoplot
