# Panel data utilities: validation and long <-> wide (N x T matrix) reshaping.

#' Reshape a long panel into outcome/treatment matrices
#'
#' Converts a long `data.frame` into the `N x T` outcome (`Y`) and treatment
#' (`W`) matrices used throughout the package, after validating that the panel
#' is well formed.
#'
#' @param data A long `data.frame`/`tibble`.
#' @param outcome,treatment,unit,time Column names (strings).
#' @return A list with matrices `Y`, `W` (rows = units, columns = times), and
#'   the sorted `units` and `times` labels.
#' @keywords internal
#' @noRd
.panel_to_matrices <- function(data, outcome, treatment, unit, time) {
  cols <- c(outcome, treatment, unit, time)
  miss <- setdiff(cols, names(data))
  if (length(miss)) {
    stop("Columns not found in `data`: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  data <- as.data.frame(data)
  u <- data[[unit]]
  ti <- data[[time]]
  units <- sort(unique(u))
  times <- sort(unique(ti))
  N <- length(units); Tt <- length(times)

  if (anyDuplicated(data[c(unit, time)])) {
    stop("`data` has duplicate unit-time rows.", call. = FALSE)
  }

  w_vals <- data[[treatment]]
  if (!all(stats::na.omit(w_vals) %in% c(0, 1))) {
    stop("`treatment` must be a 0/1 indicator.", call. = FALSE)
  }

  ri <- match(u, units)
  ci <- match(ti, times)
  Y <- matrix(NA_real_, N, Tt, dimnames = list(units, as.character(times)))
  W <- matrix(0, N, Tt, dimnames = list(units, as.character(times)))
  idx <- cbind(ri, ci)
  Y[idx] <- as.numeric(data[[outcome]])
  W[idx] <- as.numeric(w_vals)
  W[is.na(W)] <- 0

  if (all(W == 0)) stop("No treated cells found.", call. = FALSE)
  if (all(W == 1)) stop("No control cells found.", call. = FALSE)

  list(Y = Y, W = W, units = units, times = times)
}

#' Reshape a long panel into TROP matrices
#'
#' Converts a long panel `data.frame` into the `N x T` outcome (`Y`) and
#' treatment (`W`) matrices used by [trop_matrix()], and reports the treated rows
#' and the width of the treated block. It is the data-frame-to-matrix companion
#' to [trop_matrix()]: rather than building `Y`, `W`, `treated_units` and
#' `treated_periods` by hand, call `panel_matrices()` and pass the pieces
#' straight through.
#'
#' @param data A long `data.frame`/`tibble` with one row per unit-time cell.
#' @param outcome,treatment,unit,time Column names (strings). `treatment` must be
#'   a 0/1 indicator.
#' @return A list with components:
#'   \describe{
#'     \item{`Y`}{`N x T` outcome matrix (rows = units, columns = times).}
#'     \item{`W`}{`N x T` 0/1 treatment matrix on the same grid.}
#'     \item{`units`,`times`}{the sorted unit and time labels.}
#'     \item{`treated_units`}{integer row indices of units with any treated cell.}
#'     \item{`treated_periods`}{number of periods in which any treated unit is
#'       active, i.e. the width of the post block.}
#'   }
#' @seealso [trop_matrix()], [trop()]
#' @export
#' @examples
#' df <- sim_panel(N = 14, T = 9, n_treated = 3, t0 = 6, att = 2, seed = 1)
#' pm <- panel_matrices(df, "y", "w", "id", "t")
#' str(pm)
#' # Feed straight into the matrix-in entry point (treated rows/periods inferred).
#' trop_matrix(pm$Y, pm$W, lambda_unit = 0.1, lambda_time = 0.1, lambda_nn = Inf)
panel_matrices <- function(data, outcome, treatment, unit, time) {
  pm <- .panel_to_matrices(data, outcome, treatment, unit, time)
  tu <- which(rowSums(pm$W) > 0)
  tp <- sum(colSums(pm$W[tu, , drop = FALSE]) > 0)
  list(
    Y = pm$Y, W = pm$W,
    units = pm$units, times = pm$times,
    treated_units = tu, treated_periods = tp
  )
}

#' Reshape long covariate columns into a list of N x T matrices
#'
#' Aligns each covariate column to the same unit/time grid as the outcome and
#' treatment matrices produced by `.panel_to_matrices()`.
#'
#' @param data A long `data.frame`.
#' @param covariates Character vector of covariate column names (or `NULL`).
#' @param unit,time Column names (strings).
#' @param units,times The sorted unit/time labels defining the grid.
#' @return `NULL` if no covariates, else a named list of N x T numeric matrices.
#' @keywords internal
#' @noRd
.covariate_matrices <- function(data, covariates, unit, time, units, times) {
  if (is.null(covariates) || !length(covariates)) return(NULL)
  data <- as.data.frame(data)
  miss <- setdiff(covariates, names(data))
  if (length(miss)) {
    stop("Covariate column(s) not found in `data`: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  N <- length(units); Tt <- length(times)
  idx <- cbind(match(data[[unit]], units), match(data[[time]], times))
  out <- lapply(covariates, function(cn) {
    v <- data[[cn]]
    if (!is.numeric(v)) {
      stop("Covariate `", cn, "` must be numeric.", call. = FALSE)
    }
    M <- matrix(NA_real_, N, Tt,
                dimnames = list(as.character(units), as.character(times)))
    M[idx] <- as.numeric(v)
    if (anyNA(M)) {
      stop("Covariate `", cn, "` has missing unit-time cells; covariates must ",
           "be fully observed.", call. = FALSE)
    }
    M
  })
  stats::setNames(out, covariates)
}

#' Describe the assignment pattern of a treatment matrix
#'
#' Determines whether a treatment matrix follows a synthetic-control-style
#' block design (a set of units treated from a common period onward, never
#' reverting) or a general / staggered / non-absorbing pattern.
#'
#' @param W An `N x T` 0/1 treatment matrix.
#' @return A list describing the design: `type` (`"block"` or `"general"`),
#'   treated unit indices, first treated period, and counts.
#' @keywords internal
#' @noRd
.assignment_pattern <- function(W) {
  N <- nrow(W); Tt <- ncol(W)
  treated_units <- which(rowSums(W) > 0)
  ever_treated_periods <- which(colSums(W) > 0)
  first_treat_period <- if (length(ever_treated_periods)) {
    min(ever_treated_periods)
  } else NA_integer_

  is_block <- TRUE
  t0 <- NA_integer_
  if (length(treated_units)) {
    # block: each treated unit is control before some t0 and treated from t0 on,
    # with a common t0 across treated units.
    starts <- vapply(treated_units, function(i) {
      tr <- which(W[i, ] == 1)
      st <- min(tr)
      # absorbing from st onward?
      if (!all(W[i, st:Tt] == 1)) return(NA_integer_)
      if (st > 1 && any(W[i, 1:(st - 1)] == 1)) return(NA_integer_)
      st
    }, integer(1))
    if (anyNA(starts) || length(unique(starts)) != 1L) {
      is_block <- FALSE
    } else {
      t0 <- unique(starts)
    }
  }

  list(
    type = if (is_block) "block" else "general",
    treated_units = treated_units,
    n_treated_units = length(treated_units),
    n_treated_cells = sum(W == 1),
    first_treat_period = first_treat_period,
    block_t0 = t0
  )
}

# Resolve the set of estimators to run for the comparison entry points
# (panel_compare(), panel_rmse(), rmse_curve()). De-duplicates `methods`,
# validates them against `known`, then drops anything named in `exclude`. Keeping
# this in one place means the `methods` / `exclude` semantics and error messages
# stay identical across the user-facing functions.
.resolve_methods <- function(methods, exclude, known) {
  methods <- unique(methods)
  bad <- setdiff(methods, known)
  if (length(bad)) {
    stop("Unknown method(s): ", paste(bad, collapse = ", "),
         ". Choose from ", paste(known, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.null(exclude)) {
    exclude <- unique(exclude)
    exbad <- setdiff(exclude, known)
    if (length(exbad)) {
      warning("Ignoring unknown method(s) in `exclude`: ",
              paste(exbad, collapse = ", "), call. = FALSE)
    }
    methods <- setdiff(methods, exclude)
    if (!length(methods)) {
      stop("`exclude` removed every requested method; nothing to run.",
           call. = FALSE)
    }
  }
  methods
}
