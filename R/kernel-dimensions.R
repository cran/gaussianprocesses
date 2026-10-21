# Kernels on selected input columns -------------------------------------------
#
# A select node restricts its single child to fixed input columns S:
# k_S(x, x') = k(x_S, x'_S). It has no hyperparameters and is transparent in
# parameter paths, so select_dimensions(rbf_kernel(), 2) has the paths
# `variance` and `length_scale`, exactly like rbf_kernel(). Columns index the
# inputs that the node receives, so nested selections compose.


.validate_selected_columns <- function(columns) {
  if (!is.numeric(columns) || !is.null(dim(columns)) ||
      length(columns) == 0L || any(!is.finite(columns)) ||
      any(columns < 1) || any(columns != round(columns))) {
    .gp_abort(
      paste0(
        "'columns' must be a non-empty vector of positive integer column ",
        "indices."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (anyDuplicated(columns)) {
    .gp_abort(
      "'columns' must not contain duplicate column indices.",
      "gaussianprocesses_validation_error"
    )
  }

  as.integer(columns)
}


.assert_selected_columns_in_range <- function(columns, n_features) {
  outside <- columns[columns > n_features]

  if (length(outside) > 0L) {
    .gp_abort(
      sprintf(
        paste0(
          "select_dimensions() selects column(s) %s, but its inputs have ",
          "%d column(s)."
        ),
        paste(outside, collapse = ", "),
        n_features
      ),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(columns)
}


# The selected columns of an input matrix.
.select_input_columns <- function(kernel, x) {
  .assert_selected_columns_in_range(kernel$columns, ncol(x))
  x[, kernel$columns, drop = FALSE]
}


# Checks that every selection is within the inputs that its node receives and
# that every ARD parameter holds one value or one value per input dimension of
# its own kernel. Models call this before computing any covariance, so that
# errors name the parameter path.
.assert_kernel_dimensions <- function(kernel, n_features, prefix = "") {
  if (length(kernel$children) == 0L) {
    definitions <- .kernel_entry(kernel$type)$parameters

    for (name in names(definitions)) {
      n_values <- length(kernel$parameters[[name]])

      if (identical(definitions[[name]]$shape, "ard") &&
          !(n_values %in% c(1L, n_features))) {
        .gp_abort(
          sprintf(
            paste0(
              "Kernel parameter '%s' has %d values, but its kernel acts on %d ",
              "input dimension(s). It needs one value or one per dimension."
            ),
            paste0(prefix, name),
            n_values,
            n_features
          ),
          "gaussianprocesses_validation_error"
        )
      }
    }

    return(invisible(kernel))
  }

  if (identical(kernel$type, "select")) {
    .assert_selected_columns_in_range(kernel$columns, n_features)
    n_features <- length(kernel$columns)
  }

  if (identical(kernel$type, "changepoint")) {
    .assert_changepoint_column_in_range(kernel$column, n_features)
  }

  child_prefixes <- .child_prefixes(kernel, prefix)

  for (index in seq_along(kernel$children)) {
    .assert_kernel_dimensions(
      kernel$children[[index]],
      n_features,
      child_prefixes[[index]]
    )
  }

  invisible(kernel)
}


# Input derivatives of a selected kernel in all `n_features` input dimensions:
# derivatives along unselected columns are zero, because the covariance does
# not depend on them.
.embed_input_derivatives <- function(derivatives, columns, n_features) {
  zero <- matrix(0, nrow(derivatives$value), ncol(derivatives$value))
  positions <- match(seq_len(n_features), columns)

  embed <- function(matrices) {
    lapply(
      positions,
      function(position) if (is.na(position)) zero else matrices[[position]]
    )
  }

  derivatives$gradient_x <- embed(derivatives$gradient_x)
  derivatives$gradient_y <- embed(derivatives$gradient_y)

  if (!is.null(derivatives$cross)) {
    derivatives$cross <- lapply(
      positions,
      function(position) {
        if (is.na(position)) {
          return(rep(list(zero), n_features))
        }

        embed(derivatives$cross[[position]])
      }
    )
  }

  derivatives
}


#' Restrict a kernel to selected input columns
#'
#' Makes a kernel act on some input columns only:
#' \eqn{k_S(x, x') = k(x_S, x'_S)} for the selected columns \eqn{S}. Sums of
#' restricted kernels give additive models such as
#' \eqn{f(x) = f_1(x_1) + f_2(x_2, x_3)}.
#'
#' @param kernel A Gaussian-process kernel specification.
#' @param columns Positive integer indices of the input columns the kernel
#'   acts on, without duplicates.
#'
#' @details
#' A restricted kernel is positive semidefinite whenever `kernel` is, because
#' it is `kernel` composed with a linear projection. A sum
#' \eqn{\sum_j k_j(x_{S_j}, x'_{S_j})} is the covariance of
#' \eqn{f(x) = \sum_j f_j(x_{S_j})} with independent components \eqn{f_j}.
#'
#' The selection is not a hyperparameter. It adds no parameter and no level
#' to parameter paths: `select_dimensions(rbf_kernel(), 2)` has the paths
#' `variance` and `length_scale`, like `rbf_kernel()`. ARD parameters of
#' `kernel` have one value per selected column.
#'
#' Columns index the inputs the kernel receives. The number of input columns
#' is not known when a specification is built, so column indices and ARD
#' lengths are checked when the kernel is evaluated and when a model is
#' fitted. Derivatives with respect to unselected input columns are zero.
#'
#' @return A composite kernel specification.
#' @examples
#' # An additive model f(x) = f1(x1) + f2(x2, x3).
#' kernel <- sum_kernel(
#'   select_dimensions(rbf_kernel(length_scale = 0.5), 1),
#'   select_dimensions(matern52_kernel(length_scale = c(1, 2)), 2:3)
#' )
#' kernel
#' names(kernel_parameters(kernel, flatten = TRUE))
#'
#' x <- cbind(c(0, 0.5, 1), c(1, 0, -1), c(2, 2, 0))
#' all.equal(
#'   evaluate_kernel(kernel, x),
#'   evaluate_kernel(rbf_kernel(length_scale = 0.5), x[, 1]) +
#'     evaluate_kernel(matern52_kernel(length_scale = c(1, 2)), x[, 2:3])
#' )
#' @template stable
#' @export
select_dimensions <- function(kernel, columns) {
  .assert_kernel(kernel)
  columns <- .validate_selected_columns(columns)

  selected <- .new_composite_kernel("select", list(kernel))
  selected$columns <- columns
  selected
}
