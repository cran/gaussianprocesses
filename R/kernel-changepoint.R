# Changepoint kernels ---------------------------------------------------------
#
# With the sigmoid s(x) = 1 / (1 + exp(-a (x_t - c))) along column t, which
# rises from 0 before c to 1 after it,
#   k(x, x') = (1 - s(x)) k_1(x, x') (1 - s(x')) + s(x) k_2(x, x') s(x'),
# with k_1 the kernel before the change and k_2 the kernel after it.
# Each term is g(x) k_i(x, x') g(x'), positive semidefinite, so the sum is.
# The location c is real (identity coordinate) and the steepness a positive
# (log coordinate):
#   ds / dc = -a s (1 - s),        ds / d log a = a (x_t - c) s (1 - s),
#   ds / dx_t = a s (1 - s) =: v,  dv / dc = -a (1 - 2 s) v,
#   dv / d log a = v + a (1 - 2 s) (x_t - c) v.
# The weights (1 - s(x)) (1 - s(x')) and s(x) s(x') are rank-one "kernels"
# whose derivative tuples (R/kernel-input-derivatives.R) are outer products,
# so input and hyperparameter derivatives follow from the product rule.


.validate_changepoint_column <- function(column) {
  if (!is.numeric(column) || length(column) != 1L || !is.finite(column) ||
      column < 1 || column != round(column)) {
    .gp_abort(
      "'column' must be a single positive integer column index.",
      "gaussianprocesses_validation_error"
    )
  }

  as.integer(column)
}


.assert_changepoint_column_in_range <- function(column, n_features) {
  if (column > n_features) {
    .gp_abort(
      sprintf(
        paste0(
          "changepoint_kernel() acts along column %d, but its inputs have ",
          "%d column(s)."
        ),
        column,
        n_features
      ),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(column)
}


# The sigmoid weight s of the `after` kernel at the rows of `x`, its slope
# along the changepoint column, and their derivatives with respect to the
# optimizer coordinates of the location and the steepness.
.changepoint_sigmoid <- function(kernel, x) {
  .assert_changepoint_column_in_range(kernel$column, ncol(x))
  steepness <- kernel$parameters$steepness
  offset <- x[, kernel$column] - kernel$parameters$location
  z <- steepness * offset
  value <- stats::plogis(z)
  # 1 - s and s (1 - s) without cancellation in the tails.
  complement <- stats::plogis(-z)
  spread <- value * complement
  slope <- steepness * spread
  curvature <- steepness * (1 - 2 * value)

  list(
    value = value,
    complement = complement,
    slope = slope,
    location = list(value = -slope, slope = -curvature * slope),
    steepness = list(
      value = offset * slope,
      slope = slope + curvature * offset * slope
    )
  )
}


# The derivative tuple of w_x(x) w_y(y) when both weights depend only on
# column `column`, with slopes dw / dx_column.
.outer_weight_tuple <- function(
    weight_x,
    slope_x,
    weight_y,
    slope_y,
    column,
    n_features,
    order) {
  zero <- matrix(0, length(weight_x), length(weight_y))
  gradient_x <- rep(list(zero), n_features)
  gradient_y <- rep(list(zero), n_features)
  gradient_x[[column]] <- outer(slope_x, weight_y)
  gradient_y[[column]] <- outer(weight_x, slope_y)

  tuple <- list(
    value = outer(weight_x, weight_y),
    gradient_x = gradient_x,
    gradient_y = gradient_y
  )

  if (order == 2L) {
    cross <- rep(list(rep(list(zero), n_features)), n_features)
    cross[[column]][[column]] <- outer(slope_x, slope_y)
    tuple$cross <- cross
  }

  tuple
}


# Weight tuples of the two terms, and their derivatives with respect to the
# location and the steepness.
.changepoint_weight_tuples <- function(kernel, x, y, order) {
  sigmoid_x <- .changepoint_sigmoid(kernel, x)
  sigmoid_y <- .changepoint_sigmoid(kernel, y)
  n_features <- ncol(x)

  weight <- function(sigmoid, before) {
    if (!before) {
      return(sigmoid)
    }

    # 1 - s and its derivatives.
    list(
      value = sigmoid$complement,
      slope = -sigmoid$slope,
      location = lapply(sigmoid$location, `-`),
      steepness = lapply(sigmoid$steepness, `-`)
    )
  }

  tuples <- function(before) {
    w_x <- weight(sigmoid_x, before)
    w_y <- weight(sigmoid_y, before)
    tuple <- function(a, a_slope, b, b_slope) {
      .outer_weight_tuple(
        a,
        a_slope,
        b,
        b_slope,
        kernel$column,
        n_features,
        order
      )
    }
    derivative <- function(parameter) {
      .add_derivatives(
        tuple(
          w_x[[parameter]]$value,
          w_x[[parameter]]$slope,
          w_y$value,
          w_y$slope
        ),
        tuple(
          w_x$value,
          w_x$slope,
          w_y[[parameter]]$value,
          w_y[[parameter]]$slope
        )
      )
    }

    list(
      tuple = tuple(w_x$value, w_x$slope, w_y$value, w_y$slope),
      location = derivative("location"),
      steepness = derivative("steepness")
    )
  }

  list(before = tuples(TRUE), after = tuples(FALSE))
}


.evaluate_changepoint <- function(kernel, x, y) {
  inputs <- .prepare_kernel_inputs(x, y)
  sigmoid_x <- .changepoint_sigmoid(kernel, inputs$x)
  sigmoid_y <- if (is.null(y)) {
    sigmoid_x
  } else {
    .changepoint_sigmoid(kernel, inputs$y)
  }

  outer(sigmoid_x$complement, sigmoid_y$complement) *
    evaluate_kernel(kernel$children[[1L]], x, y) +
    outer(sigmoid_x$value, sigmoid_y$value) *
      evaluate_kernel(kernel$children[[2L]], x, y)
}


.changepoint_diagonal <- function(kernel, x) {
  sigmoid <- .changepoint_sigmoid(kernel, x)

  sigmoid$complement^2 * .kernel_diagonal(kernel$children[[1L]], x) +
    sigmoid$value^2 * .kernel_diagonal(kernel$children[[2L]], x)
}


# Covariance and hyperparameter derivatives (R/kernel-gradients.R).
.changepoint_value_and_gradient <- function(kernel, x, y, prefix) {
  prefixes <- .child_prefixes(kernel, prefix)
  children <- Map(
    function(child, child_prefix) {
      .kernel_value_and_gradient(child, x, y, child_prefix)
    },
    kernel$children,
    prefixes
  )
  weights <- .changepoint_weight_tuples(kernel, x, y, order = 1L)
  before <- weights$before
  after <- weights$after
  value_before <- children[[1L]]$value
  value_after <- children[[2L]]$value

  own <- lapply(
    c(location = "location", steepness = "steepness"),
    function(parameter) {
      before[[parameter]]$value * value_before +
        after[[parameter]]$value * value_after
    }
  )
  names(own) <- c(
    .parameter_paths("location", 1L, prefix),
    .parameter_paths("steepness", 1L, prefix)
  )

  list(
    value = before$tuple$value * value_before +
      after$tuple$value * value_after,
    gradient = c(
      own,
      lapply(children[[1L]]$gradient, function(d) before$tuple$value * d),
      lapply(children[[2L]]$gradient, function(d) after$tuple$value * d)
    )
  )
}


# Derivative tuple (R/kernel-input-derivatives.R).
.changepoint_input_derivatives <- function(kernel, x, y, order) {
  weights <- .changepoint_weight_tuples(kernel, x, y, order)

  .add_derivatives(
    .multiply_derivatives(
      weights$before$tuple,
      .kernel_input_derivatives(kernel$children[[1L]], x, y, order)
    ),
    .multiply_derivatives(
      weights$after$tuple,
      .kernel_input_derivatives(kernel$children[[2L]], x, y, order)
    )
  )
}


# Derivative tuple and its hyperparameter derivatives
# (R/derivative-observations.R).
.changepoint_derivative_tuples <- function(kernel, x, y, prefix) {
  prefixes <- .child_prefixes(kernel, prefix)
  before_child <- .kernel_derivative_tuples(
    kernel$children[[1L]],
    x,
    y,
    prefixes[[1L]]
  )
  after_child <- .kernel_derivative_tuples(
    kernel$children[[2L]],
    x,
    y,
    prefixes[[2L]]
  )
  weights <- .changepoint_weight_tuples(kernel, x, y, order = 2L)
  combine <- function(before_weight, after_weight) {
    .add_derivatives(
      .multiply_derivatives(before_weight, before_child$tuple),
      .multiply_derivatives(after_weight, after_child$tuple)
    )
  }

  own <- list(
    combine(weights$before$location, weights$after$location),
    combine(weights$before$steepness, weights$after$steepness)
  )
  names(own) <- c(
    .parameter_paths("location", 1L, prefix),
    .parameter_paths("steepness", 1L, prefix)
  )

  list(
    tuple = combine(weights$before$tuple, weights$after$tuple),
    gradient = c(
      own,
      lapply(
        before_child$gradient,
        .multiply_derivatives,
        first = weights$before$tuple
      ),
      lapply(
        after_child$gradient,
        .multiply_derivatives,
        first = weights$after$tuple
      )
    )
  )
}


# Columns of the training inputs along which each changepoint location of a
# kernel acts, named by parameter path. `columns` maps the columns a node
# receives to columns of the training inputs.
.changepoint_location_columns <- function(kernel, columns, prefix = "") {
  if (length(kernel$children) == 0L) {
    return(integer(0))
  }

  own <- integer(0)

  if (identical(kernel$type, "select")) {
    columns <- columns[kernel$columns]
  }

  if (identical(kernel$type, "changepoint")) {
    own <- stats::setNames(
      columns[[kernel$column]],
      .parameter_paths("location", 1L, prefix)
    )
  }

  child_prefixes <- .child_prefixes(kernel, prefix)

  c(
    own,
    unlist(
      Map(
        function(child, child_prefix) {
          .changepoint_location_columns(child, columns, child_prefix)
        },
        kernel$children,
        child_prefixes
      )
    )
  )
}


#' Changepoint kernels
#'
#' Switches from one covariance to another along one input column, smoothly
#' around an estimable location: a structural break in a time series, or a
#' regime change in space.
#'
#' @param before,after Kernel specifications that hold before and after
#'   the change.
#' @param location Finite location \eqn{c} of the change, in the units of
#'   the input column.
#' @param steepness Positive steepness \eqn{a} of the transition, in inverse
#'   input units: the transition takes about \eqn{4 / a} from 12% to 88%.
#' @param column The input column along which the change happens.
#'
#' @details
#' With the sigmoid \eqn{s(x) = 1 / (1 + e^{-a (x_t - c)})} along column
#' \eqn{t},
#' \deqn{k(x, x') = (1 - s(x)) k_1(x, x') (1 - s(x')) +
#' s(x) k_2(x, x') s(x'),}
#' where \eqn{k_1} is `before` and \eqn{k_2} is `after`. Each term has the
#' form \eqn{g(x) k(x, x') g(x')}, so the kernel is positive semidefinite.
#' The two regimes are independent processes: as \eqn{a \to \infty}, the
#' covariance between points on opposite sides of \eqn{c} tends to 0.
#'
#' `location` is a real-valued parameter, optimized on the identity scale,
#' and `steepness` a positive one, optimized on the log scale. Parameter
#' paths are `location` and `steepness` for the change itself, and
#' `before.` and `after.` for the two kernels, so nested changepoints model
#' several changes. [optimize_gp()] bounds every changepoint location by the
#' training range of its column and, because the likelihood is often
#' multimodal in the location, starts it at evenly spaced quantiles of that
#' column.
#'
#' The kernel is as smooth as its two components, so it supports
#' [kernel_input_gradient()], [predict_gradient_gp()], and derivative
#' observations when they do.
#'
#' @return A composite kernel specification.
#' @seealso [time_series_kernel()], whose `changepoint` argument adds a
#'   change to the time-series kernel.
#' @examples
#' kernel <- changepoint_kernel(
#'   before = rbf_kernel(variance = 1, length_scale = 2),
#'   after = rbf_kernel(variance = 0.2, length_scale = 0.5),
#'   location = 5,
#'   steepness = 4
#' )
#' kernel
#' names(kernel_parameters(kernel, flatten = TRUE))
#'
#' # Points on opposite sides of a steep change are nearly uncorrelated.
#' evaluate_kernel(kernel, c(4, 6))
#' @template experimental
#' @export
changepoint_kernel <- function(
    before,
    after,
    location = 0,
    steepness = 1,
    column = 1L) {
  .assert_kernel(before, "before")
  .assert_kernel(after, "after")
  .assert_numeric_scalar(location, "location")
  .assert_numeric_scalar(steepness, "steepness", lower = 0, lower_open = TRUE)
  column <- .validate_changepoint_column(column)

  kernel <- .new_composite_kernel(
    "changepoint",
    list(before, after),
    parameters = list(location = location, steepness = steepness)
  )
  kernel$column <- column
  kernel
}
