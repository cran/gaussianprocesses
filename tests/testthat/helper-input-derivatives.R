# Helpers for checking kernel derivatives with respect to the inputs against
# central finite differences.

# Step in the inputs, and the tolerance on the discrepancy relative to the
# largest analytical derivative of the same order. The truncation error of a
# central difference is O(step^2) for smooth kernels and O(step) for the
# Matérn-3/2 cross derivative at coincident inputs; the rounding error is
# about 1e-16 / step.
input_fd_step <- 1e-7
input_fd_tolerance <- 1e-6


shift_input_column <- function(x, dimension, step) {
  x[, dimension] <- x[, dimension] + step
  x
}


# dK(x, y) / dx_d for every d.
finite_difference_input_gradient <- function(
    kernel,
    x,
    y,
    step = input_fd_step) {
  lapply(
    seq_len(ncol(x)),
    function(d) {
      (evaluate_kernel(kernel, shift_input_column(x, d, step), y) -
        evaluate_kernel(kernel, shift_input_column(x, d, -step), y)) /
        (2 * step)
    }
  )
}


# dK(x, y) / dy_e for every e.
finite_difference_input_gradient_y <- function(
    kernel,
    x,
    y,
    step = input_fd_step) {
  lapply(
    seq_len(ncol(y)),
    function(e) {
      (evaluate_kernel(kernel, x, shift_input_column(y, e, step)) -
        evaluate_kernel(kernel, x, shift_input_column(y, e, -step))) /
        (2 * step)
    }
  )
}


# d2K(x, y) / dx_d dy_e, from central differences in y_e of the analytical
# dK / dx_d, which the first-order check verifies separately.
finite_difference_input_cross <- function(
    kernel,
    x,
    y,
    step = input_fd_step) {
  dimensions <- seq_len(ncol(x))
  differences <- lapply(
    dimensions,
    function(e) {
      Map(
        function(plus, minus) (plus - minus) / (2 * step),
        kernel_input_gradient(kernel, x, shift_input_column(y, e, step)),
        kernel_input_gradient(kernel, x, shift_input_column(y, e, -step))
      )
    }
  )

  lapply(dimensions, function(d) lapply(dimensions, function(e) {
    differences[[e]][[d]]
  }))
}


# Largest discrepancy relative to the largest analytical derivative, over all
# matrices of a (possibly nested) list.
relative_derivative_discrepancy <- function(analytical, numerical) {
  analytical <- unlist(analytical)
  numerical <- unlist(numerical)

  max(abs(analytical - numerical)) / max(abs(analytical))
}


# What every differentiable kernel must satisfy: dK / dx_d, dK / dy_e, and
# d2K / dx_d dy_e match central differences, and the exported gradient equals
# the internal one.
expect_input_derivative_contract <- function(kernel, x, y, label) {
  derivatives <- gaussianprocesses:::.kernel_input_derivatives(
    kernel,
    x,
    y,
    order = 2L
  )

  expect_equal(derivatives$value, evaluate_kernel(kernel, x, y), info = label)
  expect_identical(
    kernel_input_gradient(kernel, x, y),
    derivatives$gradient_x,
    info = label
  )

  checks <- list(
    gradient_x = finite_difference_input_gradient(kernel, x, y),
    gradient_y = finite_difference_input_gradient_y(kernel, x, y),
    cross = finite_difference_input_cross(kernel, x, y)
  )

  for (name in names(checks)) {
    expect_lte(
      relative_derivative_discrepancy(derivatives[[name]], checks[[name]]),
      input_fd_tolerance,
      label = paste(label, name)
    )
  }
}


# Central differences of the derivative tuple (value, gradient_x,
# gradient_y, cross) of K(x, y) in each parameter's optimizer coordinate,
# flattened to a vector per parameter path.
finite_difference_tuple_gradient <- function(kernel, x, y, step = 1e-6) {
  parameters <- kernel_parameters(kernel, flatten = TRUE)
  constraints <- gaussianprocesses:::.flatten_kernel_constraints(kernel)
  parts <- c("value", "gradient_x", "gradient_y", "cross")

  tuple_at <- function(path, direction) {
    constraint <- constraints[[path]]
    shifted <- parameters[path]
    shifted[] <- constraint$from_coordinate(
      constraint$to_coordinate(parameters[[path]]) + direction * step
    )
    tuple <- gaussianprocesses:::.kernel_input_derivatives(
      update_kernel_parameters(kernel, shifted),
      x,
      y,
      order = 2L
    )
    unlist(tuple[parts])
  }

  derivatives <- lapply(
    names(parameters),
    function(path) (tuple_at(path, 1) - tuple_at(path, -1)) / (2 * step)
  )
  names(derivatives) <- names(parameters)
  derivatives
}


# The hyperparameter derivatives of a derivative tuple, which the likelihood
# gradient of derivative observations uses, match central differences.
expect_input_parameter_gradient_contract <- function(kernel, x, y, label) {
  analytical <- gaussianprocesses:::.kernel_derivative_tuples(
    kernel,
    x,
    y
  )$gradient
  numerical <- finite_difference_tuple_gradient(kernel, x, y)

  expect_identical(names(analytical), names(numerical), info = label)

  for (path in names(numerical)) {
    values <- unlist(
      analytical[[path]][c("value", "gradient_x", "gradient_y", "cross")]
    )

    expect_lte(
      max(abs(values - numerical[[path]])) / max(1, max(abs(values))),
      1e-6,
      label = paste(label, path)
    )
  }
}


# Covariance of (f(x), df(x) / dx_1, ..., df(x) / dx_D) stacked by blocks.
joint_value_gradient_covariance <- function(kernel, x) {
  derivatives <- gaussianprocesses:::.kernel_input_derivatives(
    kernel,
    x,
    x,
    order = 2L
  )
  dimensions <- seq_len(ncol(x))

  top <- do.call(cbind, c(list(derivatives$value), derivatives$gradient_y))
  bottom <- lapply(
    dimensions,
    function(d) {
      do.call(
        cbind,
        c(list(derivatives$gradient_x[[d]]), derivatives$cross[[d]])
      )
    }
  )

  do.call(rbind, c(list(top), bottom))
}
