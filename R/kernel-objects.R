# Kernel specifications -------------------------------------------------------

.assert_kernel <- function(kernel, name = "kernel") {
  if (!inherits(kernel, "gaussianprocesses_kernel")) {
    .gp_abort(
      sprintf("'%s' must be a Gaussian-process kernel specification.", name),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(kernel)
}


# Leaf kernels are defined by the kernel registry (R/kernel-registry.R).
# Parameters are stored in registry order; omitted ones take their defaults.
.new_leaf_kernel <- function(type, parameters = list()) {
  parameters <- .complete_kernel_parameters(type, parameters)
  .validate_leaf_kernel_parameters(type, parameters)

  structure(
    list(
      type = type,
      parameters = parameters,
      children = list()
    ),
    class = "gaussianprocesses_kernel"
  )
}


.new_composite_kernel <- function(type, children, parameters = list()) {
  if (length(children) == 0L) {
    .gp_abort(
      "A composite kernel must contain at least one child kernel.",
      "gaussianprocesses_validation_error"
    )
  }

  for (index in seq_along(children)) {
    .assert_kernel(children[[index]], sprintf("kernels[[%d]]", index))
  }

  structure(
    list(
      type = type,
      parameters = parameters,
      children = children
    ),
    class = "gaussianprocesses_kernel"
  )
}


#' Create an RBF kernel specification
#'
#' Kernel specifications are reusable objects that store covariance parameters
#' separately from the observations on which they are evaluated. They are the
#' building blocks for composite kernels and later Gaussian-process fitting.
#'
#' @param variance Positive marginal variance.
#' @param length_scale Positive scalar or numeric vector. A vector supplies one
#'   ARD length scale per input dimension and is validated against the data when
#'   the kernel is evaluated.
#'
#' @return An object of class `gaussianprocesses_kernel`.
#' @seealso [kernel_rbf()] computes the same covariance directly from
#'   parameter values, without a kernel specification.
#' @examples
#' kernel <- rbf_kernel(variance = 1.5, length_scale = c(0.5, 2))
#' kernel
#'
#' evaluate_kernel(kernel, cbind(c(0, 1), c(0, 1)))
#' @template stable
#' @export
rbf_kernel <- function(variance = 1, length_scale = 1) {
  .new_leaf_kernel(
    "rbf",
    list(variance = variance, length_scale = length_scale)
  )
}


#' Create Matérn kernel specifications
#'
#' @inheritParams rbf_kernel
#'
#' @return An object of class `gaussianprocesses_kernel`.
#' @seealso [kernel_matern12()], [kernel_matern32()], and
#'   [kernel_matern52()] compute the same covariances directly from
#'   parameter values.
#' @template stable
#' @name matern_kernel_specs
#' @examples
#' matern52_kernel(variance = 1, length_scale = 0.8)
#'
#' # ARD length scales become indexed parameter paths.
#' kernel_parameters(matern32_kernel(length_scale = c(1, 2)), flatten = TRUE)
NULL


#' @rdname matern_kernel_specs
#' @export
matern12_kernel <- function(variance = 1, length_scale = 1) {
  .new_leaf_kernel(
    "matern12",
    list(variance = variance, length_scale = length_scale)
  )
}


#' @rdname matern_kernel_specs
#' @export
matern32_kernel <- function(variance = 1, length_scale = 1) {
  .new_leaf_kernel(
    "matern32",
    list(variance = variance, length_scale = length_scale)
  )
}


#' @rdname matern_kernel_specs
#' @export
matern52_kernel <- function(variance = 1, length_scale = 1) {
  .new_leaf_kernel(
    "matern52",
    list(variance = variance, length_scale = length_scale)
  )
}


#' Create a rational-quadratic kernel specification
#'
#' @inheritParams rbf_kernel
#' @param alpha Positive scale-mixture parameter.
#'
#' @return An object of class `gaussianprocesses_kernel`.
#' @seealso [kernel_rational_quadratic()] computes the same covariance
#'   directly from parameter values, without a kernel specification.
#' @examples
#' kernel <- rational_quadratic_kernel(length_scale = 0.8, alpha = 2)
#' evaluate_kernel(kernel, c(0, 0.5, 1))
#' @template stable
#' @export
rational_quadratic_kernel <- function(
    variance = 1,
    length_scale = 1,
    alpha = 1) {
  .new_leaf_kernel(
    "rational_quadratic",
    list(
      variance = variance,
      length_scale = length_scale,
      alpha = alpha
    )
  )
}


#' Create a periodic kernel specification
#'
#' @inheritParams rbf_kernel
#' @param period Positive period.
#'
#' @return An object of class `gaussianprocesses_kernel`.
#' @seealso [kernel_periodic()] computes the same covariance directly from
#'   parameter values, without a kernel specification.
#' @examples
#' # Monthly data with a yearly cycle.
#' kernel <- periodic_kernel(length_scale = 0.8, period = 12)
#' evaluate_kernel(kernel, c(0, 6, 12))
#' @template stable
#' @export
periodic_kernel <- function(variance = 1, length_scale = 1, period = 1) {
  .new_leaf_kernel(
    "periodic",
    list(
      variance = variance,
      length_scale = length_scale,
      period = period
    )
  )
}


#' Create a linear kernel specification
#'
#' @param variance Positive variance multiplier.
#'
#' @return An object of class `gaussianprocesses_kernel`.
#' @seealso [kernel_linear()] computes the same covariance directly from
#'   parameter values, without a kernel specification.
#' @examples
#' # A linear trend component.
#' evaluate_kernel(linear_kernel(variance = 0.1), 1:4)
#' @template stable
#' @export
linear_kernel <- function(variance = 1) {
  .new_leaf_kernel("linear", list(variance = variance))
}


#' Create a white-noise kernel specification
#'
#' @param variance Positive white-noise variance.
#'
#' @return An object of class `gaussianprocesses_kernel`.
#' @seealso [kernel_white_noise()] computes the same covariance directly from
#'   parameter values, without a kernel specification.
#' @examples
#' # Add an independent nugget to a smooth kernel.
#' kernel <- sum_kernel(rbf_kernel(), white_noise_kernel(variance = 0.05))
#' evaluate_kernel(kernel, c(0, 0.5))
#' @template stable
#' @export
white_noise_kernel <- function(variance = 1) {
  .new_leaf_kernel("white_noise", list(variance = variance))
}


#' Add covariance kernels
#'
#' The covariance matrices generated by all supplied kernels are added
#' element-wise.
#'
#' @param ... Two or more Gaussian-process kernel specifications.
#'
#' @return A composite kernel specification.
#' @examples
#' kernel <- sum_kernel(linear_kernel(variance = 0.1), periodic_kernel(period = 12))
#' kernel
#'
#' names(kernel_parameters(kernel, flatten = TRUE))
#' @template stable
#' @export
sum_kernel <- function(...) {
  kernels <- list(...)

  if (length(kernels) < 2L) {
    .gp_abort(
      "'sum_kernel' requires at least two kernels.",
      "gaussianprocesses_validation_error"
    )
  }

  .new_composite_kernel("sum", kernels)
}


#' Multiply covariance kernels
#'
#' The covariance matrices generated by all supplied kernels are multiplied
#' element-wise.
#'
#' @inheritParams sum_kernel
#'
#' @return A composite kernel specification.
#' @examples
#' # A locally periodic kernel: a periodic pattern whose shape drifts slowly.
#' kernel <- product_kernel(
#'   periodic_kernel(period = 12),
#'   rbf_kernel(length_scale = 24)
#' )
#' evaluate_kernel(kernel, c(0, 12, 48))[1, ]
#' @template stable
#' @export
product_kernel <- function(...) {
  kernels <- list(...)

  if (length(kernels) < 2L) {
    .gp_abort(
      "'product_kernel' requires at least two kernels.",
      "gaussianprocesses_validation_error"
    )
  }

  .new_composite_kernel("product", kernels)
}


#' Scale a covariance kernel
#'
#' Multiplies a covariance kernel by a positive scalar. The scalar is retained
#' as a named parameter and can therefore be optimized together with parameters
#' of nested kernels.
#'
#' @param kernel A Gaussian-process kernel specification.
#' @param scale Positive scalar multiplier.
#'
#' @return A composite kernel specification.
#' @examples
#' kernel <- scale_kernel(rbf_kernel(), scale = 3)
#' all.equal(
#'   evaluate_kernel(kernel, 1:3),
#'   3 * evaluate_kernel(rbf_kernel(), 1:3)
#' )
#' @template stable
#' @export
scale_kernel <- function(kernel, scale = 1) {
  .assert_kernel(kernel)
  .assert_numeric_scalar(scale, "scale", lower = 0, lower_open = TRUE)

  .new_composite_kernel(
    "scale",
    list(kernel),
    parameters = list(scale = scale)
  )
}


.evaluate_leaf_kernel <- function(kernel, x, y) {
  .evaluate_leaf_parameters(kernel$type, kernel$parameters, x, y)
}


#' Evaluate a kernel specification
#'
#' The package offers two ways to compute a covariance matrix. Kernel
#' specifications, such as [rbf_kernel()], are objects: they can be combined
#' with [sum_kernel()], [product_kernel()], and [scale_kernel()], their
#' parameters can be read, updated, and optimized, and models store them.
#' `evaluate_kernel()` computes their covariance. The covariance functions,
#' such as [kernel_rbf()], take the parameter values as arguments and return
#' one matrix directly. Both use the same formulas.
#'
#' @param kernel A Gaussian-process kernel specification.
#' @param x Numeric vector or matrix of observations.
#' @param y Optional numeric vector or matrix. If `NULL`, evaluates the
#'   covariance of `x` with itself.
#'
#' @return A numeric covariance matrix.
#' @examples
#' kernel <- sum_kernel(rbf_kernel(length_scale = 0.5), white_noise_kernel(0.1))
#' x <- c(0, 0.5, 1)
#'
#' evaluate_kernel(kernel, x)
#'
#' # Cross-covariance between training and new inputs.
#' evaluate_kernel(kernel, x, c(0.25, 2))
#' @template stable
#' @export
evaluate_kernel <- function(kernel, x, y = NULL) {
  .assert_kernel(kernel)

  if (.is_leaf_kernel_type(kernel$type)) {
    return(.evaluate_leaf_kernel(kernel, x, y))
  }

  if (identical(kernel$type, "changepoint")) {
    return(.evaluate_changepoint(kernel, x, y))
  }

  if (identical(kernel$type, "select")) {
    inputs <- .prepare_kernel_inputs(x, y)

    return(
      evaluate_kernel(
        kernel$children[[1L]],
        .select_input_columns(kernel, inputs$x),
        if (!is.null(y)) .select_input_columns(kernel, inputs$y)
      )
    )
  }

  child_matrices <- lapply(
    kernel$children,
    evaluate_kernel,
    x = x,
    y = y
  )

  switch(
    kernel$type,
    sum = Reduce(`+`, child_matrices),
    product = Reduce(`*`, child_matrices),
    scale = kernel$parameters$scale * child_matrices[[1L]],
    .gp_abort(
      sprintf("Unknown composite kernel type '%s'.", kernel$type),
      "gaussianprocesses_validation_error"
    )
  )
}


.kernel_diagonal <- function(kernel, x) {
  .assert_kernel(kernel)
  x_matrix <- .as_input_matrix(x, name = "x")

  if (length(kernel$children) == 0L) {
    entry <- .kernel_entry(kernel$type)
    .validate_leaf_kernel_parameters(kernel$type, kernel$parameters)

    return(entry$diagonal(x_matrix, kernel$parameters))
  }

  if (identical(kernel$type, "select")) {
    return(
      .kernel_diagonal(
        kernel$children[[1L]],
        .select_input_columns(kernel, x_matrix)
      )
    )
  }

  if (identical(kernel$type, "changepoint")) {
    return(.changepoint_diagonal(kernel, x_matrix))
  }

  child_diagonals <- lapply(
    kernel$children,
    .kernel_diagonal,
    x = x_matrix
  )

  switch(
    kernel$type,
    sum = Reduce(`+`, child_diagonals),
    product = Reduce(`*`, child_diagonals),
    scale = kernel$parameters$scale * child_diagonals[[1L]],
    .gp_abort(
      sprintf("Unknown composite kernel type '%s'.", kernel$type),
      "gaussianprocesses_validation_error"
    )
  )
}


#' Evaluate only the diagonal of a kernel matrix
#'
#' Computes `diag(K(x, x))` without constructing the full covariance matrix.
#' This is particularly useful in sparse Gaussian-process approximations.
#'
#' @param kernel A Gaussian-process kernel specification.
#' @param x Numeric vector or matrix of observations.
#'
#' @return A numeric vector with one prior variance per observation.
#' @examples
#' kernel <- sum_kernel(rbf_kernel(variance = 2), linear_kernel())
#' x <- cbind(1:3, 0:2)
#'
#' kernel_diagonal(kernel, x)
#' all.equal(kernel_diagonal(kernel, x), diag(evaluate_kernel(kernel, x)))
#' @template stable
#' @export
kernel_diagonal <- function(kernel, x) {
  .kernel_diagonal(kernel, x)
}


.nested_kernel_parameters <- function(kernel) {
  if (length(kernel$children) == 0L) {
    return(kernel$parameters)
  }

  if (identical(kernel$type, "scale")) {
    return(
      list(
        scale = kernel$parameters$scale,
        kernel = .nested_kernel_parameters(kernel$children[[1L]])
      )
    )
  }

  if (identical(kernel$type, "select")) {
    return(.nested_kernel_parameters(kernel$children[[1L]]))
  }

  if (identical(kernel$type, "changepoint")) {
    return(
      list(
        location = kernel$parameters$location,
        steepness = kernel$parameters$steepness,
        before = .nested_kernel_parameters(kernel$children[[1L]]),
        after = .nested_kernel_parameters(kernel$children[[2L]])
      )
    )
  }

  children <- lapply(kernel$children, .nested_kernel_parameters)
  names(children) <- paste0("kernel", seq_along(children))
  children
}


# Flattened paths of one parameter: `name` for a scalar and `name[1]`,
# `name[2]`, ... for a vector such as ARD length scales.
.parameter_paths <- function(parameter_name, n_values, prefix) {
  if (n_values == 1L) {
    return(paste0(prefix, parameter_name))
  }

  paste0(prefix, parameter_name, "[", seq_len(n_values), "]")
}


# Path prefixes of a composite kernel's children: `kernel.` for the single
# child of a scale kernel, the kernel's own prefix for the single child of a
# select kernel, which adds no level, `before.` and `after.` for the two
# children of a changepoint kernel, and `kernel1.`, `kernel2.`, ...
# otherwise.
.child_prefixes <- function(kernel, prefix) {
  if (identical(kernel$type, "scale")) {
    return(paste0(prefix, "kernel."))
  }

  if (identical(kernel$type, "changepoint")) {
    return(paste0(prefix, c("before.", "after.")))
  }

  if (identical(kernel$type, "select")) {
    return(prefix)
  }

  paste0(prefix, "kernel", seq_along(kernel$children), ".")
}


.flatten_parameter <- function(parameter, parameter_name, prefix) {
  values <- as.numeric(parameter)
  names(values) <- .parameter_paths(parameter_name, length(values), prefix)
  values
}


.flatten_kernel_parameters <- function(kernel, prefix = "") {
  if (length(kernel$children) == 0L) {
    return(
      do.call(
        c,
        lapply(
          names(kernel$parameters),
          function(parameter_name) {
            .flatten_parameter(
              kernel$parameters[[parameter_name]],
              parameter_name,
              prefix
            )
          }
        )
      )
    )
  }

  own <- if (identical(kernel$type, "scale")) {
    .flatten_parameter(kernel$parameters$scale, "scale", prefix)
  } else if (identical(kernel$type, "changepoint")) {
    c(
      .flatten_parameter(kernel$parameters$location, "location", prefix),
      .flatten_parameter(kernel$parameters$steepness, "steepness", prefix)
    )
  }

  children <- Map(
    .flatten_kernel_parameters,
    kernel$children,
    .child_prefixes(kernel, prefix)
  )

  c(own, do.call(c, unname(children)))
}


# Constraint of every flattened parameter, named and ordered like
# kernel_parameters(kernel, flatten = TRUE). ARD values share the constraint
# of their parameter; the scale of a scale kernel is positive.
.flatten_kernel_constraints <- function(kernel, prefix = "") {
  if (length(kernel$children) == 0L) {
    definitions <- .kernel_entry(kernel$type)$parameters

    return(
      do.call(
        c,
        lapply(
          names(kernel$parameters),
          function(parameter_name) {
            paths <- .parameter_paths(
              parameter_name,
              length(kernel$parameters[[parameter_name]]),
              prefix
            )
            constraints <- rep(
              list(definitions[[parameter_name]]$constraint),
              length(paths)
            )
            names(constraints) <- paths
            constraints
          }
        )
      )
    )
  }

  own <- if (identical(kernel$type, "scale")) {
    constraint <- list(.positive_constraint())
    names(constraint) <- .parameter_paths("scale", 1L, prefix)
    constraint
  } else if (identical(kernel$type, "changepoint")) {
    # The location is real-valued and the steepness positive.
    constraints <- list(.real_constraint(), .positive_constraint())
    names(constraints) <- c(
      .parameter_paths("location", 1L, prefix),
      .parameter_paths("steepness", 1L, prefix)
    )
    constraints
  }

  children <- Map(
    .flatten_kernel_constraints,
    kernel$children,
    .child_prefixes(kernel, prefix)
  )

  c(own, do.call(c, unname(children)))
}


#' Extract kernel parameters
#'
#' @param kernel A Gaussian-process kernel specification.
#' @param flatten If `FALSE`, return a nested list matching the kernel
#'   structure. If `TRUE`, return a named numeric vector with stable parameter
#'   paths suitable for optimization. Vector parameters use indexed names such
#'   as `length_scale[1]`.
#'
#' @details
#' Flat parameter names encode the position of each nested kernel. Scalar
#' parameters use paths such as `kernel1.variance`; vector-valued ARD
#' parameters use paths such as `kernel1.length_scale[1]` and
#' `kernel1.length_scale[2]`. These names remain stable when parameter values
#' are updated and can be passed directly to optimization routines.
#'
#' A large ARD length scale means covariance changes slowly along that input
#' dimension, while a small value means covariance changes more rapidly. This is
#' a model sensitivity parameterization and should not be interpreted as a
#' causal importance measure.
#'
#' @return A nested list or named numeric vector.
#' @examples
#' kernel <- sum_kernel(
#'   rbf_kernel(length_scale = c(1, 2)),
#'   white_noise_kernel(variance = 0.1)
#' )
#'
#' str(kernel_parameters(kernel))
#' kernel_parameters(kernel, flatten = TRUE)
#' @template stable
#' @export
kernel_parameters <- function(kernel, flatten = FALSE) {
  .assert_kernel(kernel)

  if (!is.logical(flatten) || length(flatten) != 1L || is.na(flatten)) {
    .gp_abort(
      "'flatten' must be TRUE or FALSE.",
      "gaussianprocesses_validation_error"
    )
  }

  if (flatten) {
    return(.flatten_kernel_parameters(kernel))
  }

  .nested_kernel_parameters(kernel)
}


.rebuild_leaf_parameters <- function(parameters, values, prefix) {
  rebuilt <- parameters

  for (parameter_name in names(parameters)) {
    paths <- .parameter_paths(
      parameter_name,
      length(parameters[[parameter_name]]),
      prefix
    )
    rebuilt[[parameter_name]] <- unname(values[paths])
  }

  rebuilt
}


.rebuild_kernel_parameters <- function(kernel, values, prefix = "") {
  if (length(kernel$children) == 0L) {
    parameters <- .rebuild_leaf_parameters(
      kernel$parameters,
      values,
      prefix
    )

    .validate_leaf_kernel_parameters(kernel$type, parameters)
    kernel$parameters <- parameters
    return(kernel)
  }

  if (identical(kernel$type, "scale")) {
    scale_value <- unname(values[[.parameter_paths("scale", 1L, prefix)]])
    .assert_numeric_scalar(
      scale_value,
      "scale",
      lower = 0,
      lower_open = TRUE
    )
    kernel$parameters$scale <- scale_value
  }

  if (identical(kernel$type, "changepoint")) {
    location <- unname(values[[.parameter_paths("location", 1L, prefix)]])
    steepness <- unname(values[[.parameter_paths("steepness", 1L, prefix)]])
    .assert_numeric_scalar(location, "location")
    .assert_numeric_scalar(steepness, "steepness", lower = 0, lower_open = TRUE)
    kernel$parameters$location <- location
    kernel$parameters$steepness <- steepness
  }

  kernel$children <- Map(
    function(child, child_prefix) {
      .rebuild_kernel_parameters(child, values, child_prefix)
    },
    kernel$children,
    .child_prefixes(kernel, prefix)
  )

  kernel
}


#' Update parameters in a kernel specification
#'
#' Parameters are identified by the names returned by
#' `kernel_parameters(kernel, flatten = TRUE)`. Partial updates are supported,
#' which allows an optimizer or tuning routine to modify selected parameters
#' without manually unpacking nested kernel objects.
#'
#' @param kernel A Gaussian-process kernel specification.
#' @param parameters Named numeric vector of parameter updates.
#'
#' @return A new kernel specification with updated parameters.
#' @examples
#' kernel <- sum_kernel(
#'   rbf_kernel(length_scale = c(1, 2)),
#'   white_noise_kernel(variance = 0.1)
#' )
#'
#' updated <- update_kernel_parameters(
#'   kernel,
#'   c("kernel1.length_scale[2]" = 5, "kernel2.variance" = 0.01)
#' )
#' kernel_parameters(updated, flatten = TRUE)
#' @template stable
#' @export
update_kernel_parameters <- function(kernel, parameters) {
  .assert_kernel(kernel)

  if (!is.numeric(parameters) || !is.null(dim(parameters)) ||
      length(parameters) == 0L) {
    .gp_abort(
      "'parameters' must be a non-empty named numeric vector.",
      "gaussianprocesses_validation_error"
    )
  }

  parameter_names <- names(parameters)

  if (is.null(parameter_names) || any(parameter_names == "")) {
    .gp_abort(
      "'parameters' must have non-empty names.",
      "gaussianprocesses_validation_error"
    )
  }

  if (anyDuplicated(parameter_names)) {
    .gp_abort(
      "'parameters' contains duplicate names.",
      "gaussianprocesses_validation_error"
    )
  }

  if (any(!is.finite(parameters))) {
    .gp_abort(
      "'parameters' must contain only finite values.",
      "gaussianprocesses_validation_error"
    )
  }

  current <- kernel_parameters(kernel, flatten = TRUE)
  unknown <- setdiff(parameter_names, names(current))

  if (length(unknown) > 0L) {
    .gp_abort(
      sprintf(
        "Unknown kernel parameter path(s): %s.",
        paste(unknown, collapse = ", ")
      ),
      "gaussianprocesses_validation_error"
    )
  }

  current[parameter_names] <- parameters

  .rebuild_kernel_parameters(kernel, current)
}


.kernel_display_name <- function(type) {
  if (.is_leaf_kernel_type(type)) {
    return(.kernel_entry(type)$label)
  }

  switch(
    type,
    sum = "SumKernel",
    product = "ProductKernel",
    scale = "ScaleKernel",
    select = "SelectDimensions",
    changepoint = "ChangepointKernel",
    type
  )
}


.format_kernel_parameter_value <- function(value) {
  # Format elements separately: format() on a vector pads every element to a
  # shared number of decimals, printing c(0.5, 2) as c(0.5, 2.0).
  formatted <- vapply(
    signif(value, 6),
    format,
    character(1),
    trim = TRUE
  )

  if (length(formatted) == 1L) {
    return(formatted)
  }

  paste0("c(", paste(formatted, collapse = ", "), ")")
}


.format_kernel_parameters <- function(parameters) {
  paste(
    vapply(
      names(parameters),
      function(parameter_name) {
        paste0(
          parameter_name,
          "=",
          .format_kernel_parameter_value(
            parameters[[parameter_name]]
          )
        )
      },
      character(1)
    ),
    collapse = ", "
  )
}


.format_kernel <- function(kernel, level = 0L) {
  indentation <- paste(rep("  ", level), collapse = "")
  display_name <- .kernel_display_name(kernel$type)

  if (length(kernel$children) == 0L) {
    return(
      paste0(
        indentation,
        display_name,
        "(",
        .format_kernel_parameters(kernel$parameters),
        ")"
      )
    )
  }

  if (identical(kernel$type, "scale")) {
    header <- paste0(
      indentation,
      display_name,
      "(scale=",
      format(signif(kernel$parameters$scale, 6), trim = TRUE),
      ","
    )
  } else if (identical(kernel$type, "select")) {
    header <- paste0(
      indentation,
      display_name,
      "(columns=",
      .format_kernel_parameter_value(kernel$columns),
      ","
    )
  } else if (identical(kernel$type, "changepoint")) {
    header <- paste0(
      indentation,
      display_name,
      "(",
      .format_kernel_parameters(kernel$parameters),
      ", column=",
      kernel$column,
      ","
    )
  } else {
    header <- paste0(indentation, display_name, "(")
  }

  child_lines <- unlist(
    lapply(
      kernel$children,
      .format_kernel,
      level = level + 1L
    ),
    use.names = FALSE
  )

  c(header, child_lines, paste0(indentation, ")"))
}


#' @export
print.gaussianprocesses_kernel <- function(x, ...) {
  cat(paste(.format_kernel(x), collapse = "\n"), "\n", sep = "")
  invisible(x)
}
