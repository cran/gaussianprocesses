# Kernel registry -------------------------------------------------------------
#
# Every leaf kernel is defined by exactly one registry entry. Constructors,
# evaluate_kernel(), kernel_diagonal(), kernel_gradient(),
# kernel_input_gradient(), parameter handling, printing, the FITC diagonal,
# and optimize_gp() all read leaf information from here. Sum, product, scale,
# select, and changepoint kernels stay structural (see R/kernel-objects.R,
# R/kernel-dimensions.R, and R/kernel-changepoint.R) and are not registry
# entries.
#
# An entry is a list with these fields:
#
#   type            Identifier stored in `kernel$type`.
#   label           Display name used by print().
#   parameters      Named list of parameter definitions from
#                   .kernel_parameter(), in constructor order: default value,
#                   shape ("scalar", or "ard" for a scalar or one value per
#                   input dimension), and constraint.
#   stationary      TRUE if k(x, y) depends only on x - y.
#   evaluate        function(x, y, parameters): the covariance matrix between
#                   the rows of input matrices `x` and `y`.
#   diagonal        function(x, parameters): diag(evaluate(x, x, parameters))
#                   without forming the matrix.
#   gradient        function(x, y, parameters, covariance): named list with
#                   dK / d eta for every parameter, where eta is the
#                   parameter's optimizer coordinate (see the constraints
#                   below) and `covariance` is evaluate(x, y, parameters). An
#                   ARD parameter gives one matrix per input dimension, or a
#                   single matrix when one value is shared by all dimensions.
#   smoothness      Number of mean-square derivatives of the process: 0
#                   (not differentiable), a positive integer, or Inf.
#   input_gradient  function(x, y, parameters): list over input dimensions d
#                   of dk(x_i, y_j) / dx_d, or NULL when smoothness is 0.
#   input_hessian   function(x, y, parameters): list over d of lists over e
#                   of d2k(x_i, y_j) / dx_d dy_e, or NULL when smoothness
#                   is 0. See R/kernel-input-derivatives.R.
#   input_parameter_gradient
#                   function(x, y, parameters): the derivatives of the two
#                   lists above with respect to each parameter's optimizer
#                   coordinate, for derivative observations; NULL when
#                   smoothness is 0. See R/kernel-input-derivatives.R.
#   state_space     function(parameters): the stochastic differential
#                   equation of a one-dimensional process, or NULL when the
#                   kernel has no exact finite-dimensional one. See
#                   .matern_state_space() in R/gp-state-space.R for the
#                   fields.
#
# Parameter values are validated against their definitions before `evaluate`,
# `diagonal`, or `gradient` is called, so these functions can assume valid
# parameters and input matrices with matching columns.
#
# The registry is internal. Adding a kernel is described in CONTRIBUTING.md.


# Parameter constraints --------------------------------------------------------
#
# A constraint maps a natural parameter theta to the unconstrained coordinate
# eta used by the optimizer. `jacobian` is d theta / d eta, which converts a
# derivative with respect to theta into one with respect to eta.

.new_parameter_constraint <- function(
    type,
    coordinate,
    vector_requirement,
    scalar_requirement,
    contains,
    to_coordinate,
    from_coordinate,
    jacobian,
    coordinate_bounds) {
  structure(
    list(
      type = type,
      coordinate = coordinate,
      vector_requirement = vector_requirement,
      scalar_requirement = scalar_requirement,
      contains = contains,
      to_coordinate = to_coordinate,
      from_coordinate = from_coordinate,
      jacobian = jacobian,
      coordinate_bounds = coordinate_bounds
    ),
    class = "gaussianprocesses_parameter_constraint"
  )
}


# theta > 0 with eta = log(theta). `coordinate_bounds` receives the
# natural-scale bounds of optimize_gp().
.positive_constraint <- function() {
  .new_parameter_constraint(
    type = "positive",
    coordinate = "log",
    vector_requirement = "finite, strictly positive values",
    scalar_requirement = "greater than 0",
    contains = function(theta) is.finite(theta) & theta > 0,
    to_coordinate = function(theta) log(theta),
    from_coordinate = function(eta) exp(eta),
    jacobian = function(theta) theta,
    coordinate_bounds = function(lower, upper) c(log(lower), log(upper))
  )
}


# theta in R with eta = theta, unbounded.
.real_constraint <- function() {
  .new_parameter_constraint(
    type = "real",
    coordinate = "identity",
    vector_requirement = "finite values",
    scalar_requirement = "finite",
    contains = function(theta) is.finite(theta),
    to_coordinate = function(theta) theta,
    from_coordinate = function(eta) eta,
    jacobian = function(theta) rep(1, length(theta)),
    coordinate_bounds = function(lower, upper) c(-Inf, Inf)
  )
}


# Largest |eta| the optimizer may use for a bounded parameter. At |eta| = 30
# the parameter lies about 1e-13 of the interval width inside its bounds, so
# it never rounds onto a bound for intervals of moderate magnitude.
.logit_coordinate_limit <- 30


# lower < theta < upper with eta = logit((theta - lower) / (upper - lower)),
# so d theta / d eta = (theta - lower) (upper - theta) / (upper - lower).
.bounded_constraint <- function(lower, upper) {
  if (!is.numeric(lower) || !is.numeric(upper) ||
      length(lower) != 1L || length(upper) != 1L ||
      !is.finite(lower) || !is.finite(upper) || lower >= upper) {
    .gp_abort(
      "A bounded constraint needs finite bounds with 'lower' < 'upper'.",
      "gaussianprocesses_validation_error"
    )
  }

  width <- upper - lower
  interval <- sprintf(
    "strictly between %s and %s",
    format(lower),
    format(upper)
  )

  .new_parameter_constraint(
    type = "bounded",
    coordinate = "logit",
    vector_requirement = paste("values", interval),
    scalar_requirement = interval,
    contains = function(theta) {
      is.finite(theta) & theta > lower & theta < upper
    },
    to_coordinate = function(theta) stats::qlogis((theta - lower) / width),
    from_coordinate = function(eta) lower + width * stats::plogis(eta),
    jacobian = function(theta) (theta - lower) * (upper - theta) / width,
    coordinate_bounds = function(lower_bound, upper_bound) {
      c(-.logit_coordinate_limit, .logit_coordinate_limit)
    }
  )
}


.kernel_parameter <- function(
    default,
    shape = c("scalar", "ard"),
    constraint = .positive_constraint()) {
  list(
    default = default,
    shape = match.arg(shape),
    constraint = constraint
  )
}


# Registry ---------------------------------------------------------------------

.kernel_registry <- new.env(parent = emptyenv())

.kernel_registry_state <- new.env(parent = emptyenv())


.ensure_builtin_kernels <- function() {
  if (!isTRUE(.kernel_registry_state$builtins_registered)) {
    .register_builtin_kernels()
    .kernel_registry_state$builtins_registered <- TRUE
  }

  invisible(NULL)
}

.composite_kernel_types <- c(
  "sum", "product", "scale", "select", "changepoint"
)


.validate_kernel_entry <- function(entry) {
  fields <- c(
    "type", "label", "parameters", "stationary", "smoothness", "evaluate",
    "diagonal", "gradient", "input_gradient", "input_hessian",
    "input_parameter_gradient", "state_space"
  )
  missing_fields <- setdiff(fields, names(entry))

  if (length(missing_fields) > 0L) {
    .gp_abort(
      sprintf(
        "A kernel registry entry is missing: %s.",
        paste(missing_fields, collapse = ", ")
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (!is.character(entry$type) || length(entry$type) != 1L ||
      !nzchar(entry$type) || entry$type %in% .composite_kernel_types) {
    .gp_abort(
      "A kernel registry entry needs a new, non-empty 'type'.",
      "gaussianprocesses_validation_error"
    )
  }

  if (!is.list(entry$parameters) || length(entry$parameters) == 0L ||
      is.null(names(entry$parameters)) ||
      any(!nzchar(names(entry$parameters))) ||
      anyDuplicated(names(entry$parameters))) {
    .gp_abort(
      sprintf(
        "Kernel '%s' needs uniquely named parameter definitions.",
        entry$type
      ),
      "gaussianprocesses_validation_error"
    )
  }

  for (name in names(entry$parameters)) {
    definition <- entry$parameters[[name]]

    if (!inherits(
      definition$constraint,
      "gaussianprocesses_parameter_constraint"
    ) ||
      !all(definition$constraint$contains(definition$default))) {
      .gp_abort(
        sprintf(
          "Default of parameter '%s' of kernel '%s' violates its constraint.",
          name,
          entry$type
        ),
        "gaussianprocesses_validation_error"
      )
    }
  }

  for (field in c("evaluate", "diagonal", "gradient")) {
    if (!is.function(entry[[field]])) {
      .gp_abort(
        sprintf(
          "Field '%s' of kernel '%s' must be a function.",
          field,
          entry$type
        ),
        "gaussianprocesses_validation_error"
      )
    }
  }

  if (!is.numeric(entry$smoothness) || length(entry$smoothness) != 1L ||
      is.na(entry$smoothness) || entry$smoothness < 0) {
    .gp_abort(
      sprintf("Kernel '%s' needs a non-negative 'smoothness'.", entry$type),
      "gaussianprocesses_validation_error"
    )
  }

  differentiable <- is.function(entry$input_gradient) &&
    is.function(entry$input_hessian) &&
    is.function(entry$input_parameter_gradient)

  if (entry$smoothness >= 1 && !differentiable) {
    .gp_abort(
      sprintf(
        paste0(
          "Kernel '%s' is differentiable, so 'input_gradient', ",
          "'input_hessian', and 'input_parameter_gradient' must be ",
          "functions."
        ),
        entry$type
      ),
      "gaussianprocesses_validation_error"
    )
  }

  for (field in c(
    "input_gradient", "input_hessian", "input_parameter_gradient",
    "state_space"
  )) {
    if (!is.null(entry[[field]]) && !is.function(entry[[field]])) {
      .gp_abort(
        sprintf(
          "Field '%s' of kernel '%s' must be NULL or a function.",
          field,
          entry$type
        ),
        "gaussianprocesses_validation_error"
      )
    }
  }

  invisible(entry)
}


.register_kernel <- function(entry) {
  .validate_kernel_entry(entry)
  assign(entry$type, entry, envir = .kernel_registry)
  invisible(entry)
}


.unregister_kernel <- function(type) {
  if (exists(type, envir = .kernel_registry, inherits = FALSE)) {
    rm(list = type, envir = .kernel_registry)
  }

  invisible(NULL)
}


.registered_kernel_types <- function() {
  .ensure_builtin_kernels()
  sort(ls(.kernel_registry, all.names = TRUE))
}


.is_leaf_kernel_type <- function(type) {
  .ensure_builtin_kernels()
  is.character(type) && length(type) == 1L &&
    exists(type, envir = .kernel_registry, inherits = FALSE)
}


.kernel_entry <- function(type) {
  if (!.is_leaf_kernel_type(type)) {
    .gp_abort(
      sprintf("Unknown kernel type '%s'.", paste(type, collapse = ", ")),
      "gaussianprocesses_validation_error"
    )
  }

  get(type, envir = .kernel_registry, inherits = FALSE)
}


# Parameter validation ---------------------------------------------------------

.validate_kernel_parameter <- function(value, name, definition) {
  constraint <- definition$constraint

  if (identical(definition$shape, "ard")) {
    if (!is.numeric(value) ||
        !is.null(dim(value)) ||
        length(value) == 0L ||
        !all(constraint$contains(value))) {
      .gp_abort(
        sprintf(
          "'%s' must be a non-empty numeric vector containing only %s.",
          name,
          constraint$vector_requirement
        ),
        "gaussianprocesses_validation_error"
      )
    }

    return(invisible(value))
  }

  if (identical(constraint$type, "positive")) {
    .assert_numeric_scalar(value, name, lower = 0, lower_open = TRUE)
  } else {
    .assert_numeric_scalar(value, name)

    if (!constraint$contains(value)) {
      .gp_abort(
        sprintf("'%s' must be %s.", name, constraint$scalar_requirement),
        "gaussianprocesses_validation_error"
      )
    }
  }

  invisible(value)
}


.validate_leaf_kernel_parameters <- function(type, parameters) {
  definitions <- .kernel_entry(type)$parameters

  for (name in names(definitions)) {
    .validate_kernel_parameter(parameters[[name]], name, definitions[[name]])
  }

  invisible(parameters)
}


# ARD parameters hold one value or one value per input dimension.
.validate_kernel_parameter_dimensions <- function(
    type,
    parameters,
    n_features) {
  definitions <- .kernel_entry(type)$parameters

  for (name in names(definitions)) {
    if (identical(definitions[[name]]$shape, "ard") &&
        !(length(parameters[[name]]) %in% c(1L, n_features))) {
      .gp_abort(
        sprintf(
          paste0(
            "'%s' must contain either one value or one value per ",
            "input dimension. Received %d value(s) for %d dimension(s)."
          ),
          name,
          length(parameters[[name]]),
          n_features
        ),
        "gaussianprocesses_validation_error"
      )
    }
  }

  invisible(parameters)
}


# Parameters in registry order, with defaults for parameters not supplied.
.complete_kernel_parameters <- function(type, parameters) {
  definitions <- .kernel_entry(type)$parameters
  unknown <- setdiff(names(parameters), names(definitions))

  if (length(unknown) > 0L) {
    .gp_abort(
      sprintf(
        "Unknown parameter(s) for kernel '%s': %s.",
        type,
        paste(unknown, collapse = ", ")
      ),
      "gaussianprocesses_validation_error"
    )
  }

  completed <- vector("list", length(definitions))
  names(completed) <- names(definitions)

  for (name in names(definitions)) {
    completed[name] <- if (name %in% names(parameters)) {
      parameters[name]
    } else {
      list(definitions[[name]]$default)
    }
  }

  completed
}


# Covariance of one leaf kernel. Parameters are validated on every call, as
# the exported kernel_*() functions always have.
.evaluate_leaf_parameters <- function(type, parameters, x, y = NULL) {
  entry <- .kernel_entry(type)
  .validate_leaf_kernel_parameters(type, parameters)
  inputs <- .prepare_kernel_inputs(x, y)
  .validate_kernel_parameter_dimensions(type, parameters, ncol(inputs$x))

  entry$evaluate(inputs$x, inputs$y, parameters)
}


# Shared building blocks for entries ------------------------------------------

.constant_diagonal <- function(x, parameters) {
  rep(parameters$variance, nrow(x))
}


# One derivative matrix for a scalar (isotropic) length scale, or one per
# input dimension for ARD length scales.
.length_scale_derivatives <- function(factor, components, length_scale) {
  if (length(length_scale) == 1L) {
    return(factor * Reduce(`+`, components))
  }

  lapply(components, function(component) factor * component)
}


# Gradient of a radial kernel k = variance * f(r^2) with
# r^2 = sum_d s_d and s_d = ((x_d - y_d) / l_d)^2. Since
# d s_d / d log(l_d) = -2 s_d, the length-scale derivative is
# dk / d log(l_d) = g * s_d with g = -2 * variance * f'(r^2).
# `length_scale_factor(covariance, squared_distance, parameters)` returns g;
# `extra(covariance, squared_distance, parameters)` returns the derivatives of
# any further parameters.
.radial_gradient <- function(length_scale_factor, extra = NULL) {
  function(x, y, parameters, covariance) {
    components <- .scaled_squared_difference_components(
      x,
      y,
      parameters$length_scale
    )
    squared_distance <- Reduce(`+`, components)

    derivatives <- list(
      variance = covariance,
      length_scale = .length_scale_derivatives(
        length_scale_factor(covariance, squared_distance, parameters),
        components,
        parameters$length_scale
      )
    )

    if (!is.null(extra)) {
      derivatives <- c(
        derivatives,
        extra(covariance, squared_distance, parameters)
      )
    }

    derivatives
  }
}


.proportional_gradient <- function(x, y, parameters, covariance) {
  list(variance = covariance)
}


.radial_kernel_parameters <- function() {
  list(
    variance = .kernel_parameter(1),
    length_scale = .kernel_parameter(1, shape = "ard")
  )
}


# k = variance * exp(-2 * sum_d q_d) with
# q_d = sin^2(pi * (x_d - y_d) / p) / l_d^2.
.periodic_gradient <- function(x, y, parameters, covariance) {
  period <- parameters$period
  resolved <- .resolve_length_scale(parameters$length_scale, ncol(x))

  angles <- lapply(
    seq_len(ncol(x)),
    function(dimension) {
      pi * outer(x[, dimension], y[, dimension], FUN = "-") / period
    }
  )

  # dk / d log(l_d) = 4 k q_d.
  components <- Map(
    function(angle, length_scale) sin(angle)^2 / length_scale^2,
    angles,
    resolved
  )
  length_scale <- if (length(parameters$length_scale) == 1L) {
    4 * covariance * Reduce(`+`, components)
  } else {
    lapply(components, function(component) 4 * covariance * component)
  }

  # dk / d log(p) = 2 k sum_d (angle_d * sin(2 angle_d) / l_d^2).
  period_terms <- Map(
    function(angle, length_scale) angle * sin(2 * angle) / length_scale^2,
    angles,
    resolved
  )

  list(
    variance = covariance,
    length_scale = length_scale,
    period = 2 * covariance * Reduce(`+`, period_terms)
  )
}


# Built-in kernels -------------------------------------------------------------
#
# Covariance formulas and derivations are in vignette("v02-kernels-and-
# composition") and vignette("v07-kernel-derivatives").
#
# The built-in kernels are registered on the first use of the registry
# rather than when the package is built, so the entries' functions are
# created from the loaded namespace; coverage tools can then trace them.
.register_builtin_kernels <- function() {
  .register_kernel(
    list(
      type = "rbf",
      label = "RBF",
      parameters = .radial_kernel_parameters(),
      stationary = TRUE,
      smoothness = Inf,
      evaluate = function(x, y, parameters) {
        scaled_squared_distance <- .scaled_pairwise_squared_distance(
          x,
          y,
          parameters$length_scale
        )

        parameters$variance * exp(-0.5 * scaled_squared_distance)
      },
      diagonal = .constant_diagonal,
      gradient = .radial_gradient(
        function(covariance, squared_distance, parameters) covariance
      ),
      input_gradient = .radial_input_gradient(.rbf_input_coefficients),
      input_hessian = .radial_input_hessian(.rbf_input_coefficients),
      input_parameter_gradient = .radial_input_parameter_gradient(
        .rbf_input_coefficients
      ),
      state_space = NULL
    )
  )

  .register_kernel(
    list(
      type = "matern12",
      label = "Matern-1/2",
      parameters = .radial_kernel_parameters(),
      stationary = TRUE,
      smoothness = 0,
      evaluate = function(x, y, parameters) {
        scaled_distance <- sqrt(
          .scaled_pairwise_squared_distance(x, y, parameters$length_scale)
        )

        parameters$variance * exp(-scaled_distance)
      },
      diagonal = .constant_diagonal,
      gradient = .radial_gradient(
        function(covariance, squared_distance, parameters) {
          distance <- sqrt(squared_distance)
          factor <- covariance / distance
          # g * s_d -> 0 as r -> 0 because s_d <= r^2.
          factor[distance == 0] <- 0
          factor
        }
      ),
      input_gradient = NULL,
      input_hessian = NULL,
      input_parameter_gradient = NULL,
      state_space = .matern_state_space(0L)
    )
  )

  .register_kernel(
    list(
      type = "matern32",
      label = "Matern-3/2",
      parameters = .radial_kernel_parameters(),
      stationary = TRUE,
      smoothness = 1,
      evaluate = function(x, y, parameters) {
        distance <- sqrt(
          .scaled_pairwise_squared_distance(x, y, parameters$length_scale)
        )
        scaled_distance <- sqrt(3) * distance

        parameters$variance * (1 + scaled_distance) * exp(-scaled_distance)
      },
      diagonal = .constant_diagonal,
      gradient = .radial_gradient(
        function(covariance, squared_distance, parameters) {
          3 * covariance / (1 + sqrt(3) * sqrt(squared_distance))
        }
      ),
      input_gradient = .radial_input_gradient(.matern32_input_coefficients),
      input_hessian = .radial_input_hessian(.matern32_input_coefficients),
      input_parameter_gradient = .radial_input_parameter_gradient(
        .matern32_input_coefficients
      ),
      state_space = .matern_state_space(1L)
    )
  )

  .register_kernel(
    list(
      type = "matern52",
      label = "Matern-5/2",
      parameters = .radial_kernel_parameters(),
      stationary = TRUE,
      smoothness = 2,
      evaluate = function(x, y, parameters) {
        distance <- sqrt(
          .scaled_pairwise_squared_distance(x, y, parameters$length_scale)
        )
        scaled_distance <- sqrt(5) * distance

        parameters$variance *
          (1 + scaled_distance + 5 * distance^2 / 3) *
          exp(-scaled_distance)
      },
      diagonal = .constant_diagonal,
      gradient = .radial_gradient(
        function(covariance, squared_distance, parameters) {
          distance <- sqrt(squared_distance)

          (5 / 3) * covariance * (1 + sqrt(5) * distance) /
            (1 + sqrt(5) * distance + 5 * squared_distance / 3)
        }
      ),
      input_gradient = .radial_input_gradient(.matern52_input_coefficients),
      input_hessian = .radial_input_hessian(.matern52_input_coefficients),
      input_parameter_gradient = .radial_input_parameter_gradient(
        .matern52_input_coefficients
      ),
      state_space = .matern_state_space(2L)
    )
  )

  .register_kernel(
    list(
      type = "rational_quadratic",
      label = "RationalQuadratic",
      parameters = c(
        .radial_kernel_parameters(),
        list(alpha = .kernel_parameter(1))
      ),
      stationary = TRUE,
      smoothness = Inf,
      evaluate = function(x, y, parameters) {
        scaled_squared_distance <- .scaled_pairwise_squared_distance(
          x,
          y,
          parameters$length_scale
        )

        alpha <- parameters$alpha

        parameters$variance *
          (1 + scaled_squared_distance / (2 * alpha))^(-alpha)
      },
      diagonal = .constant_diagonal,
      gradient = .radial_gradient(
        function(covariance, squared_distance, parameters) {
          covariance / (1 + squared_distance / (2 * parameters$alpha))
        },
        extra = function(covariance, squared_distance, parameters) {
          alpha <- parameters$alpha
          half_scaled <- squared_distance / (2 * alpha)

          # d log(k) / d log(alpha) = -alpha * log(u) + r^2 / (2u),
          # with u = 1 + r^2 / (2 alpha).
          list(
            alpha = covariance * (
              squared_distance / (2 * (1 + half_scaled)) -
                alpha * log1p(half_scaled)
            )
          )
        }
      ),
      input_gradient = .radial_input_gradient(
        .rational_quadratic_input_coefficients
      ),
      input_hessian = .radial_input_hessian(
        .rational_quadratic_input_coefficients
      ),
      input_parameter_gradient = .radial_input_parameter_gradient(
        .rational_quadratic_input_coefficients
      ),
      state_space = NULL
    )
  )

  .register_kernel(
    list(
      type = "periodic",
      label = "Periodic",
      parameters = c(
        .radial_kernel_parameters(),
        list(period = .kernel_parameter(1))
      ),
      stationary = TRUE,
      smoothness = Inf,
      evaluate = function(x, y, parameters) {
        periodic_distance <- .pairwise_periodic_scaled_distance(
          x,
          y,
          length_scale = parameters$length_scale,
          period = parameters$period
        )

        parameters$variance * exp(-2 * periodic_distance)
      },
      diagonal = .constant_diagonal,
      gradient = .periodic_gradient,
      input_gradient = function(x, y, parameters) {
        .periodic_input_derivatives(x, y, parameters, 1L)
      },
      input_hessian = function(x, y, parameters) {
        .periodic_input_derivatives(x, y, parameters, 2L)
      },
      input_parameter_gradient = .periodic_input_parameter_derivatives,
      state_space = NULL
    )
  )

  .register_kernel(
    list(
      type = "linear",
      label = "Linear",
      parameters = list(variance = .kernel_parameter(1)),
      stationary = FALSE,
      smoothness = Inf,
      evaluate = function(x, y, parameters) {
        parameters$variance * tcrossprod(x, y)
      },
      diagonal = function(x, parameters) {
        parameters$variance * rowSums(x^2)
      },
      gradient = .proportional_gradient,
      input_gradient = function(x, y, parameters) {
        .linear_input_derivatives(x, y, parameters, 1L)
      },
      input_hessian = function(x, y, parameters) {
        .linear_input_derivatives(x, y, parameters, 2L)
      },
      input_parameter_gradient = function(x, y, parameters) {
        .linear_input_derivatives(x, y, parameters, 3L)
      },
      state_space = NULL
    )
  )

  .register_kernel(
    list(
      type = "white_noise",
      label = "WhiteNoise",
      parameters = list(variance = .kernel_parameter(1)),
      stationary = TRUE,
      smoothness = 0,
      evaluate = function(x, y, parameters) {
        parameters$variance * .pairwise_identical_rows(x, y)
      },
      diagonal = .constant_diagonal,
      gradient = .proportional_gradient,
      input_gradient = NULL,
      input_hessian = NULL,
      input_parameter_gradient = NULL,
      state_space = NULL
    )
  )

  .register_kernel(.spectral_component_entry())

  invisible(NULL)
}
