# The kernel registry: built-in entries, parameter constraints, and a kernel
# registered only inside this file to show that one registry entry is enough
# for a kernel to work throughout the package.

built_in_constructors <- list(
  rbf = rbf_kernel,
  matern12 = matern12_kernel,
  matern32 = matern32_kernel,
  matern52 = matern52_kernel,
  rational_quadratic = rational_quadratic_kernel,
  periodic = periodic_kernel,
  linear = linear_kernel,
  white_noise = white_noise_kernel,
  # spectral_mixture_kernel() combines components of this type.
  spectral_component = gaussianprocesses:::.spectral_component_kernel
)

# Kernels with a covariance-matrix function kernel_name().
built_in_matrix_functions <- list(
  rbf = kernel_rbf,
  matern12 = kernel_matern12,
  matern32 = kernel_matern32,
  matern52 = kernel_matern52,
  rational_quadratic = kernel_rational_quadratic,
  periodic = kernel_periodic,
  linear = kernel_linear,
  white_noise = kernel_white_noise
)

registry_x <- cbind(
  c(-1.2, -0.4, 0, 0, 0.35, 1.1, 2),
  c(0.5, -0.3, 1.2, 1.2, 0, -0.8, 0.4),
  c(2, 0.1, -0.7, -0.7, 1.4, 0.3, -1)
)
registry_y <- registry_x[c(2, 5, 7), ] + 0.3

# Number of mean-square derivatives of each built-in process.
expected_smoothness <- list(
  rbf = Inf,
  matern12 = 0,
  matern32 = 1,
  matern52 = 2,
  rational_quadratic = Inf,
  periodic = Inf,
  linear = Inf,
  white_noise = 0,
  spectral_component = Inf
)


# What every registered kernel must satisfy: a symmetric positive
# semidefinite Gram matrix, a diagonal that matches it, and analytical
# derivatives that match central finite differences, including derivatives
# with respect to the inputs when the kernel is differentiable.
expect_kernel_contract <- function(kernel, x, y, label) {
  if (gaussianprocesses:::.kernel_smoothness(kernel) >= 1) {
    expect_input_derivative_contract(kernel, x, y, label)
    expect_input_parameter_gradient_contract(kernel, x, y, label)
  } else {
    expect_error(
      kernel_input_gradient(kernel, x, y),
      class = "gaussianprocesses_smoothness_error"
    )
  }

  covariance <- evaluate_kernel(kernel, x)
  scale <- max(abs(diag(covariance)))

  expect_equal(covariance, t(covariance), tolerance = 1e-14, info = label)
  expect_gte(
    min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values),
    -1e-10 * scale,
    label = paste(label, "smallest eigenvalue")
  )
  expect_equal(
    kernel_diagonal(kernel, x),
    diag(covariance),
    tolerance = 1e-14,
    info = label
  )

  for (inputs in list(list(x = x, y = NULL), list(x = x, y = y))) {
    analytical <- kernel_gradient(kernel, inputs$x, inputs$y)
    numerical <- coordinate_fd_kernel_gradient(kernel, inputs$x, inputs$y)

    expect_identical(names(analytical), names(numerical), info = label)

    for (path in names(analytical)) {
      expect_lte(
        max(abs(analytical[[path]] - numerical[[path]])) /
          max(1, max(abs(analytical[[path]]))),
        1e-7,
        label = paste(label, path)
      )
    }
  }
}


# A kernel with non-default values of every parameter: each coordinate is
# shifted from the default, with distinct ARD values in three dimensions.
kernel_with_shifted_parameters <- function(type) {
  definitions <- gaussianprocesses:::.kernel_entry(type)$parameters
  parameters <- lapply(
    definitions,
    function(definition) {
      default <- definition$default
      offset <- 0.2

      if (identical(definition$shape, "ard")) {
        default <- rep(default, ncol(registry_x))
        offset <- c(-0.3, 0.2, 0.6)
      }

      definition$constraint$from_coordinate(
        definition$constraint$to_coordinate(default) + offset
      )
    }
  )

  gaussianprocesses:::.new_leaf_kernel(type, parameters)
}


test_that("the registry defines exactly the built-in leaf kernels", {
  expect_identical(
    gaussianprocesses:::.registered_kernel_types(),
    sort(names(built_in_constructors))
  )

  for (type in names(built_in_constructors)) {
    entry <- gaussianprocesses:::.kernel_entry(type)

    expect_identical(entry$type, type)
    expect_true(is.character(entry$label) && nzchar(entry$label))
    expect_identical(entry$stationary, !identical(type, "linear"), info = type)
    expect_identical(
      entry$smoothness,
      expected_smoothness[[type]],
      info = type
    )
    differentiable <- entry$smoothness >= 1
    expect_identical(is.function(entry$input_gradient), differentiable)
    expect_identical(is.function(entry$input_hessian), differentiable)
    expect_identical(
      is.function(entry$input_parameter_gradient),
      differentiable
    )
    expect_identical(
      is.function(entry$state_space),
      type %in% c("matern12", "matern32", "matern52"),
      info = type
    )
    expect_true(
      all(names(entry) %in% c(
        "type", "label", "parameters", "stationary", "smoothness", "evaluate",
        "diagonal", "gradient", "input_gradient", "input_hessian",
        "input_parameter_gradient", "state_space"
      ))
    )
  }
})


test_that("registry entries declare a valid smoothness", {
  entry <- gaussianprocesses:::.kernel_entry("rbf")

  for (smoothness in list(NULL, -1, NA_real_, c(1, 2), "Inf")) {
    invalid <- entry
    invalid$type <- "invalid_smoothness"
    invalid["smoothness"] <- list(smoothness)
    expect_error(
      gaussianprocesses:::.register_kernel(invalid),
      class = "gaussianprocesses_validation_error"
    )
  }

  missing_derivatives <- entry
  missing_derivatives$type <- "missing_derivatives"
  missing_derivatives["input_hessian"] <- list(NULL)
  expect_error(
    gaussianprocesses:::.register_kernel(missing_derivatives),
    "must be functions",
    class = "gaussianprocesses_validation_error"
  )
  expect_false(
    "missing_derivatives" %in% gaussianprocesses:::.registered_kernel_types()
  )
})


test_that("constructor and matrix-function defaults match the registry", {
  for (type in names(built_in_constructors)) {
    definitions <- gaussianprocesses:::.kernel_entry(type)$parameters
    defaults <- lapply(definitions, `[[`, "default")

    constructor_defaults <- lapply(
      formals(built_in_constructors[[type]]),
      eval
    )
    expect_identical(constructor_defaults, defaults, info = type)

    if (!is.null(built_in_matrix_functions[[type]])) {
      matrix_defaults <- lapply(
        formals(built_in_matrix_functions[[type]])[names(defaults)],
        eval
      )
      expect_identical(matrix_defaults, defaults, info = type)
    }

    expect_identical(
      built_in_constructors[[type]](),
      gaussianprocesses:::.new_leaf_kernel(type),
      info = type
    )
  }
})


test_that("every built-in kernel satisfies the kernel contract", {
  for (type in names(built_in_constructors)) {
    expect_kernel_contract(
      gaussianprocesses:::.new_leaf_kernel(type),
      registry_x,
      registry_y,
      paste(type, "defaults")
    )
    expect_kernel_contract(
      kernel_with_shifted_parameters(type),
      registry_x,
      registry_y,
      paste(type, "shifted")
    )
  }
})


test_that("parameter constraints transform coordinates consistently", {
  constraints <- list(
    positive = list(
      constraint = gaussianprocesses:::.positive_constraint(),
      values = c(1e-3, 0.4, 1, 7.5)
    ),
    real = list(
      constraint = gaussianprocesses:::.real_constraint(),
      values = c(-3, -0.2, 0, 1.6)
    ),
    bounded = list(
      constraint = gaussianprocesses:::.bounded_constraint(-1, 3),
      values = c(-0.99, -0.5, 1, 2.7)
    )
  )
  step <- 1e-6

  for (name in names(constraints)) {
    constraint <- constraints[[name]]$constraint
    values <- constraints[[name]]$values
    coordinates <- constraint$to_coordinate(values)

    expect_true(all(constraint$contains(values)), info = name)
    expect_equal(
      constraint$from_coordinate(coordinates),
      values,
      tolerance = 1e-14
    )

    # d theta / d eta against central differences of the inverse transform.
    numerical <- (
      constraint$from_coordinate(coordinates + step) -
        constraint$from_coordinate(coordinates - step)
    ) / (2 * step)
    expect_equal(
      constraint$jacobian(values),
      numerical,
      tolerance = 1e-8,
      info = name
    )
  }

  positive <- constraints$positive$constraint
  expect_identical(positive$to_coordinate(2), log(2))
  expect_identical(positive$jacobian(2), 2)
  expect_false(any(positive$contains(c(0, -1, Inf, NA))))

  real <- constraints$real$constraint
  expect_identical(real$to_coordinate(-2.5), -2.5)
  expect_identical(real$jacobian(c(-2, 5)), c(1, 1))
  expect_identical(real$coordinate_bounds(1e-8, 1e8), c(-Inf, Inf))

  bounded <- constraints$bounded$constraint
  theta <- c(-0.5, 1, 2.7)
  expect_equal(bounded$jacobian(theta), (theta + 1) * (3 - theta) / 4)
  expect_false(any(bounded$contains(c(-1, 3, -2, NA))))
  expect_true(all(bounded$contains(bounded$from_coordinate(c(-30, 30)))))

  expect_error(
    gaussianprocesses:::.bounded_constraint(1, 1),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("the registry rejects malformed entries and unknown types", {
  entry <- gaussianprocesses:::.kernel_entry("linear")
  register <- gaussianprocesses:::.register_kernel
  invalid_default <- list(
    variance = gaussianprocesses:::.kernel_parameter(-1)
  )

  expect_error(
    register(entry[setdiff(names(entry), "gradient")]),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    register(utils::modifyList(entry, list(type = "sum"))),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    register(utils::modifyList(entry, list(parameters = invalid_default))),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    register(utils::modifyList(entry, list(evaluate = 1))),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    gaussianprocesses:::.new_leaf_kernel("not_a_kernel"),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    gaussianprocesses:::.new_leaf_kernel("rbf", list(period = 2)),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    gaussianprocesses:::.kernel_entry(NULL),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("built-in kernel parameters are optimized on the log scale", {
  x <- seq(-2, 2, length.out = 12)
  model <- optimize_gp(
    x,
    sin(2 * x),
    kernel = sum_kernel(rbf_kernel(), scale_kernel(linear_kernel(), 2)),
    noise_variance = 0.05,
    n_starts = 1
  )

  expect_identical(model$optimization$parameterization, "log")
  expect_identical(
    model$optimization$coordinates,
    c(
      "kernel1.variance" = "log",
      "kernel1.length_scale" = "log",
      "kernel2.scale" = "log",
      "kernel2.kernel.variance" = "log",
      noise_variance = "log"
    )
  )
})


# A test-only kernel ----------------------------------------------------------
#
# k(x, y) = variance * (weight * sum_d (x_d - center) (y_d - center)
#                       + (1 - weight) * exp(-r^2 / 2)),
# with r^2 the squared distance scaled by ARD length scales. It is positive
# semidefinite as a convex combination of an offset linear kernel and an RBF
# kernel. Its parameters cover every constraint: positive variance and length
# scales (log coordinate), a real-valued center (identity), and a weight in
# (0, 1) (logit). With u_d = (x_d - y_d) / l_d^2 and s = exp(-r^2 / 2), its
# input derivatives are
#   dk / dx_d       = variance (weight (y_d - center) - (1 - weight) s u_d),
#   d2k / dx_d dy_e = variance (weight delta_de
#                     + (1 - weight) s (delta_de / l_d^2 - u_d u_e)).
# Their hyperparameter derivatives use ds / d log l_c = s s_c with
# s_c = (x_c - y_c)^2 / l_c^2, du_d / d log l_c = -2 u_d delta_dc, and the
# logit Jacobian weight (1 - weight) of the weight.

test_kernel_type <- "test_offset_mixture"

test_kernel_entry <- function() {
  weight_constraint <- gaussianprocesses:::.bounded_constraint(0, 1)

  scaled_components <- function(x, y, length_scale) {
    length_scale <- rep_len(length_scale, ncol(x))
    lapply(
      seq_len(ncol(x)),
      function(d) {
        outer(x[, d] / length_scale[d], y[, d] / length_scale[d], "-")^2
      }
    )
  }

  # s and u_d of the input derivatives.
  input_parts <- function(x, y, parameters) {
    length_scale <- rep_len(parameters$length_scale, ncol(x))

    list(
      length_scale = length_scale,
      smooth = exp(
        -0.5 * Reduce(`+`, scaled_components(x, y, length_scale))
      ),
      scaled = lapply(
        seq_len(ncol(x)),
        function(d) outer(x[, d], y[, d], "-") / length_scale[d]^2
      )
    )
  }

  list(
    type = test_kernel_type,
    label = "TestOffsetMixture",
    parameters = list(
      variance = gaussianprocesses:::.kernel_parameter(1),
      length_scale = gaussianprocesses:::.kernel_parameter(1, shape = "ard"),
      center = gaussianprocesses:::.kernel_parameter(
        0,
        constraint = gaussianprocesses:::.real_constraint()
      ),
      weight = gaussianprocesses:::.kernel_parameter(
        0.5,
        constraint = weight_constraint
      )
    ),
    stationary = FALSE,
    smoothness = Inf,
    evaluate = function(x, y, parameters) {
      linear <- tcrossprod(x - parameters$center, y - parameters$center)
      smooth <- exp(
        -0.5 * Reduce(`+`, scaled_components(x, y, parameters$length_scale))
      )

      parameters$variance *
        (parameters$weight * linear + (1 - parameters$weight) * smooth)
    },
    diagonal = function(x, parameters) {
      parameters$variance * (
        parameters$weight * rowSums((x - parameters$center)^2) +
          1 - parameters$weight
      )
    },
    gradient = function(x, y, parameters, covariance) {
      variance <- parameters$variance
      weight <- parameters$weight
      center <- parameters$center
      components <- scaled_components(x, y, parameters$length_scale)
      linear <- tcrossprod(x - center, y - center)
      smooth <- exp(-0.5 * Reduce(`+`, components))

      # d smooth / d log(l_d) = smooth * s_d.
      length_scale <- lapply(
        components,
        function(component) variance * (1 - weight) * smooth * component
      )
      if (length(parameters$length_scale) == 1L) {
        length_scale <- Reduce(`+`, length_scale)
      }

      list(
        variance = covariance,
        length_scale = length_scale,
        # Identity coordinate: the natural derivative.
        center = -variance * weight *
          outer(rowSums(x - center), rowSums(y - center), "+"),
        # Logit coordinate: the natural derivative times d weight / d eta.
        weight = variance * (linear - smooth) *
          weight_constraint$jacobian(weight)
      )
    },
    input_gradient = function(x, y, parameters) {
      parts <- input_parts(x, y, parameters)

      lapply(
        seq_len(ncol(x)),
        function(d) {
          parameters$variance * (
            parameters$weight *
              outer(rep(1, nrow(x)), y[, d] - parameters$center) -
              (1 - parameters$weight) * parts$smooth * parts$scaled[[d]]
          )
        }
      )
    },
    input_hessian = function(x, y, parameters) {
      parts <- input_parts(x, y, parameters)
      dimensions <- seq_len(ncol(x))

      lapply(
        dimensions,
        function(d) {
          lapply(
            dimensions,
            function(e) {
              smooth_part <- -parts$scaled[[d]] * parts$scaled[[e]]

              if (d == e) {
                smooth_part <- smooth_part + 1 / parts$length_scale[[d]]^2
              }

              parameters$variance * (
                parameters$weight * (d == e) +
                  (1 - parameters$weight) * parts$smooth * smooth_part
              )
            }
          )
        }
      )
    },
    input_parameter_gradient = function(x, y, parameters) {
      parts <- input_parts(x, y, parameters)
      dimensions <- seq_len(ncol(x))
      variance <- parameters$variance
      weight <- parameters$weight
      smooth <- parts$smooth
      u <- parts$scaled
      inverse_squared <- 1 / parts$length_scale^2
      components <- scaled_components(x, y, parts$length_scale)
      offset <- lapply(
        dimensions,
        function(d) outer(rep(1, nrow(x)), y[, d] - parameters$center)
      )
      curvature <- function(d, e) {
        (d == e) * inverse_squared[[d]] - u[[d]] * u[[e]]
      }
      terms <- function(gradient_x, cross) {
        list(
          gradient_x = lapply(dimensions, gradient_x),
          cross = lapply(
            dimensions,
            function(d) lapply(dimensions, function(e) cross(d, e))
          )
        )
      }

      length_scale <- lapply(
        dimensions,
        function(c_index) {
          s_c <- components[[c_index]]

          terms(
            function(d) {
              -variance * (1 - weight) * smooth * u[[d]] *
                (s_c - 2 * (d == c_index))
            },
            function(d, e) {
              variance * (1 - weight) * smooth * (
                s_c * curvature(d, e) -
                  2 * (d == e) * (d == c_index) * inverse_squared[[d]] +
                  2 * ((d == c_index) + (e == c_index)) * u[[d]] * u[[e]]
              )
            }
          )
        }
      )

      if (length(parameters$length_scale) == 1L) {
        length_scale <- Reduce(
          function(first, second) {
            list(
              gradient_x = Map(`+`, first$gradient_x, second$gradient_x),
              cross = Map(
                function(a, b) Map(`+`, a, b),
                first$cross,
                second$cross
              )
            )
          },
          length_scale
        )
      }

      jacobian <- weight_constraint$jacobian(weight)

      list(
        variance = terms(
          function(d) {
            variance * (weight * offset[[d]] - (1 - weight) * smooth * u[[d]])
          },
          function(d, e) {
            variance * (weight * (d == e) + (1 - weight) * smooth *
              curvature(d, e))
          }
        ),
        length_scale = length_scale,
        center = terms(
          function(d) matrix(-variance * weight, nrow(x), nrow(y)),
          function(d, e) matrix(0, nrow(x), nrow(y))
        ),
        weight = terms(
          function(d) variance * jacobian * (offset[[d]] + smooth * u[[d]]),
          function(d, e) {
            variance * jacobian * ((d == e) - smooth * curvature(d, e))
          }
        )
      )
    },
    state_space = NULL
  )
}

test_kernel <- function(...) {
  gaussianprocesses:::.new_leaf_kernel(test_kernel_type, list(...))
}

# Evaluates `code` with the test kernel registered, and always removes it.
with_test_kernel <- function(code) {
  gaussianprocesses:::.register_kernel(test_kernel_entry())
  on.exit(gaussianprocesses:::.unregister_kernel(test_kernel_type), add = TRUE)
  force(code)
}


test_that("a registry entry alone makes a kernel work in kernel functions", {
  with_test_kernel({
    defaults <- test_kernel()
    expect_identical(
      defaults$parameters,
      list(variance = 1, length_scale = 1, center = 0, weight = 0.5)
    )

    kernel <- test_kernel(
      variance = 1.2,
      length_scale = c(0.7, 1.4, 0.9),
      center = -0.5,
      weight = 0.3
    )

    expect_kernel_contract(kernel, registry_x, registry_y, "test kernel")
    expect_kernel_contract(
      test_kernel(length_scale = 0.8, center = 1.2, weight = 0.85),
      registry_x,
      registry_y,
      "test kernel, shared length scale"
    )

    expect_identical(
      names(kernel_parameters(kernel, flatten = TRUE)),
      c("variance", "length_scale[1]", "length_scale[2]", "length_scale[3]",
        "center", "weight")
    )

    updated <- update_kernel_parameters(kernel, c(center = 0.25, weight = 0.9))
    expect_identical(updated$parameters$center, 0.25)
    expect_identical(updated$parameters$weight, 0.9)

    expect_error(
      update_kernel_parameters(kernel, c(weight = 1)),
      "'weight' must be strictly between 0 and 1",
      class = "gaussianprocesses_validation_error"
    )
    expect_error(
      test_kernel(center = Inf),
      class = "gaussianprocesses_validation_error"
    )
    expect_error(
      test_kernel(length_scale = c(1, -2)),
      "containing only finite, strictly positive values",
      class = "gaussianprocesses_validation_error"
    )
    expect_error(
      evaluate_kernel(kernel, registry_x[, 1:2]),
      "Received 3 value\\(s\\) for 2 dimension\\(s\\)",
      class = "gaussianprocesses_validation_error"
    )

    expect_identical(
      utils::capture.output(print(kernel)),
      paste0(
        "TestOffsetMixture(variance=1.2, length_scale=c(0.7, 1.4, 0.9), ",
        "center=-0.5, weight=0.3)"
      )
    )

    composite <- scale_kernel(sum_kernel(kernel, white_noise_kernel(0.1)), 2)
    constraints <- gaussianprocesses:::.flatten_kernel_constraints(composite)
    expect_identical(
      names(constraints),
      names(kernel_parameters(composite, flatten = TRUE))
    )
    expect_identical(
      vapply(constraints, `[[`, character(1), "coordinate")[
        c("scale", "kernel.kernel1.center", "kernel.kernel1.weight")
      ],
      c(
        scale = "log",
        kernel.kernel1.center = "identity",
        kernel.kernel1.weight = "logit"
      )
    )
    expect_kernel_contract(
      composite,
      registry_x,
      registry_y,
      "test kernel composite"
    )
  })
})


test_that("models fit and optimize a registry-only kernel", {
  with_test_kernel({
    x <- cbind(seq(-1, 2, length.out = 25), sin(1:25))
    truth <- test_kernel(
      length_scale = c(0.8, 1.5),
      center = 0.8,
      weight = 0.6
    )
    y <- simulate_gp_data(x, truth, noise_variance = 0.01, seed = 11)$observed
    start <- test_kernel(length_scale = c(1, 1), center = -0.5, weight = 0.3)

    # The marginal-likelihood gradient in each coordinate, including the
    # identity coordinate of the center and the logit coordinate of the weight.
    model <- fit_gp(x, y, start, noise_variance = 0.05)
    gradient <- log_marginal_likelihood_gradient(model)
    constraints <- gaussianprocesses:::.flatten_kernel_constraints(start)
    parameters <- kernel_parameters(start, flatten = TRUE)
    step <- 1e-5

    for (path in names(parameters)) {
      coordinate <- constraints[[path]]$to_coordinate(parameters[[path]])
      lml_at <- function(direction) {
        shifted <- parameters[path]
        shifted[] <- constraints[[path]]$from_coordinate(
          coordinate + direction * step
        )
        kernel <- update_kernel_parameters(start, shifted)
        log_marginal_likelihood(fit_gp(x, y, kernel, noise_variance = 0.05))
      }
      numerical <- (lml_at(1) - lml_at(-1)) / (2 * step)

      expect_equal(gradient[[path]], numerical, tolerance = 1e-6, info = path)
    }

    optimized <- optimize_gp(
      x,
      y,
      kernel = start,
      noise_variance = 0.05,
      n_starts = 2
    )
    optimization <- optimized$optimization

    expect_identical(optimization$parameterization, "mixed")
    expect_identical(
      optimization$coordinates[
        c("variance", "length_scale[1]", "center", "weight", "noise_variance")
      ],
      c(
        variance = "log",
        "length_scale[1]" = "log",
        center = "identity",
        weight = "logit",
        noise_variance = "log"
      )
    )

    # The center starts at its own (negative) value: an identity coordinate.
    expect_identical(optimization$starts[1, "center"], -0.5)
    expect_identical(optimization$starts[1, "weight"], stats::qlogis(0.3))
    expect_true(optimization$converged)
    expect_gt(
      log_marginal_likelihood(optimized),
      log_marginal_likelihood(model)
    )

    # The optimum is a local maximum along the center.
    center <- kernel_parameters(optimized$kernel, flatten = TRUE)[["center"]]
    for (offset in c(-0.05, 0.05)) {
      perturbed <- fit_gp(
        x,
        y,
        update_kernel_parameters(optimized$kernel, c(center = center + offset)),
        noise_variance = optimized$noise_variance
      )
      expect_lt(
        log_marginal_likelihood(perturbed),
        log_marginal_likelihood(optimized)
      )
    }

    # Prediction and FITC both use the registry's evaluation and diagonal. With
    # inducing points at the data, FITC reproduces the exact posterior.
    x_new <- cbind(c(-0.5, 0.4, 1.7), c(0.2, -0.6, 0.9))
    exact <- predict_gp(model, x_new)
    sparse <- fit_sparse_gp(
      x,
      y,
      start,
      noise_variance = 0.05,
      inducing_points = x
    )
    sparse_prediction <- predict_sparse_gp(sparse, x_new)

    expect_equal(sparse_prediction$mean, exact$mean, tolerance = 1e-6)
    expect_equal(
      sparse_prediction$latent_variance,
      exact$latent_variance,
      tolerance = 1e-6
    )

    reduced <- fit_sparse_gp(x, y, start, noise_variance = 0.05, n_inducing = 8)
    expect_true(all(is.finite(predict_sparse_gp(reduced, x_new)$mean)))
  })
})


test_that("the test kernel is removed from the registry after its tests", {
  expect_false(gaussianprocesses:::.is_leaf_kernel_type(test_kernel_type))
})
