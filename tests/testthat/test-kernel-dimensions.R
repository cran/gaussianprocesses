# Kernels restricted to selected input columns.

dimension_x <- cbind(
  c(-1.2, -0.4, 0, 0, 0.35, 1.1, 2),
  c(0.5, -0.3, 1.2, 1.2, 0, -0.8, 0.4),
  c(2, 0.1, -0.7, -0.7, 1.4, 0.3, -1),
  c(0.3, 0.9, -1.5, -1.5, 0.2, 0.6, -0.1)
)
dimension_y <- dimension_x[c(2, 5, 7), ] + 0.3

# Kernels to restrict, with ARD parameters sized for `n_columns` columns.
restrictable_kernels <- function(n_columns) {
  ard <- function(values) values[seq_len(n_columns)]

  list(
    rbf = rbf_kernel(1.3, ard(c(0.7, 1.4, 2.1))),
    matern12 = matern12_kernel(0.8, ard(c(1.1, 0.6, 1.7))),
    matern32 = matern32_kernel(1.1, ard(c(0.9, 1.6, 0.5))),
    matern52 = matern52_kernel(0.9, 1.2),
    rational_quadratic = rational_quadratic_kernel(
      0.7,
      ard(c(1.5, 0.8, 1)),
      alpha = 0.6
    ),
    periodic = periodic_kernel(1.2, ard(c(0.8, 1.3, 1.9)), period = 1.7),
    linear = linear_kernel(0.4),
    white_noise = white_noise_kernel(0.2),
    sum = sum_kernel(
      rbf_kernel(0.6, ard(c(1, 2, 3))),
      linear_kernel(0.3)
    ),
    product = product_kernel(
      periodic_kernel(period = 2.3),
      matern52_kernel(length_scale = ard(c(1.4, 0.7, 1)))
    ),
    scale = scale_kernel(rational_quadratic_kernel(alpha = 2), 1.8)
  )
}

column_selections <- list(1L, c(3L, 1L), c(2L, 4L, 1L))


test_that("a restricted kernel equals the kernel on the selected columns", {
  for (columns in column_selections) {
    kernels <- restrictable_kernels(length(columns))
    x_selected <- dimension_x[, columns, drop = FALSE]
    y_selected <- dimension_y[, columns, drop = FALSE]

    for (name in names(kernels)) {
      kernel <- kernels[[name]]
      restricted <- select_dimensions(kernel, columns)
      label <- sprintf("%s, columns %s", name, toString(columns))

      expect_identical(
        evaluate_kernel(restricted, dimension_x),
        evaluate_kernel(kernel, x_selected),
        info = label
      )
      expect_identical(
        evaluate_kernel(restricted, dimension_x, dimension_y),
        evaluate_kernel(kernel, x_selected, y_selected),
        info = label
      )
      expect_identical(
        kernel_diagonal(restricted, dimension_x),
        kernel_diagonal(kernel, x_selected),
        info = label
      )
      expect_identical(
        kernel_gradient(restricted, dimension_x, dimension_y),
        kernel_gradient(kernel, x_selected, y_selected),
        info = label
      )
    }
  }
})


test_that("selections nest and compose with other kernels", {
  inner <- matern32_kernel(length_scale = c(0.8, 1.5))

  # Inner columns index the outer selection: columns 4 and 2 of the inputs.
  nested <- select_dimensions(select_dimensions(inner, c(3, 1)), c(2, 3, 4))
  expect_identical(
    evaluate_kernel(nested, dimension_x),
    evaluate_kernel(inner, dimension_x[, c(4, 2)])
  )

  additive <- sum_kernel(
    select_dimensions(rbf_kernel(length_scale = 0.6), 1),
    select_dimensions(inner, c(2, 3)),
    scale_kernel(select_dimensions(linear_kernel(), 4), 0.5)
  )
  expect_equal(
    evaluate_kernel(additive, dimension_x, dimension_y),
    evaluate_kernel(rbf_kernel(length_scale = 0.6), dimension_x[, 1],
                    dimension_y[, 1]) +
      evaluate_kernel(inner, dimension_x[, 2:3], dimension_y[, 2:3]) +
      0.5 * evaluate_kernel(linear_kernel(), dimension_x[, 4],
                            dimension_y[, 4]),
    tolerance = 1e-15
  )

  covariance <- evaluate_kernel(additive, dimension_x)
  expect_gte(
    min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values),
    -1e-10 * max(diag(covariance))
  )
})


test_that("selections add no parameter path level", {
  kernel <- rbf_kernel(length_scale = c(0.5, 2))
  restricted <- select_dimensions(kernel, c(1, 3))

  expect_identical(
    kernel_parameters(restricted, flatten = TRUE),
    kernel_parameters(kernel, flatten = TRUE)
  )
  expect_identical(kernel_parameters(restricted), kernel_parameters(kernel))

  additive <- sum_kernel(
    select_dimensions(rbf_kernel(), 1),
    scale_kernel(select_dimensions(periodic_kernel(), 2), 2)
  )
  expect_identical(
    names(kernel_parameters(additive, flatten = TRUE)),
    c(
      "kernel1.variance", "kernel1.length_scale",
      "kernel2.scale", "kernel2.kernel.variance",
      "kernel2.kernel.length_scale", "kernel2.kernel.period"
    )
  )

  updated <- update_kernel_parameters(
    additive,
    c("kernel1.length_scale" = 0.3, "kernel2.kernel.period" = 4)
  )
  expect_s3_class(updated, "gaussianprocesses_kernel")
  expect_identical(updated$children[[1]]$columns, 1L)
  expect_identical(
    updated$children[[1]]$children[[1]]$parameters$length_scale,
    0.3
  )
  expect_equal(
    evaluate_kernel(updated, dimension_x),
    evaluate_kernel(rbf_kernel(length_scale = 0.3), dimension_x[, 1]) +
      2 * evaluate_kernel(periodic_kernel(period = 4), dimension_x[, 2]),
    tolerance = 1e-15
  )
  expect_identical(
    names(gaussianprocesses:::.flatten_kernel_constraints(additive)),
    names(kernel_parameters(additive, flatten = TRUE))
  )
})


test_that("hyperparameter gradients of restricted kernels are correct", {
  kernels <- list(
    single = select_dimensions(matern52_kernel(1.2, c(0.7, 1.9)), c(4, 2)),
    additive = sum_kernel(
      select_dimensions(rbf_kernel(0.8, 0.6), 1),
      select_dimensions(
        product_kernel(
          periodic_kernel(1.1, c(0.9, 1.4), period = 1.6),
          rational_quadratic_kernel(0.5, 1.3, alpha = 0.9)
        ),
        c(2, 3)
      ),
      scale_kernel(select_dimensions(linear_kernel(0.3), 4), 1.5)
    )
  )

  for (name in names(kernels)) {
    for (y in list(NULL, dimension_y)) {
      analytical <- kernel_gradient(kernels[[name]], dimension_x, y)
      numerical <- finite_difference_kernel_gradient(
        kernels[[name]],
        dimension_x,
        y
      )

      expect_identical(names(analytical), names(numerical), info = name)

      for (path in names(analytical)) {
        expect_lte(
          scaled_gradient_discrepancy(analytical[[path]], numerical[[path]]),
          gradient_fd_tolerance,
          label = paste0(name, ": ", path)
        )
      }
    }
  }
})


test_that("input derivatives along unselected columns are zero", {
  inner <- product_kernel(
    matern52_kernel(length_scale = c(0.9, 1.6)),
    linear_kernel(0.5)
  )
  restricted <- select_dimensions(inner, c(3, 1))
  derivatives <- gaussianprocesses:::.kernel_input_derivatives(
    restricted,
    dimension_x,
    dimension_y,
    order = 2L
  )
  selected <- gaussianprocesses:::.kernel_input_derivatives(
    inner,
    dimension_x[, c(3, 1)],
    dimension_y[, c(3, 1)],
    order = 2L
  )
  # Position of each input column in the selection, NA when unselected.
  positions <- c(2L, NA, 1L, NA)
  zero <- matrix(0, nrow(dimension_x), nrow(dimension_y))

  for (d in 1:4) {
    for (side in c("gradient_x", "gradient_y")) {
      expected <- zero

      if (!is.na(positions[d])) {
        expected <- selected[[side]][[positions[d]]]
      }

      expect_identical(derivatives[[side]][[d]], expected, info = c(side, d))
    }

    for (e in 1:4) {
      expected <- zero

      if (!is.na(positions[d]) && !is.na(positions[e])) {
        expected <- selected$cross[[positions[d]]][[positions[e]]]
      }

      expect_identical(derivatives$cross[[d]][[e]], expected, info = c(d, e))
    }
  }

  additive <- sum_kernel(
    select_dimensions(rbf_kernel(length_scale = 0.7), 2),
    select_dimensions(
      periodic_kernel(period = 2, length_scale = c(1, 1.5)),
      c(1, 4)
    )
  )

  for (kernel in list(restricted, additive)) {
    expect_input_derivative_contract(
      kernel,
      dimension_x,
      dimension_y,
      "restricted kernel"
    )
  }

  expect_error(
    kernel_input_gradient(select_dimensions(matern12_kernel(), 2), dimension_x),
    class = "gaussianprocesses_smoothness_error"
  )
})


test_that("invalid selections raise classed errors", {
  for (columns in list(0, -1, 1.5, NA, Inf, c(1, 1), "1", numeric(0),
                       matrix(1:2, 1), TRUE)) {
    expect_error(
      select_dimensions(rbf_kernel(), columns),
      class = "gaussianprocesses_validation_error"
    )
  }

  expect_error(
    select_dimensions(list(), 1),
    class = "gaussianprocesses_validation_error"
  )
  expect_identical(select_dimensions(rbf_kernel(), c(2, 1))$columns, c(2L, 1L))

  out_of_range <- select_dimensions(rbf_kernel(), c(2, 5))
  for (call in list(
    function() evaluate_kernel(out_of_range, dimension_x),
    function() {
      evaluate_kernel(out_of_range, dimension_x[, 1:2], dimension_y[, 1:2])
    },
    function() kernel_diagonal(out_of_range, dimension_x),
    function() kernel_gradient(out_of_range, dimension_x),
    function() kernel_input_gradient(out_of_range, dimension_x)
  )) {
    expect_error(
      call(),
      "selects column(s) 5, but its inputs have",
      fixed = TRUE,
      class = "gaussianprocesses_validation_error"
    )
  }

  # ARD lengths are checked against the number of selected columns.
  expect_error(
    evaluate_kernel(
      select_dimensions(rbf_kernel(length_scale = c(1, 2, 3)), 1:2),
      dimension_x
    ),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("models check selections and ARD lengths before fitting", {
  x <- dimension_x[, 1:3]
  y <- sin(x[, 1]) + x[, 2]
  ard_mismatch <- sum_kernel(
    select_dimensions(rbf_kernel(), 1),
    select_dimensions(rbf_kernel(length_scale = c(1, 2, 3)), 2:3)
  )
  out_of_range <- select_dimensions(rbf_kernel(), c(1, 4))

  fitters <- list(
    fit_gp = function(kernel) fit_gp(x, y, kernel, noise_variance = 0.1),
    optimize_gp = function(kernel) {
      optimize_gp(x, y, kernel, noise_variance = 0.1, n_starts = 1)
    },
    fit_sparse_gp = function(kernel) {
      fit_sparse_gp(x, y, kernel, noise_variance = 0.1, n_inducing = 3)
    },
    fit_heteroscedastic_gp = function(kernel) {
      fit_heteroscedastic_gp(x, y, kernel, noise_kernel = rbf_kernel())
    },
    noise_kernel = function(kernel) {
      fit_heteroscedastic_gp(x, y, rbf_kernel(), noise_kernel = kernel)
    },
    sample_gp_prior = function(kernel) sample_gp_prior(x, kernel, seed = 1),
    simulate_gp_data = function(kernel) {
      simulate_gp_data(x, kernel, noise_variance = 0.1, seed = 1)
    }
  )

  for (name in names(fitters)) {
    expect_error(
      fitters[[name]](ard_mismatch),
      "Kernel parameter 'kernel2.length_scale' has 3 values",
      fixed = TRUE,
      class = "gaussianprocesses_validation_error",
      info = name
    )
    expect_error(
      fitters[[name]](out_of_range),
      "selects column(s) 4",
      fixed = TRUE,
      class = "gaussianprocesses_validation_error",
      info = name
    )
  }

  # Kernels without selections get the same path-aware check.
  expect_error(
    fit_gp(x, y, scale_kernel(matern32_kernel(length_scale = c(1, 2)), 2)),
    "Kernel parameter 'kernel.length_scale' has 2 values",
    fixed = TRUE,
    class = "gaussianprocesses_validation_error"
  )
})


test_that("models with restricted kernels equal models on selected columns", {
  set.seed(25)
  x <- cbind(runif(30, -2, 2), runif(30, -2, 2), runif(30, -2, 2))
  y <- sin(2 * x[, 1]) + 0.5 * x[, 3]^2 + rnorm(30, sd = 0.1)
  x_new <- cbind(c(-1, 0.2, 1.5), c(0, 1, -1), c(0.5, -0.5, 1.8))
  inner <- matern52_kernel(length_scale = c(0.8, 1.5))
  restricted <- select_dimensions(inner, c(1, 3))

  model <- fit_gp(x, y, restricted, noise_variance = 0.01)
  reference <- fit_gp(x[, c(1, 3)], y, inner, noise_variance = 0.01)
  prediction <- predict_gp(model, x_new, include_covariance = TRUE)
  expected <- predict_gp(
    reference,
    x_new[, c(1, 3)],
    include_covariance = TRUE
  )

  for (field in c("mean", "latent_variance", "latent_covariance")) {
    expect_equal(prediction[[field]], expected[[field]], tolerance = 1e-12)
  }
  expect_equal(
    log_marginal_likelihood(model),
    log_marginal_likelihood(reference),
    tolerance = 1e-12
  )
  expect_equal(
    log_marginal_likelihood_gradient(model),
    log_marginal_likelihood_gradient(reference),
    tolerance = 1e-12
  )

  draws <- sample_gp_posterior(model, x_new, n_draws = 4, seed = 3)
  expected_draws <- sample_gp_posterior(
    reference,
    x_new[, c(1, 3)],
    n_draws = 4,
    seed = 3
  )
  expect_true(is.matrix(draws$latent))
  expect_equal(draws$latent, expected_draws$latent, tolerance = 1e-10)

  gradient <- predict_gradient_gp(model, x_new)
  expected_gradient <- predict_gradient_gp(reference, x_new[, c(1, 3)])
  expect_equal(unname(gradient$mean[, c(1, 3)]), unname(expected_gradient$mean))
  expect_identical(unname(gradient$mean[, 2]), c(0, 0, 0))
  expect_identical(unname(gradient$variance[, 2]), c(0, 0, 0))

  inducing <- x[c(2, 9, 15, 22, 28), ]
  sparse <- fit_sparse_gp(
    x,
    y,
    restricted,
    noise_variance = 0.01,
    inducing_points = inducing
  )
  sparse_reference <- fit_sparse_gp(
    x[, c(1, 3)],
    y,
    inner,
    noise_variance = 0.01,
    inducing_points = inducing[, c(1, 3)]
  )
  expect_equal(
    predict_sparse_gp(sparse, x_new)$mean,
    predict_sparse_gp(sparse_reference, x_new[, c(1, 3)])$mean,
    tolerance = 1e-12
  )

  heteroscedastic <- suppressWarnings(
    fit_heteroscedastic_gp(
      x,
      y,
      restricted,
      noise_kernel = select_dimensions(rbf_kernel(length_scale = 2), 2),
      max_iterations = 3
    )
  )
  expect_length(predict_heteroscedastic_gp(heteroscedastic, x_new)$mean, 3L)
})


test_that("optimization estimates the parameters of an additive kernel", {
  set.seed(52)
  x <- cbind(runif(60, -2, 2), runif(60, -2, 2), runif(60, -2, 2))
  y <- sin(2 * x[, 1]) + 0.5 * x[, 2] * x[, 3] + rnorm(60, sd = 0.1)
  kernel <- sum_kernel(
    select_dimensions(rbf_kernel(), 1),
    select_dimensions(rbf_kernel(length_scale = c(1, 1)), 2:3)
  )

  model <- optimize_gp(x, y, kernel, noise_variance = 0.05, n_starts = 2)

  expect_true(model$optimization$converged)
  expect_identical(
    names(model$optimization$optimized_parameters),
    c(
      "kernel1.variance", "kernel1.length_scale", "kernel2.variance",
      "kernel2.length_scale[1]", "kernel2.length_scale[2]", "noise_variance"
    )
  )
  expect_gt(
    log_marginal_likelihood(model),
    log_marginal_likelihood(fit_gp(x, y, kernel, noise_variance = 0.05))
  )
  expect_identical(model$kernel$children[[2]]$columns, 2:3)
  expect_lt(model$noise_variance, 0.05)
})


test_that("restricted kernels print their columns", {
  kernel <- sum_kernel(
    select_dimensions(rbf_kernel(length_scale = 0.5), 1),
    select_dimensions(linear_kernel(), c(3, 2))
  )

  output <- utils::capture.output(print(kernel))

  expect_identical(
    output,
    c(
      "SumKernel(",
      "  SelectDimensions(columns=1,",
      "    RBF(variance=1, length_scale=0.5)",
      "  )",
      "  SelectDimensions(columns=c(3, 2),",
      "    Linear(variance=1)",
      "  )",
      ")"
    )
  )
})
