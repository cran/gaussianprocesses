# Central finite differences of log_marginal_likelihood() with respect to
# log kernel parameters and, for homoscedastic models, the log noise
# variance, refitting the model at each perturbed parameter.
finite_difference_lml_gradient <- function(model, step = 1e-5) {
  kernel_values <- kernel_parameters(model$kernel, flatten = TRUE)
  homoscedastic <- identical(model$noise_structure, "homoscedastic")
  parameters <- if (homoscedastic) {
    c(kernel_values, noise_variance = model$noise_variance)
  } else {
    kernel_values
  }

  lml_at <- function(index, direction) {
    shifted <- parameters
    shifted[index] <- shifted[index] * exp(direction * step)
    noise <- if (homoscedastic) {
      unname(shifted[["noise_variance"]])
    } else {
      model$noise_variance
    }

    log_marginal_likelihood(
      fit_gp(
        model$x,
        model$y,
        kernel = update_kernel_parameters(
          model$kernel,
          shifted[names(kernel_values)]
        ),
        noise_variance = noise,
        mean = model$mean
      )
    )
  }

  derivatives <- vapply(
    seq_along(parameters),
    function(index) (lml_at(index, 1) - lml_at(index, -1)) / (2 * step),
    numeric(1)
  )
  names(derivatives) <- names(parameters)
  derivatives
}


lml_gradient_models <- function() {
  x_1d <- seq(-2, 2, length.out = 18)
  y_1d <- sin(2 * x_1d) + 0.2 * cos(5 * x_1d)
  x_2d <- cbind(seq(-1.5, 1.5, length.out = 20), cos(1:20))
  y_2d <- sin(2 * x_2d[, 1]) + 0.3 * x_2d[, 2]
  x_3d <- cbind(
    seq(-1.5, 1.5, length.out = 24),
    sin(1:24),
    cos(2 * (1:24))
  )
  y_3d <- sin(2 * x_3d[, 1]) + 0.4 * x_3d[, 2] - 0.2 * x_3d[, 3]

  list(
    rbf = fit_gp(x_1d, y_1d, rbf_kernel(1.2, 0.6), noise_variance = 0.02),
    matern12 = fit_gp(
      x_1d,
      y_1d,
      matern12_kernel(0.8, 1.3),
      noise_variance = 0.05
    ),
    rational_quadratic = fit_gp(
      x_1d,
      y_1d,
      rational_quadratic_kernel(1, 0.7, alpha = 0.6),
      noise_variance = 0.01,
      mean = constant_mean(0.2)
    ),
    ard_matern52 = fit_gp(
      x_2d,
      y_2d,
      matern52_kernel(0.9, c(0.6, 2)),
      noise_variance = 0.03
    ),
    composite = fit_gp(
      x_2d,
      y_2d,
      sum_kernel(
        scale_kernel(rbf_kernel(0.8, c(0.7, 1.5)), scale = 1.4),
        product_kernel(
          periodic_kernel(0.5, 1, period = 1.7),
          linear_kernel(0.3)
        )
      ),
      noise_variance = 0.02
    ),
    known_noise = fit_gp(
      x_1d,
      y_1d,
      matern32_kernel(1.1, 0.9),
      noise_variance = seq(0.01, 0.1, length.out = length(x_1d))
    ),
    matern32 = fit_gp(
      x_1d,
      y_1d,
      matern32_kernel(0.7, 0.5),
      noise_variance = 0.04
    ),
    periodic = fit_gp(
      x_1d,
      y_1d,
      periodic_kernel(1, 0.9, period = 1.4),
      noise_variance = 0.03
    ),
    shared_matern12_3d = fit_gp(
      x_3d,
      y_3d,
      matern12_kernel(1.3, 1.1),
      noise_variance = 0.05
    ),
    ard_rbf_3d = fit_gp(
      x_3d,
      y_3d,
      rbf_kernel(1, c(0.6, 1.4, 2.5)),
      noise_variance = 0.02
    ),
    ard_rational_quadratic = fit_gp(
      x_2d,
      y_2d,
      rational_quadratic_kernel(0.8, c(0.5, 1.7), alpha = 1.3),
      noise_variance = 0.02
    ),
    rbf_plus_white_noise = fit_gp(
      x_2d,
      y_2d,
      sum_kernel(rbf_kernel(1.1, c(0.8, 1.2)), white_noise_kernel(0.05)),
      noise_variance = 0.01
    ),
    linear_plus_rbf = fit_gp(
      x_1d,
      y_1d,
      sum_kernel(linear_kernel(0.2), rbf_kernel(0.9, 0.5)),
      noise_variance = 0.02
    ),
    scaled_product = fit_gp(
      x_2d,
      y_2d,
      scale_kernel(
        product_kernel(
          rbf_kernel(1, c(0.9, 2)),
          periodic_kernel(0.8, 1.1, period = 1.5)
        ),
        scale = 0.7
      ),
      noise_variance = 0.03
    )
  )
}


simulated_ard_data <- function() {
  x <- cbind(seq(-2, 2, length.out = 40), cos(1:40))
  y <- simulate_gp_data(
    x,
    kernel = rbf_kernel(variance = 1, length_scale = c(0.8, 1.5)),
    noise_variance = 0.02,
    seed = 11
  )$observed

  list(x = x, y = y)
}


test_that("log marginal likelihood gradient matches finite differences", {
  models <- lml_gradient_models()

  for (name in names(models)) {
    model <- models[[name]]
    analytical <- log_marginal_likelihood_gradient(model)
    numerical <- finite_difference_lml_gradient(model)

    expect_identical(names(analytical), names(numerical), info = name)
    expect_lte(
      max(abs(analytical - numerical)) / max(1, max(abs(analytical))),
      1e-6,
      label = name
    )
  }
})


test_that("the trace term matches Cholesky solves for each parameter", {
  # The gradient obtains tr(C^{-1} dC/d eta_j) from the precision matrix
  # computed once from the Cholesky factor. Solving C X_j = dC/d eta_j with
  # the same factor for each parameter must give the same gradient.
  model <- lml_gradient_models()$scaled_product
  factor <- model$cholesky
  derivatives <- kernel_gradient(model$kernel, model$x)
  derivatives$noise_variance <- diag(model$noise_variance, nrow(model$x))

  by_solves <- vapply(
    derivatives,
    function(derivative) {
      solved <- backsolve(factor, forwardsolve(t(factor), derivative))
      0.5 * (
        drop(crossprod(model$alpha, derivative %*% model$alpha)) -
          sum(diag(solved))
      )
    },
    numeric(1)
  )

  expect_equal(
    log_marginal_likelihood_gradient(model),
    by_solves,
    tolerance = 1e-10
  )
})


test_that("gradient names follow kernel paths and the noise parameter", {
  models <- lml_gradient_models()

  expect_identical(
    names(log_marginal_likelihood_gradient(models$composite)),
    c(
      names(kernel_parameters(models$composite$kernel, flatten = TRUE)),
      "noise_variance"
    )
  )

  # Observation-specific noise is fixed data, not a parameter.
  expect_identical(
    names(log_marginal_likelihood_gradient(models$known_noise)),
    names(kernel_parameters(models$known_noise$kernel, flatten = TRUE))
  )
})


test_that("gradient is finite when the fit needed numerical jitter", {
  model <- fit_gp(
    c(0, 0, 1, 2),
    c(1, 1, 0.5, -0.2),
    kernel = rbf_kernel(),
    noise_variance = 0
  )

  expect_gt(model$jitter, 0)
  expect_true(all(is.finite(log_marginal_likelihood_gradient(model))))
})


test_that("analytical and numerical gradients reach the same optimum", {
  # Noisy data from a known ARD kernel, so the optimum is interior and
  # identifiable rather than on a bound or a flat ridge.
  data <- simulated_ard_data()
  x <- data$x
  y <- data$y
  kernel <- rbf_kernel(variance = 1, length_scale = c(1, 1))

  analytical <- optimize_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = 0.05,
    n_starts = 1
  )
  numerical <- optimize_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = 0.05,
    n_starts = 1,
    gradient = "numerical"
  )

  expect_true(analytical$optimization$converged)

  expect_identical(analytical$optimization$gradient, "analytical")
  expect_identical(numerical$optimization$gradient, "numerical")
  expect_equal(
    analytical$optimization$log_marginal_likelihood,
    numerical$optimization$log_marginal_likelihood,
    tolerance = 1e-6
  )
  expect_equal(
    analytical$optimization$optimized_parameters,
    numerical$optimization$optimized_parameters,
    tolerance = 1e-3
  )

  # At an interior optimum the gradient vanishes.
  expect_lt(max(abs(log_marginal_likelihood_gradient(analytical))), 1e-2)
})


test_that("evaluation counts include finite-difference objective calls", {
  data <- simulated_ard_data()
  n_parameters <- 4L

  for (gradient in c("analytical", "numerical")) {
    optimization <- optimize_gp(
      data$x,
      data$y,
      kernel = rbf_kernel(variance = 1, length_scale = c(1, 1)),
      noise_variance = 0.05,
      n_starts = 2,
      gradient = gradient
    )$optimization

    for (start in optimization$start_results) {
      counts <- start$counts

      if (identical(gradient, "analytical")) {
        # Value and gradient are requested together and share one fit.
        expect_identical(
          start$function_evaluations,
          unname(counts[["function"]])
        )
      } else {
        # optim() does not count the two objective calls per parameter it
        # makes for each finite-difference gradient.
        expect_identical(
          start$function_evaluations,
          unname(counts[["function"]] + 2L * n_parameters * counts[["gradient"]])
        )
      }

      expect_identical(
        start$gradient_evaluations,
        unname(counts[["gradient"]])
      )
      expect_lte(start$model_fits, start$function_evaluations)
    }

    for (field in c(
      "function_evaluations",
      "gradient_evaluations",
      "model_fits"
    )) {
      expect_identical(
        optimization[[field]],
        sum(vapply(optimization$start_results, `[[`, integer(1), field)),
        info = paste(gradient, field)
      )
    }
  }
})


test_that("multi-start optimization is deterministic for both gradients", {
  data <- simulated_ard_data()

  run <- function(gradient) {
    optimize_gp(
      data$x,
      data$y,
      kernel = rbf_kernel(variance = 1, length_scale = c(1, 1)),
      noise_variance = 0.05,
      n_starts = 3,
      gradient = gradient
    )$optimization
  }

  first <- run("analytical")
  second <- run("analytical")
  numerical <- run("numerical")

  expect_identical(first$starts, numerical$starts)
  expect_identical(first$optimized_parameters, second$optimized_parameters)
  expect_identical(first$function_evaluations, second$function_evaluations)
  expect_identical(first$best_start, second$best_start)
})


test_that("failed fits give a finite penalty and a zero gradient", {
  # Duplicated inputs without noise or jitter cannot be factorized.
  problem <- gaussianprocesses:::.optimization_problem(
    x = matrix(c(0, 0, 1), ncol = 1),
    y = c(1, 1, 0),
    kernel = rbf_kernel(),
    mean = zero_mean(),
    optimize_noise = FALSE,
    fixed_noise_variance = 0,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 1L,
    tolerance = sqrt(.Machine$double.eps)
  )
  theta <- c(variance = 0, length_scale = 0)

  expect_identical(problem$objective(theta), 1e100)
  expect_identical(problem$gradient(theta), c(0, 0))
})


test_that("optimized model printout names the gradient method", {
  model <- optimize_gp(
    seq(-1, 1, length.out = 10),
    sin(seq(-1, 1, length.out = 10)),
    kernel = rbf_kernel(),
    noise_variance = 0.01,
    n_starts = 1
  )

  output <- capture.output(print(model))

  expect_true(any(grepl("analytical gradient", output, fixed = TRUE)))
  expect_true(any(grepl("evaluations:", output, fixed = TRUE)))
  expect_error(
    optimize_gp(1:5, 1:5, kernel = rbf_kernel(), gradient = "exact")
  )
})
