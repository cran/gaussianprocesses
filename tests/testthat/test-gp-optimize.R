test_that("log marginal likelihood matches the direct Gaussian formula", {
  x <- c(0, 1, 2)
  y <- c(1, -0.5, 0.25)
  kernel <- rbf_kernel(variance = 1.4, length_scale = 0.8)
  noise_variance <- 0.2

  model <- fit_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = noise_variance
  )

  covariance <- evaluate_kernel(kernel, x) +
    diag(noise_variance, length(x))
  direct <- -0.5 * drop(crossprod(y, solve(covariance, y))) -
    0.5 * as.numeric(determinant(covariance, logarithm = TRUE)$modulus) -
    0.5 * length(y) * log(2 * pi)

  expect_equal(
    log_marginal_likelihood(model),
    direct,
    tolerance = 1e-10
  )
})


test_that("optimization does not reduce the best log marginal likelihood", {
  x <- seq(-1.5, 1.5, length.out = 9)
  y <- sin(2 * x) + 0.15 * x
  initial_kernel <- rbf_kernel(
    variance = 0.3,
    length_scale = 3
  )

  initial_model <- fit_gp(
    x,
    y,
    kernel = initial_kernel,
    noise_variance = 0.05
  )

  optimized <- optimize_gp(
    x,
    y,
    kernel = initial_kernel,
    noise_variance = 0.05,
    optimize_noise = FALSE,
    n_starts = 3,
    control = list(maxit = 40)
  )

  expect_gte(
    log_marginal_likelihood(optimized),
    log_marginal_likelihood(initial_model) - 1e-8
  )
  expect_true(!is.null(optimized$optimization))
  expect_equal(optimized$optimization$n_starts, 3L)
})


test_that("noise variance can be optimized separately from kernel parameters", {
  x <- seq(0, 1, length.out = 8)
  y <- c(0.05, 0.31, 0.51, 0.79, 0.91, 1.06, 0.96, 0.84)

  model <- optimize_gp(
    x,
    y,
    kernel = matern32_kernel(variance = 1, length_scale = 0.4),
    noise_variance = 0.1,
    optimize_noise = TRUE,
    n_starts = 2,
    control = list(maxit = 30)
  )

  expect_gt(model$noise_variance, 0)
  expect_true(
    "noise_variance" %in%
      names(model$optimization$optimized_parameters)
  )
  expect_true(is.finite(model$optimization$log_marginal_likelihood))
})


test_that("zero initial noise is lifted to the positive optimization bound", {
  x <- c(0, 1, 2, 3)
  y <- c(0, 1, 0.5, 1.2)

  model <- optimize_gp(
    x,
    y,
    kernel = rbf_kernel(),
    noise_variance = 0,
    optimize_noise = TRUE,
    n_starts = 1,
    lower = 1e-6,
    control = list(maxit = 5)
  )

  expect_gte(
    model$optimization$initial_parameters[["noise_variance"]],
    1e-6
  )
  expect_gt(model$noise_variance, 0)
})


test_that("nested kernel parameter paths survive optimization", {
  x <- c(0, 0.5, 1, 1.5, 2)
  y <- c(0.1, 0.4, 0.2, -0.1, 0.05)

  kernel <- sum_kernel(
    rbf_kernel(variance = 1, length_scale = 0.8),
    white_noise_kernel(variance = 0.1)
  )

  model <- optimize_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = 0.05,
    optimize_noise = FALSE,
    n_starts = 1,
    control = list(maxit = 10)
  )

  paths <- names(kernel_parameters(model$kernel, flatten = TRUE))

  expect_equal(
    paths,
    c(
      "kernel1.variance",
      "kernel1.length_scale",
      "kernel2.variance"
    )
  )
})


test_that("deterministic multi-start construction is reproducible", {
  initial <- log(c(variance = 1, length_scale = 2))

  first <- gaussianprocesses:::.make_optimization_starts(
    initial,
    n_starts = 5,
    spread = 0.5,
    lower = log(1e-4),
    upper = log(1e4)
  )
  second <- gaussianprocesses:::.make_optimization_starts(
    initial,
    n_starts = 5,
    spread = 0.5,
    lower = log(1e-4),
    upper = log(1e4)
  )

  expect_equal(first, second)
  expect_equal(first[1, ], initial)
})


test_that("optimization arguments and bounds are validated", {
  x <- c(0, 1, 2)
  y <- c(0, 1, 0)

  expect_error(
    optimize_gp(
      x,
      y,
      kernel = rbf_kernel(),
      lower = 1,
      upper = 1
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    optimize_gp(
      x,
      y,
      kernel = rbf_kernel(),
      n_starts = 0
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    optimize_gp(
      x,
      y,
      kernel = rbf_kernel(),
      optimize_noise = NA
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    optimize_gp(
      x,
      y,
      kernel = rbf_kernel(variance = 1e-10),
      lower = 1e-8
    ),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("the default iteration limit lets periodic composites converge", {
  expect_identical(eval(formals(optimize_gp)$control), list(maxit = 500))

  # This start needs more than 100 L-BFGS-B iterations; the earlier default
  # of 100 stopped it on the iteration limit (convergence code 1).
  x <- cbind(seq(-2, 2, length.out = 60), cos(1:60))
  y <- simulate_gp_data(
    x,
    kernel = sum_kernel(
      rbf_kernel(0.8, c(0.6, 2)),
      periodic_kernel(0.4, 1, period = 1.3)
    ),
    noise_variance = 0.02,
    seed = 7
  )$observed

  model <- optimize_gp(
    x,
    y,
    kernel = sum_kernel(
      rbf_kernel(1, c(1, 1)),
      periodic_kernel(0.5, 1, period = 1.2)
    ),
    noise_variance = 0.05,
    n_starts = 1
  )

  expect_true(model$optimization$converged)
  expect_gt(model$optimization$counts[[1]], 100)
})


test_that("optimized model printout reports optimization state", {
  model <- optimize_gp(
    c(0, 1, 2, 3),
    c(0, 1, 0, 1),
    kernel = rbf_kernel(),
    noise_variance = 0.1,
    optimize_noise = FALSE,
    n_starts = 1,
    control = list(maxit = 5)
  )

  output <- capture.output(print(model))

  expect_true(any(grepl("optimization:", output, fixed = TRUE)))
  expect_true(any(grepl("log marginal likelihood:", output, fixed = TRUE)))
})
