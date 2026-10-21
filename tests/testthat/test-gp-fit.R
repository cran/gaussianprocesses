test_that("exact GP fit matches a direct linear solve on a toy problem", {
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

  expected_covariance <- evaluate_kernel(kernel, x) +
    diag(noise_variance, length(x))
  expected_alpha <- solve(expected_covariance, y)

  expect_s3_class(model, "gaussianprocesses_model")
  expect_equal(model$alpha, expected_alpha, tolerance = 1e-10)
  expect_equal(model$jitter, 0)
  expect_equal(model$n_observations, 3L)
  expect_equal(model$n_features, 1L)
})


test_that("constant prior mean is removed before conditioning", {
  x <- c(0, 1, 2)
  y <- c(3, 4, 5)
  mean <- constant_mean(3)
  kernel <- matern32_kernel(variance = 1, length_scale = 2)

  model <- fit_gp(
    x,
    y,
    kernel = kernel,
    mean = mean,
    noise_variance = 0.1
  )

  expected_centered <- c(0, 1, 2)
  expected_covariance <- evaluate_kernel(kernel, x) + diag(0.1, 3)

  expect_equal(model$prior_mean, rep(3, 3))
  expect_equal(model$centered_y, expected_centered)
  expect_equal(
    model$alpha,
    solve(expected_covariance, expected_centered),
    tolerance = 1e-10
  )
})


test_that("exact GP accepts multidimensional matrix inputs", {
  x <- matrix(
    c(
      0, 0,
      1, 0,
      0, 1,
      1, 1
    ),
    ncol = 2,
    byrow = TRUE
  )
  y <- c(0, 1, 1, 2)

  model <- fit_gp(
    x,
    y,
    kernel = rbf_kernel(length_scale = 1.5),
    noise_variance = 0.05
  )

  expect_equal(dim(model$x), c(4L, 2L))
  expect_equal(model$n_features, 2L)
  expect_length(model$alpha, 4L)
})


test_that("observation noise and numerical jitter remain separate", {
  x <- c(0, 0)
  y <- c(1, 1)

  model <- fit_gp(
    x,
    y,
    kernel = rbf_kernel(),
    noise_variance = 0,
    fallback_jitter = 1e-10
  )

  expect_equal(model$noise_variance, 0)
  expect_gt(model$jitter, 0)

  # The statistical covariance excludes the jitter; the factor includes it.
  covariance <- gaussianprocesses:::.model_observation_covariance(model)
  expect_equal(covariance, matrix(1, 2, 2))
  expect_equal(
    crossprod(model$cholesky),
    covariance + diag(model$jitter, 2)
  )
})


test_that("positive observation noise resolves duplicate observations statistically", {
  x <- c(0, 0)
  y <- c(1, 1)

  model <- fit_gp(
    x,
    y,
    kernel = rbf_kernel(),
    noise_variance = 0.01
  )

  expect_equal(model$jitter, 0)
  expect_equal(
    diag(gaussianprocesses:::.model_observation_covariance(model)),
    c(1.01, 1.01)
  )
})


test_that("training data and model arguments are validated", {
  expect_error(
    fit_gp(
      c(0, 1, 2),
      c(1, 2),
      kernel = rbf_kernel()
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    fit_gp(
      c(0, 1),
      c(1, 2),
      kernel = "rbf"
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    fit_gp(
      c(0, 1),
      c(1, 2),
      kernel = rbf_kernel(),
      noise_variance = -0.1
    ),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("fitted model exposes the kernel and numerical state", {
  model <- fit_gp(
    c(0, 1, 2),
    c(1, 0, 1),
    kernel = sum_kernel(
      rbf_kernel(variance = 2, length_scale = 0.5),
      white_noise_kernel(variance = 0.05)
    ),
    noise_variance = 0.1
  )

  output <- capture.output(print(model))

  expect_true(any(grepl("Exact Gaussian-process", output, fixed = TRUE)))
  expect_true(any(grepl("observations: 3", output, fixed = TRUE)))
  expect_true(any(grepl("SumKernel", output, fixed = TRUE)))
  expect_true(any(grepl("noise variance", output, fixed = TRUE)))
  expect_true(any(grepl("numerical jitter", output, fixed = TRUE)))
})
