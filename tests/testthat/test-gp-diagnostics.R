test_that("LOO summaries match brute-force refitting", {
  x <- c(0, 0.7, 1.4, 2.1)
  y <- c(0.2, 0.8, -0.1, 0.4)
  kernel <- rbf_kernel(
    variance = 1.3,
    length_scale = 0.9
  )
  noise_variance <- 0.15

  model <- fit_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = noise_variance
  )
  loo <- loo_gp(model, interval_level = 0.9)

  for (index in seq_along(x)) {
    training_index <- setdiff(seq_along(x), index)
    reduced_model <- fit_gp(
      x[training_index],
      y[training_index],
      kernel = kernel,
      noise_variance = noise_variance
    )
    prediction <- predict_gp(
      reduced_model,
      x[index],
      interval_level = 0.9
    )

    expect_equal(
      loo$mean[index],
      prediction$mean,
      tolerance = 1e-9
    )
    expect_equal(
      loo$variance[index],
      prediction$observation_variance,
      tolerance = 1e-9
    )
  }
})


test_that("standardized residuals and PIT are internally consistent", {
  model <- fit_gp(
    c(0, 1, 2, 3),
    c(0.2, 0.9, 0.1, 1.1),
    kernel = matern32_kernel(
      variance = 1,
      length_scale = 1
    ),
    noise_variance = 0.1
  )
  loo <- loo_gp(model)

  expect_equal(
    loo$standardized_residual,
    loo$residual / loo$sd,
    tolerance = 1e-12
  )
  expect_equal(
    loo$pit,
    stats::pnorm(
      loo$observed,
      mean = loo$mean,
      sd = loo$sd
    ),
    tolerance = 1e-12
  )
  expect_true(all(loo$pit >= 0 & loo$pit <= 1))
})


test_that("coverage diagnostics agree with pointwise LOO coverage", {
  model <- fit_gp(
    c(0, 1, 2, 3, 4),
    c(0, 0.7, 0.3, 1, 0.6),
    kernel = rbf_kernel(length_scale = 1.2),
    noise_variance = 0.2
  )

  loo <- loo_gp(model, interval_level = 0.8)
  calibration <- gp_calibration_diagnostics(
    model,
    interval_level = 0.8
  )

  expect_equal(
    calibration$empirical_coverage,
    mean(loo$covered)
  )
  expect_equal(
    calibration$coverage_error,
    mean(loo$covered) - 0.8
  )
  expect_equal(
    calibration$mean_standardized_residual,
    mean(loo$standardized_residual)
  )
})


test_that("numerical diagnostics keep noise and jitter separate", {
  model <- fit_gp(
    c(0, 0),
    c(1, 1),
    kernel = rbf_kernel(),
    noise_variance = 0,
    fallback_jitter = 1e-10
  )

  diagnostics <- gp_numerical_diagnostics(model)

  expect_equal(diagnostics$observation_noise_variance, 0)
  expect_gt(diagnostics$numerical_jitter, 0)
  expect_true(diagnostics$used_numerical_jitter)
  expect_gte(diagnostics$cholesky_attempts, 2L)
  expect_true(
    is.finite(diagnostics$reciprocal_condition_number)
  )
})


test_that("well-conditioned fits report no numerical jitter", {
  model <- fit_gp(
    c(0, 1, 2),
    c(0, 1, 0),
    kernel = rbf_kernel(
      variance = 1,
      length_scale = 0.5
    ),
    noise_variance = 0.2
  )

  diagnostics <- gp_numerical_diagnostics(model)

  expect_false(diagnostics$used_numerical_jitter)
  expect_equal(diagnostics$numerical_jitter, 0)
  expect_gt(diagnostics$reciprocal_condition_number, 0)
  expect_true(is.finite(diagnostics$condition_number))
})


test_that("correctly specified simulations give plausible calibration summaries", {
  set.seed(1729)

  x <- seq(-2, 2, length.out = 60)
  kernel <- rbf_kernel(
    variance = 1.2,
    length_scale = 0.7
  )
  noise_variance <- 0.15
  covariance <- evaluate_kernel(kernel, x) +
    diag(noise_variance, length(x))

  y <- drop(
    t(chol(covariance)) %*% stats::rnorm(length(x))
  )

  model <- fit_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = noise_variance
  )
  calibration <- gp_calibration_diagnostics(
    model,
    interval_level = 0.9
  )

  expect_lt(
    abs(calibration$mean_standardized_residual),
    0.4
  )
  expect_gt(calibration$sd_standardized_residual, 0.5)
  expect_lt(calibration$sd_standardized_residual, 1.5)
  expect_gt(calibration$empirical_coverage, 0.75)
  expect_lte(calibration$empirical_coverage, 1)
  expect_gt(calibration$pit_mean, 0.3)
  expect_lt(calibration$pit_mean, 0.7)
})


test_that("comprehensive diagnostics keep model and numerical summaries separate", {
  model <- fit_gp(
    c(0, 1, 2, 3),
    c(0.1, 0.8, 0.2, 0.9),
    kernel = matern52_kernel(
      variance = 1,
      length_scale = 1.1
    ),
    noise_variance = 0.1
  )

  diagnostics <- gp_diagnostics(model)

  expect_s3_class(
    diagnostics,
    "gaussianprocesses_diagnostics"
  )
  expect_true(is.data.frame(diagnostics$loo))
  expect_true(is.list(diagnostics$calibration))
  expect_true(is.list(diagnostics$numerical))
  expect_true(
    "empirical_coverage" %in% names(diagnostics$calibration)
  )
  expect_true(
    "reciprocal_condition_number" %in%
      names(diagnostics$numerical)
  )
})


test_that("diagnostic arguments are validated", {
  model <- fit_gp(
    c(0, 1, 2),
    c(0, 1, 0),
    kernel = rbf_kernel(),
    noise_variance = 0.1
  )

  expect_error(
    loo_gp(model, interval_level = 1),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    gp_numerical_diagnostics(
      model,
      condition_tolerance = 0
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    gp_diagnostics("not a model"),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("diagnostics have a concise print method", {
  model <- fit_gp(
    c(0, 1, 2, 3),
    c(0, 1, 0, 1),
    kernel = rbf_kernel(),
    noise_variance = 0.1
  )

  diagnostics <- gp_diagnostics(model)
  output <- capture.output(print(diagnostics))

  expect_true(
    any(grepl("Gaussian-process diagnostics", output, fixed = TRUE))
  )
  expect_true(
    any(grepl("LOO empirical coverage", output, fixed = TRUE))
  )
  expect_true(
    any(grepl("reciprocal condition number", output, fixed = TRUE))
  )
})


test_that("gp_diagnostics() computes leave-one-out once", {
  x <- seq(-2, 2, length.out = 15)
  model <- fit_gp(x, sin(2 * x), rbf_kernel(1, 0.6), noise_variance = 0.01)
  original <- gaussianprocesses::loo_gp
  calls <- 0L

  local_mocked_bindings(
    loo_gp = function(...) {
      calls <<- calls + 1L
      original(...)
    }
  )
  diagnostics <- gp_diagnostics(model, interval_level = 0.9)

  expect_identical(calls, 1L)
  expect_identical(
    diagnostics$calibration,
    gp_calibration_diagnostics(model, interval_level = 0.9)
  )
  expect_identical(diagnostics$loo, loo_gp(model, interval_level = 0.9))
})
