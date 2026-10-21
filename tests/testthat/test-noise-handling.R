# Observation noise in optimize_gp(), rolling_origin_gp(), and
# benchmark_sparse_gp(): a scalar variance, a fixed per-observation vector,
# or a known relative pattern with an estimated scale.

noise_x <- seq(0, 10, length.out = 60)
noise_pattern <- ifelse(noise_x > 5, 0.08, 0.02)


test_that("a fixed noise vector is kept by optimize_gp()", {
  y <- simulate_gp_data(
    noise_x,
    rbf_kernel(1, 1.2),
    noise_variance = noise_pattern,
    seed = 3
  )$observed

  model <- optimize_gp(
    noise_x,
    y,
    kernel = rbf_kernel(1, 1),
    noise_variance = noise_pattern,
    optimize_noise = FALSE,
    n_starts = 1
  )

  expect_identical(model$noise_variance, noise_pattern)
  expect_identical(model$noise_structure, "known_heteroscedastic")
  expect_identical(model$optimization$noise_parameterization, "fixed")
  optimized <- names(model$optimization$optimized_parameters)
  expect_false("noise_variance" %in% optimized)

  # At the optimum the model is the exact GP with the same parameters.
  refit <- fit_gp(noise_x, y, model$kernel, noise_variance = noise_pattern)
  expect_identical(refit$cholesky, model$cholesky)
})


test_that("the scale of a known noise pattern has an exact gradient", {
  y <- simulate_gp_data(
    noise_x,
    rbf_kernel(1, 1.2),
    noise_variance = 2 * noise_pattern,
    seed = 4
  )$observed
  kernel <- rbf_kernel(1.1, 1.3)
  lml_at_scale <- function(scale) {
    log_marginal_likelihood(
      fit_gp(noise_x, y, kernel, noise_variance = scale * noise_pattern)
    )
  }

  for (scale in c(0.5, 2, 7)) {
    model <- fit_gp(noise_x, y, kernel, noise_variance = scale * noise_pattern)
    analytical <- gaussianprocesses:::.log_marginal_likelihood_gradient(
      model,
      include_noise = TRUE,
      noise_parameter = "noise_scale"
    )[["noise_scale"]]

    # Central difference in log(scale).
    step <- 1e-5
    numerical <- (
      lml_at_scale(scale * exp(step)) - lml_at_scale(scale * exp(-step))
    ) / (2 * step)

    expect_equal(analytical, numerical, tolerance = 1e-6, info = scale)
  }
})


test_that("optimize_gp() estimates the scale of a known noise pattern", {
  # Simulated with noise 2.5 times the pattern. Over 30 seeds the estimated
  # log scale had mean bias 0.020 (Monte Carlo standard error 0.021) and
  # standard deviation 0.11, so one replicate lies within 0.5 of the truth.
  x <- seq(0, 10, length.out = 150)
  pattern <- ifelse(x > 5, 0.08, 0.02)
  y <- simulate_gp_data(
    x,
    rbf_kernel(1, 1.2),
    noise_variance = 2.5 * pattern,
    seed = 1
  )$observed

  model <- optimize_gp(
    x,
    y,
    kernel = rbf_kernel(1, 1),
    noise_variance = pattern,
    n_starts = 1
  )
  optimization <- model$optimization
  scale <- optimization$optimized_parameters[["noise_scale"]]

  expect_identical(optimization$noise_parameterization, "scale")
  expect_identical(optimization$coordinates[["noise_scale"]], "log")
  expect_identical(optimization$starts[1, "noise_scale"], 0)
  expect_true(optimization$converged)
  expect_lt(abs(log(scale / 2.5)), 0.5)
  expect_equal(model$noise_variance, scale * pattern)
  expect_identical(model$noise_structure, "known_heteroscedastic")

  # The noise-scale gradient is zero at the optimum.
  gradient <- gaussianprocesses:::.log_marginal_likelihood_gradient(
    model,
    include_noise = TRUE,
    noise_parameter = "noise_scale"
  )
  expect_lt(abs(gradient[["noise_scale"]]), 1e-2)

  # A scalar still means a homoscedastic noise variance.
  scalar <- optimize_gp(
    x,
    y,
    rbf_kernel(1, 1),
    noise_variance = 0.05,
    n_starts = 1
  )
  expect_identical(scalar$optimization$noise_parameterization, "variance")
  expect_true(
    "noise_variance" %in% names(scalar$optimization$optimized_parameters)
  )
})


test_that("noise arguments are validated", {
  y <- sin(noise_x)

  expect_error(
    optimize_gp(noise_x, y, rbf_kernel(), noise_variance = c(0.1, 0.2)),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    optimize_gp(noise_x, y, rbf_kernel(), noise_variance = rep(0, 60)),
    "needs a positive variance",
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    rolling_origin_gp(
      1:20,
      sin(1:20),
      rbf_kernel(),
      10,
      noise_variance = c(1, 2)
    ),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    benchmark_sparse_gp(noise_x, y, rbf_kernel(), 5, noise_variance = c(1, 2)),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("rolling-origin evaluation subsets per-observation noise", {
  time <- 1:30
  y <- sin(time / 4)
  noise <- seq(0.01, 0.3, length.out = 30)
  kernel <- rbf_kernel(1, 4)

  backtest <- rolling_origin_gp(
    time,
    y,
    kernel,
    initial_window = 20,
    horizon = 2,
    step = 4,
    noise_variance = noise
  )

  # The split at origin 24 must use rows 1-24 for training and the noise of
  # rows 25-26 for the forecast.
  model <- fit_time_series_gp(
    time[1:24],
    y[1:24],
    kernel,
    noise_variance = noise[1:24]
  )
  forecast <- forecast_gp(
    model,
    time[25:26],
    observation_noise_variance = noise[25:26]
  )
  rows <- backtest[backtest$origin_index == 24, ]

  expect_equal(rows$mean, forecast$mean)
  expect_equal(rows$observation_variance, forecast$observation_variance)

  # Index-dependent noise changes the forecast variance, so the subsetting
  # matters.
  homoscedastic <- rolling_origin_gp(
    time,
    y,
    kernel,
    initial_window = 20,
    horizon = 2,
    step = 4,
    noise_variance = 0.01
  )
  expect_true(
    all(backtest$observation_variance > homoscedastic$observation_variance)
  )
})


test_that("sparse benchmarks give both models the same noise", {
  x <- seq(0, 6, length.out = 40)
  y <- sin(x)
  noise <- seq(0.01, 0.1, length.out = 40)
  kernel <- rbf_kernel(1, 0.8)

  benchmark <- benchmark_sparse_gp(
    x,
    y,
    kernel,
    n_inducing = 10,
    noise_variance = noise,
    observation_noise_variance = noise
  )

  exact <- predict_gp(
    fit_gp(x, y, kernel, noise_variance = noise),
    x,
    observation_noise_variance = noise
  )
  sparse <- predict_sparse_gp(
    fit_sparse_gp(x, y, kernel, noise_variance = noise, n_inducing = 10),
    x,
    observation_noise_variance = noise
  )

  expect_equal(benchmark$exact_prediction$mean, exact$mean)
  expect_equal(benchmark$sparse_prediction$mean, sparse$mean)
  expect_equal(benchmark$exact_prediction$observation_noise_variance, noise)

  # Vector noise needs prediction-time noise, as in predict_gp().
  expect_error(
    benchmark_sparse_gp(x, y, kernel, n_inducing = 10, noise_variance = noise),
    class = "gaussianprocesses_validation_error"
  )
})
