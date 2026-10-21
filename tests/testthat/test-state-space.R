# State-space inference for one-dimensional Matern kernels --------------------

max_relative_error <- function(actual, expected) {
  max(abs(actual - expected) / abs(expected))
}

# Irregular inputs with ties among the training inputs and between training
# and prediction inputs.
state_space_data <- function(n = 120L, seed = 1L) {
  set.seed(seed)
  x <- sort(stats::runif(n, 0, 20))
  x[21:23] <- x[[20]]
  x[[n %/% 2L + 1L]] <- x[[n %/% 2L]]
  y <- sin(x) + 0.5 * cos(3 * x) + stats::rnorm(n, sd = 0.3)
  query <- c(25, sort(stats::runif(15, 0, 20)), -2, x[[40]], x[[20]])

  list(x = x, y = y, query = query, noise = stats::runif(n, 0.02, 0.2))
}

state_space_kernels <- list(
  matern12 = matern12_kernel(variance = 1.3, length_scale = 0.7),
  matern32 = matern32_kernel(variance = 0.8, length_scale = 1.9),
  matern52 = matern52_kernel(variance = 2.1, length_scale = 1.1),
  sum = sum_kernel(
    matern12_kernel(variance = 0.3, length_scale = 5),
    matern52_kernel(variance = 1, length_scale = 1.5),
    matern32_kernel(variance = 0.2, length_scale = 0.4)
  ),
  scaled = scale_kernel(
    sum_kernel(
      matern32_kernel(length_scale = 2),
      scale_kernel(matern12_kernel(length_scale = 0.3), scale = 0.2)
    ),
    scale = 1.7
  ),
  selected = select_dimensions(matern52_kernel(length_scale = 0.8), 1)
)


test_that("the Matern state-space forms solve the Lyapunov equation and reproduce the kernel", {
  constructors <- list(matern12_kernel, matern32_kernel, matern52_kernel)
  lags <- c(0, 1e-3, 0.05, 0.4, 1, 2.5, 7, 30)

  for (order in 0:2) {
    for (parameters in list(c(1, 1), c(1.7, 0.8), c(0.05, 12))) {
      kernel <- constructors[[order + 1L]](
        variance = parameters[[1L]],
        length_scale = parameters[[2L]]
      )
      form <- .matern_state_space(order)(kernel$parameters)
      stationary <- form$stationary
      lyapunov <- form$feedback %*% stationary +
        stationary %*% t(form$feedback) +
        form$spectral_density * tcrossprod(form$noise_effect)

      expect_equal(form$dimension, order + 1L)
      expect_lt(max(abs(lyapunov)) / max(abs(stationary)), 1e-12)

      representation <- .state_space_representation(kernel)
      covariance <- vapply(lags, function(lag) {
        drop(form$observation %*% representation$transition(lag) %*%
               stationary %*% form$observation)
      }, numeric(1))

      expect_lt(
        max(abs(covariance - drop(evaluate_kernel(kernel, matrix(lags), 0)))),
        1e-12 * parameters[[1L]]
      )

      # exp(F t) is a semigroup and solves dA / dt = F A.
      transition <- representation$transition
      expect_equal(transition(0), diag(order + 1L))
      expect_equal(transition(0.3) %*% transition(0.5), transition(0.8),
                   tolerance = 1e-12)
      step <- 1e-5
      derivative <- (transition(0.7 + step) - transition(0.7 - step)) /
        (2 * step)
      expect_equal(derivative, form$feedback %*% transition(0.7),
                   tolerance = 1e-7)
    }
  }
})


test_that("sums and scaled kernels stack block-diagonal states", {
  kernel <- state_space_kernels$scaled
  representation <- .state_space_representation(kernel)
  lags <- c(0, 0.1, 1, 4)
  covariance <- vapply(lags, function(lag) {
    drop(representation$observation %*% representation$transition(lag) %*%
           representation$stationary %*% representation$observation)
  }, numeric(1))

  expect_equal(representation$dimension, 3L)
  expect_lt(
    max(abs(covariance - drop(evaluate_kernel(kernel, matrix(lags), 0)))),
    1e-12
  )
  expect_equal(
    representation$feedback %*% representation$stationary +
      representation$stationary %*% t(representation$feedback) +
      representation$process_noise,
    matrix(0, 3, 3),
    tolerance = 1e-12
  )
})


test_that("Kalman smoothing reproduces the exact GP on irregular and tied inputs", {
  data <- state_space_data()

  for (name in names(state_space_kernels)) {
    kernel <- state_space_kernels[[name]]

    for (noise in list(0.09, data$noise)) {
      label <- sprintf("%s, %d noise value(s)", name, length(noise))
      exact <- fit_gp(data$x, data$y, kernel, noise_variance = noise,
                      mean = constant_mean(0.3))
      kalman <- .fit_state_space_model(data$x, data$y, kernel, noise,
                                       constant_mean(0.3))
      expected <- predict_gp(exact, data$query, include_covariance = TRUE,
                             observation_noise_variance = 0.05)
      actual <- .predict_state_space(kalman, data$query,
                                     include_covariance = TRUE,
                                     observation_noise_variance = 0.05)
      at_training <- predict_gp(exact, data$x, observation_noise_variance = 0)

      expect_identical(exact$jitter, 0, label = label)
      expect_lt(
        max_relative_error(kalman$log_marginal_likelihood,
                           log_marginal_likelihood(exact)),
        1e-8,
        label = paste(label, "log marginal likelihood")
      )
      expect_lt(max_relative_error(kalman$fitted_mean, at_training$mean),
                1e-8, label = paste(label, "fitted mean"))
      expect_lt(
        max_relative_error(kalman$fitted_variance, at_training$latent_variance),
        1e-8,
        label = paste(label, "fitted variance")
      )
      expect_lt(max_relative_error(actual$mean, expected$mean), 1e-8,
                label = paste(label, "mean"))
      expect_lt(
        max_relative_error(actual$latent_variance, expected$latent_variance),
        1e-8,
        label = paste(label, "variance")
      )
      expect_lt(
        max(abs(actual$latent_covariance - expected$latent_covariance)) /
          max(abs(expected$latent_covariance)),
        1e-8,
        label = paste(label, "covariance")
      )
      expect_equal(actual$prediction_interval, expected$prediction_interval,
                   tolerance = 1e-8)
      expect_named(actual, names(expected))
    }
  }
})


test_that("state-space inference matches at 500 points, fixed parametric means, and unsorted inputs", {
  data <- state_space_data(n = 500L, seed = 2L)
  kernel <- sum_kernel(matern52_kernel(length_scale = 2),
                       matern12_kernel(variance = 0.1, length_scale = 0.5))
  mean <- linear_mean(c(0.3, -0.05))
  exact <- fit_gp(data$x, data$y, kernel, noise_variance = 0.05, mean = mean)
  kalman <- .fit_state_space_model(data$x, data$y, kernel, 0.05, mean)

  expect_lt(
    max_relative_error(kalman$log_marginal_likelihood,
                       log_marginal_likelihood(exact)),
    1e-8
  )
  expect_lt(
    max_relative_error(.predict_state_space(kalman, data$query)$mean,
                       predict_gp(exact, data$query)$mean),
    1e-8
  )

  order <- sample(length(data$x))
  shuffled <- .fit_state_space_model(data$x[order], data$y[order], kernel,
                                     0.05, mean)

  expect_equal(shuffled$log_marginal_likelihood,
               kalman$log_marginal_likelihood, tolerance = 1e-12)
  expect_equal(shuffled$fitted_mean, kalman$fitted_mean[order],
               tolerance = 1e-10)
})


test_that("forward-filtering backward-sampling draws have the exact posterior moments", {
  data <- state_space_data(n = 40L, seed = 3L)
  kernel <- sum_kernel(matern32_kernel(length_scale = 1.5),
                       matern12_kernel(variance = 0.2, length_scale = 3))
  query <- c(-1, 2.5, data$x[[10]], data$x[[20]], 9.2, 9.25, 15, 22)
  exact <- fit_gp(data$x, data$y, kernel, noise_variance = 0.1,
                  mean = constant_mean(0.2))
  kalman <- .fit_state_space_model(data$x, data$y, kernel, 0.1,
                                   constant_mean(0.2))
  expected <- predict_gp(exact, query, include_covariance = TRUE,
                         observation_noise_variance = 0.1)
  latent <- sample_gp_posterior(kalman, query, n_draws = n_statistical_draws,
                                seed = 11)
  observation <- sample_gp_posterior(kalman, query,
                                     n_draws = n_statistical_draws,
                                     type = "observation", seed = 12)

  expect_s3_class(latent, "gaussianprocesses_samples")
  expect_identical(dim(latent$latent), c(length(query), n_statistical_draws))
  expect_equal(latent$mean, expected$mean, tolerance = 1e-10)
  expect_gaussian_moments(latent$latent, expected$mean,
                          expected$latent_covariance, "latent")
  expect_gaussian_moments(observation$observation, expected$mean,
                          expected$observation_covariance, "observation")
  expect_identical(
    sample_gp_posterior(kalman, query, n_draws = 3, seed = 11)$latent,
    sample_gp_posterior(kalman, query, n_draws = 3, seed = 11)$latent
  )
  expect_identical(dim(sample_gp_posterior(kalman, 1, seed = 1)$latent),
                   c(1L, 1L))
  expect_error(
    sample_gp_posterior(kalman, query, observation_noise_variance = 0.1),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("kernels without an exact state-space form raise classed errors", {
  x <- seq(0, 5, length.out = 10)
  y <- sin(x)
  unsupported <- list(
    rbf_kernel(),
    rational_quadratic_kernel(),
    periodic_kernel(),
    linear_kernel(),
    white_noise_kernel(),
    sum_kernel(matern32_kernel(), rbf_kernel()),
    product_kernel(matern32_kernel(), matern12_kernel()),
    changepoint_kernel(matern32_kernel(), matern12_kernel())
  )

  for (kernel in unsupported) {
    expect_error(
      .fit_state_space_model(x, y, kernel, 0.1),
      class = "gaussianprocesses_state_space_error"
    )
  }

  expect_error(
    fit_time_series_gp(x, y, rbf_kernel(), noise_variance = 0.1,
                       method = "state_space"),
    "no exact state-space representation",
    class = "gaussianprocesses_state_space_error"
  )
  expect_error(
    .fit_state_space_model(cbind(x, x), y, matern32_kernel(), 0.1),
    class = "gaussianprocesses_state_space_error"
  )
  expect_error(
    .state_space_representation(select_dimensions(matern32_kernel(), 2)),
    "must select column 1",
    class = "gaussianprocesses_state_space_error"
  )
  expect_error(
    .state_space_representation(matern32_kernel(length_scale = c(1, 2))),
    "single length scale",
    class = "gaussianprocesses_state_space_error"
  )
  expect_error(
    .fit_state_space_model(x, y, matern32_kernel(), 0),
    "positive noise",
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    .fit_state_space_model(x, y, matern32_kernel(), 0.1,
                           mean = constant_mean(estimate_coefficients())),
    "fixed coefficients",
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    fit_time_series_gp(x, y, matern32_kernel(), noise_variance = 0.1,
                       method = "state_space", initial_jitter = 1e-8),
    "'initial_jitter' apply only to method = \"exact\"",
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    fit_time_series_gp(x, y, matern32_kernel(), "auto", "state_space",
                       0.1, zero_mean(), 1),
    "unnamed",
    class = "gaussianprocesses_validation_error"
  )
})


test_that("state-space time-series models forecast like exact ones", {
  set.seed(4)
  time <- as.Date("2026-01-01") + sort(sample(0:400, 150))
  y <- sin(as.numeric(time) / 20) + stats::rnorm(150, sd = 0.2)
  kernel <- sum_kernel(matern52_kernel(length_scale = 30),
                       matern12_kernel(variance = 0.1, length_scale = 5))
  future <- as.Date("2025-12-29") + c(0, 13.5, 210, 455)
  exact <- fit_time_series_gp(time, y, kernel, noise_variance = 0.04)
  kalman <- fit_time_series_gp(time, y, kernel, noise_variance = 0.04,
                               method = "state_space")
  expected <- forecast_gp(exact, future, include_covariance = TRUE)
  actual <- forecast_gp(kalman, future, include_covariance = TRUE)

  expect_s3_class(kalman, "gaussianprocesses_state_space_model")
  expect_identical(kalman$time[c("origin", "unit", "input_type")],
                   exact$time[c("origin", "unit", "input_type")])
  expect_s3_class(actual, "gaussianprocesses_time_forecast")
  expect_identical(actual$region, expected$region)
  expect_identical(actual$time, expected$time)

  for (field in c("mean", "latent_variance", "observation_variance",
                  "latent_interval", "prediction_interval",
                  "latent_covariance")) {
    expect_equal(actual[[field]], expected[[field]], tolerance = 1e-8,
                 label = field)
  }

  expect_equal(predict(kalman)$mean, predict(exact)$mean, tolerance = 1e-8)
  expect_equal(predict(kalman, future)$mean, expected$mean, tolerance = 1e-8)
  expect_equal(as.numeric(logLik(kalman)), as.numeric(logLik(exact)),
               tolerance = 1e-10)
  expect_identical(attr(logLik(kalman), "df"), 0L)
  expect_equal(log_marginal_likelihood(kalman), as.numeric(logLik(exact)),
               tolerance = 1e-10)
  expect_identical(nobs(kalman), 150L)
  expect_equal(fitted(kalman), fitted(exact), tolerance = 1e-8)
  expect_equal(residuals(kalman), residuals(exact), tolerance = 1e-8)

  summary_output <- capture.output(print(summary(kalman)))
  expect_match(summary_output, "State-space", all = FALSE)
  expect_match(summary_output, "state dimension: 4", all = FALSE)
  expect_match(summary_output, "time unit: days", all = FALSE)
  expect_match(capture.output(print(kalman)), "Kalman smoothing",
               all = FALSE)

  heteroscedastic <- fit_time_series_gp(time, y, kernel,
                                        noise_variance = rep(0.04, 150),
                                        method = "state_space")
  expect_error(forecast_gp(heteroscedastic, future),
               class = "gaussianprocesses_validation_error")
  expect_equal(
    forecast_gp(heteroscedastic, future,
                observation_noise_variance = 0.04)$mean,
    expected$mean,
    tolerance = 1e-8
  )
})


test_that("state-space models without time metadata predict at numeric inputs", {
  data <- state_space_data(n = 60L, seed = 5L)
  kalman <- .fit_state_space_model(data$x, data$y, matern32_kernel(), 0.1)
  exact <- fit_gp(data$x, data$y, matern32_kernel(), noise_variance = 0.1)

  expect_equal(predict(kalman)$mean, predict(exact)$mean, tolerance = 1e-8)
  expect_equal(predict(kalman, c(1, 2))$mean, predict(exact, c(1, 2))$mean,
               tolerance = 1e-8)
  expect_equal(predict_gp(kalman, c(1, 2), include_covariance = TRUE),
               predict_gp(exact, c(1, 2), include_covariance = TRUE),
               tolerance = 1e-8)
  expect_error(predict_gp(kalman, 1, include_covariance = NA),
               class = "gaussianprocesses_validation_error")
  expect_error(log_marginal_likelihood_gradient(kalman),
               "finite differences",
               class = "gaussianprocesses_validation_error")
  expect_error(forecast_gp(kalman, 1),
               class = "gaussianprocesses_validation_error")
  expect_error(.prepare_state_space_model(exact),
               class = "gaussianprocesses_validation_error")
  expect_identical(kalman$schema_version, .model_schema_version)
})


test_that("rolling-origin evaluation and optimization accept the state-space method", {
  set.seed(6)
  time <- sort(stats::runif(80, 0, 40))
  y <- sin(time / 3) + stats::rnorm(80, sd = 0.2)
  kernel <- matern52_kernel(length_scale = 3)

  exact_backtest <- rolling_origin_gp(time, y, kernel, initial_window = 60,
                                      horizon = 4, step = 8,
                                      noise_variance = 0.04)
  kalman_backtest <- rolling_origin_gp(time, y, kernel, initial_window = 60,
                                       horizon = 4, step = 8,
                                       noise_variance = 0.04,
                                       method = "state_space")
  expect_equal(kalman_backtest$mean, exact_backtest$mean, tolerance = 1e-8)
  expect_equal(kalman_backtest$observation_variance,
               exact_backtest$observation_variance, tolerance = 1e-8)

  exact <- optimize_time_series_gp(time, y, kernel, noise_variance = 0.1,
                                   n_starts = 1)
  kalman <- optimize_time_series_gp(time, y, kernel, noise_variance = 0.1,
                                    n_starts = 1, method = "state_space")

  expect_s3_class(kalman, "gaussianprocesses_state_space_model")
  expect_identical(kalman$optimization$gradient, "numerical")
  expect_true(kalman$optimization$converged)
  expect_named(kalman$optimization$optimized_parameters,
               c("variance", "length_scale", "noise_variance"))
  expect_equal(kalman$log_marginal_likelihood,
               log_marginal_likelihood(exact), tolerance = 1e-6)
  expect_equal(kalman$optimization$optimized_parameters,
               exact$optimization$optimized_parameters, tolerance = 1e-3)
  expect_identical(attr(logLik(kalman), "df"), 3L)
  expect_false(is.null(kalman$time))
  expect_match(capture.output(print(kalman)), "optimization: converged",
               all = FALSE)

  # A known noise pattern has its scale estimated.
  scaled <- optimize_time_series_gp(time, y, kernel,
                                    noise_variance = rep(c(0.1, 0.2), 40),
                                    n_starts = 1, method = "state_space")
  expect_named(scaled$optimization$optimized_parameters,
               c("variance", "length_scale", "noise_scale"))
  fixed <- optimize_time_series_gp(time, y, kernel, noise_variance = 0.04,
                                   optimize_noise = FALSE, n_starts = 1,
                                   method = "state_space")
  expect_identical(fixed$noise_variance, 0.04)
  expect_error(
    optimize_time_series_gp(time, y, kernel, method = "state_space",
                            gradient = "analytical"),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    optimize_time_series_gp(time, y, kernel, method = "state_space",
                            optimize_noise = NA),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    optimize_time_series_gp(time, y, kernel, noise_variance = 0,
                            method = "state_space"),
    "positive noise",
    class = "gaussianprocesses_validation_error"
  )

  # A breakdown of the filter is reported, not propagated as NaN.
  expect_error(
    .kalman_filter(.state_space_representation(kernel), time[1:3], y[1:3],
                   c(0.1, NaN, 0.1), rep(TRUE, 3)),
    "innovation variance",
    class = "gaussianprocesses_numerical_error"
  )
})
