test_that("posterior prediction matches direct Gaussian conditioning", {
  x <- c(0, 1, 2)
  y <- c(1, -0.5, 0.25)
  x_new <- c(0.5, 1.5)
  kernel <- rbf_kernel(variance = 1.4, length_scale = 0.8)
  noise_variance <- 0.2

  model <- fit_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = noise_variance
  )
  prediction <- predict_gp(
    model,
    x_new,
    include_covariance = TRUE
  )

  training_covariance <- evaluate_kernel(kernel, x) +
    diag(noise_variance, length(x))
  cross_covariance <- evaluate_kernel(kernel, x, x_new)
  test_covariance <- evaluate_kernel(kernel, x_new)

  expected_mean <- drop(
    crossprod(
      cross_covariance,
      solve(training_covariance, y)
    )
  )
  expected_covariance <- test_covariance -
    crossprod(
      cross_covariance,
      solve(training_covariance, cross_covariance)
    )

  expect_equal(prediction$mean, expected_mean, tolerance = 1e-10)
  expect_equal(
    prediction$latent_covariance,
    expected_covariance,
    tolerance = 1e-10
  )
  expect_equal(
    prediction$latent_variance,
    diag(expected_covariance),
    tolerance = 1e-10
  )
})


test_that("constant prior mean is restored in posterior predictions", {
  x <- c(0, 1, 2)
  y <- c(3, 4, 3.5)
  x_new <- c(0.5, 1.5)
  kernel <- matern32_kernel(variance = 1, length_scale = 1.2)
  mean <- constant_mean(3)

  model <- fit_gp(
    x,
    y,
    kernel = kernel,
    mean = mean,
    noise_variance = 0.1
  )
  prediction <- predict_gp(model, x_new)

  covariance <- evaluate_kernel(kernel, x) + diag(0.1, 3)
  cross_covariance <- evaluate_kernel(kernel, x, x_new)
  expected <- 3 + drop(
    crossprod(
      cross_covariance,
      solve(covariance, y - 3)
    )
  )

  expect_equal(prediction$mean, expected, tolerance = 1e-10)
})


test_that("latent and observation uncertainty remain distinct", {
  model <- fit_gp(
    c(0, 1, 2),
    c(1, 0, 1),
    kernel = rbf_kernel(length_scale = 0.7),
    noise_variance = 0.25
  )

  prediction <- predict_gp(model, c(0.5, 1.5))

  expect_equal(
    prediction$observation_variance - prediction$latent_variance,
    rep(0.25, 2),
    tolerance = 1e-12
  )
  expect_true(all(prediction$observation_sd >= prediction$latent_sd))
})


test_that("credible and prediction intervals use the requested level", {
  model <- fit_gp(
    c(0, 1, 2),
    c(0, 1, 0),
    kernel = rbf_kernel(),
    noise_variance = 0.1
  )

  prediction <- predict_gp(
    model,
    0.5,
    interval_level = 0.8
  )

  quantile <- stats::qnorm(0.9)
  latent_half_width <- quantile * prediction$latent_sd
  observation_half_width <- quantile * prediction$observation_sd

  # A single prediction point makes `[, "upper"]` keep its column name.
  expect_equal(
    unname(prediction$latent_interval[, "upper"] - prediction$mean),
    latent_half_width,
    tolerance = 1e-12
  )
  expect_equal(
    unname(prediction$prediction_interval[, "upper"] - prediction$mean),
    observation_half_width,
    tolerance = 1e-12
  )
  expect_true(
    all(
      prediction$prediction_interval[, "upper"] -
        prediction$prediction_interval[, "lower"] >=
        prediction$latent_interval[, "upper"] -
        prediction$latent_interval[, "lower"]
    )
  )
})


test_that("noise-free GP interpolates unique training observations", {
  x <- c(-1, 0, 1)
  y <- c(1.5, -0.25, 0.75)

  model <- fit_gp(
    x,
    y,
    kernel = rbf_kernel(variance = 2, length_scale = 0.8),
    noise_variance = 0
  )
  prediction <- predict_gp(model, x)

  expect_equal(prediction$mean, y, tolerance = 1e-9)
  expect_equal(
    prediction$latent_variance,
    rep(0, length(x)),
    tolerance = 1e-9
  )
})


test_that("full posterior covariance is returned only when requested", {
  model <- fit_gp(
    c(0, 1, 2),
    c(1, 0, 1),
    kernel = matern52_kernel(),
    noise_variance = 0.05
  )

  marginal <- predict_gp(model, c(0.25, 0.75))
  full <- predict_gp(
    model,
    c(0.25, 0.75),
    include_covariance = TRUE
  )

  expect_false("latent_covariance" %in% names(marginal))
  expect_false("observation_covariance" %in% names(marginal))

  expect_equal(dim(full$latent_covariance), c(2L, 2L))
  expect_equal(dim(full$observation_covariance), c(2L, 2L))
  expect_equal(
    diag(full$observation_covariance) -
      diag(full$latent_covariance),
    rep(0.05, 2),
    tolerance = 1e-12
  )
})


test_that("small negative variance is clipped but material failure is rejected", {
  expect_equal(
    gaussianprocesses:::.sanitize_posterior_variance(
      c(1, -1e-12),
      reference_scale = 1,
      tolerance = 1e-10
    ),
    c(1, 0)
  )

  expect_error(
    gaussianprocesses:::.sanitize_posterior_variance(
      c(1, -1e-3),
      reference_scale = 1,
      tolerance = 1e-10
    ),
    class = "gaussianprocesses_numerical_error"
  )
})


test_that("prediction inputs and options are validated", {
  x <- matrix(
    c(
      0, 0,
      1, 1,
      2, 2
    ),
    ncol = 2,
    byrow = TRUE
  )
  model <- fit_gp(
    x,
    c(0, 1, 2),
    kernel = rbf_kernel(),
    noise_variance = 0.1
  )

  expect_error(
    predict_gp(model, c(0, 1)),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    predict_gp(model, matrix(c(0, 1), nrow = 1), interval_level = 1),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    predict_gp(model, matrix(c(0, 1), nrow = 1), include_covariance = NA),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    predict_gp("not a model", matrix(c(0, 1), nrow = 1)),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("prediction objects have a concise print method", {
  model <- fit_gp(
    c(0, 1),
    c(1, 2),
    kernel = rbf_kernel(),
    noise_variance = 0.1
  )

  prediction <- predict_gp(model, c(0.25, 0.75))
  output <- capture.output(print(prediction))

  expect_true(any(grepl("posterior prediction", output, fixed = TRUE)))
  expect_true(any(grepl("prediction points: 2", output, fixed = TRUE)))
  expect_true(any(grepl("interval level: 95%", output, fixed = TRUE)))
})


# Marginal prediction in blocks -------------------------------------------------

# Largest absolute difference relative to the largest reference magnitude.
relative_difference <- function(actual, expected) {
  max(abs(actual - expected)) / max(abs(expected), .Machine$double.xmin)
}


blocked_prediction_cases <- function() {
  x_1d <- seq(-2, 2, length.out = 18)
  x_2d <- cbind(seq(-1.5, 1.5, length.out = 20), cos(1:20))

  list(
    "1D RBF" = list(
      model = fit_gp(x_1d, sin(2 * x_1d), rbf_kernel(1.2, 0.6), 0.01),
      x = c(seq(-3, 3, length.out = 21), x_1d[1:2])
    ),
    "1D Matern 3/2" = list(
      model = fit_gp(x_1d, cos(x_1d), matern32_kernel(0.8, 0.9), 0.02),
      x = seq(-2.5, 2.5, length.out = 23)
    ),
    "2D ARD" = list(
      model = fit_gp(
        x_2d,
        sin(x_2d[, 1]) + x_2d[, 2],
        rbf_kernel(1, c(0.7, 1.8)),
        0.01
      ),
      x = cbind(seq(-2, 2, length.out = 23), sin(1:23))
    ),
    "composite" = list(
      model = fit_gp(
        x_2d,
        sin(2 * x_2d[, 1]),
        sum_kernel(
          product_kernel(
            matern52_kernel(1, c(0.8, 1.2)),
            periodic_kernel(0.7, 1, period = 1.5)
          ),
          linear_kernel(0.1)
        ),
        0.05
      ),
      x = cbind(seq(-2, 2, length.out = 23), cos(2:24))
    ),
    "known heteroscedastic noise" = list(
      model = fit_gp(
        x_1d,
        sin(x_1d),
        rbf_kernel(1, 0.8),
        noise_variance = seq(0.01, 0.1, length.out = length(x_1d))
      ),
      x = seq(-2.5, 2.5, length.out = 23),
      observation_noise_variance = seq(0.02, 0.2, length.out = 23)
    ),
    "constant mean" = list(
      model = fit_gp(
        x_1d,
        2 + sin(x_1d),
        rbf_kernel(1, 0.7),
        0.01,
        mean = constant_mean(2)
      ),
      x = seq(-3, 3, length.out = 23)
    )
  )
}


test_that("marginal prediction in blocks matches the full-covariance path", {
  fields <- c(
    "mean", "latent_variance", "observation_variance", "latent_sd",
    "observation_sd", "latent_interval", "prediction_interval"
  )

  for (name in names(blocked_prediction_cases())) {
    case <- blocked_prediction_cases()[[name]]
    m <- NROW(case$x)
    full <- predict_gp(
      case$model,
      case$x,
      include_covariance = TRUE,
      observation_noise_variance = case$observation_noise_variance
    )

    for (block_size in c(1L, 7L, m)) {
      marginal <- predict_gp(
        case$model,
        case$x,
        observation_noise_variance = case$observation_noise_variance,
        block_size = block_size
      )

      expect_null(marginal$latent_covariance)
      expect_identical(
        marginal$observation_noise_variance,
        full$observation_noise_variance
      )

      for (field in fields) {
        expect_lte(
          relative_difference(marginal[[field]], full[[field]]),
          1e-12,
          label = paste(name, "block", block_size, field)
        )
      }
    }
  }
})


# Records the row counts of every evaluate_kernel() call made by `code`.
record_kernel_evaluations <- function(code) {
  original <- gaussianprocesses::evaluate_kernel
  calls <- list()

  local_mocked_bindings(
    evaluate_kernel = function(kernel, x, y = NULL) {
      calls[[length(calls) + 1L]] <<- c(
        x_rows = NROW(x),
        y_rows = if (is.null(y)) NA_integer_ else NROW(y)
      )
      original(kernel, x, y)
    }
  )
  force(code)

  do.call(rbind, calls)
}


test_that("marginal prediction never forms the prediction-input covariance", {
  x <- seq(-2, 2, length.out = 20)
  model <- fit_gp(x, sin(x), rbf_kernel(), noise_variance = 0.01)
  x_new <- seq(-3, 3, length.out = 50)

  blocked <- record_kernel_evaluations(
    predict_gp(model, x_new, block_size = 7)
  )

  # Only cross-covariances between the 20 training inputs and blocks of at
  # most 7 prediction inputs: ceiling(50 / 7) = 8 blocks.
  expect_identical(nrow(blocked), 8L)
  expect_true(all(blocked[, "x_rows"] == 20L))
  expect_false(anyNA(blocked[, "y_rows"]))
  expect_true(all(blocked[, "y_rows"] <= 7L))
  expect_identical(sum(blocked[, "y_rows"]), 50L)

  # The opt-in full covariance does evaluate the 50 x 50 prior covariance.
  full <- record_kernel_evaluations(
    predict_gp(model, x_new, include_covariance = TRUE)
  )
  expect_true(any(full[, "x_rows"] == 50L & is.na(full[, "y_rows"])))
})


test_that("prediction wrappers use the marginal path", {
  simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 30)
  heteroscedastic <- fit_heteroscedastic_gp(
    simulation$x,
    simulation$observed,
    kernel = matern52_kernel(length_scale = 0.2),
    noise_kernel = rbf_kernel(length_scale = 0.4),
    max_iterations = 30
  )
  x_new <- seq(0, 1, length.out = 40)

  calls <- record_kernel_evaluations(
    predict_heteroscedastic_gp(heteroscedastic, x_new)
  )
  expect_false(anyNA(calls[, "y_rows"]))

  time_model <- fit_time_series_gp(
    1:24,
    sin((1:24) / 3),
    kernel = rbf_kernel(length_scale = 4),
    noise_variance = 0.01
  )
  calls <- record_kernel_evaluations(forecast_gp(time_model, 25:60))
  expect_false(anyNA(calls[, "y_rows"]))

  # Sampling needs the joint covariance, so it still requests it.
  calls <- record_kernel_evaluations(
    sample_gp_posterior(time_model, 25:30, n_draws = 2, seed = 1)
  )
  expect_true(any(is.na(calls[, "y_rows"])))
})


test_that("the block size is validated and chosen from a memory budget", {
  model <- fit_gp(c(0, 1, 2), c(0, 1, 0), rbf_kernel(), noise_variance = 0.1)

  for (block_size in list(0, -1, 2.5, "3", c(2, 3), NA)) {
    expect_error(
      predict_gp(model, c(0.5, 1.5), block_size = block_size),
      class = "gaussianprocesses_validation_error"
    )
  }

  # Blocks larger than the number of prediction inputs are allowed.
  expect_equal(
    predict_gp(model, c(0.5, 1.5), block_size = 100)$latent_variance,
    predict_gp(model, c(0.5, 1.5), block_size = 1)$latent_variance,
    tolerance = 1e-12
  )

  block_size <- gaussianprocesses:::.prediction_block_size
  expect_identical(block_size(NULL, n_training = 500, n_predictions = 1e5), 2097L)
  expect_identical(block_size(NULL, n_training = 500, n_predictions = 100), 100L)
  expect_identical(block_size(NULL, n_training = 1e7, n_predictions = 10), 1L)
  expect_identical(block_size(5, n_training = 500, n_predictions = 100), 5L)
})
