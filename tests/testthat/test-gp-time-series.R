test_that("numeric time index preserves order and origin explicitly", {
  time <- c(10, 11.5, 15)
  index <- gp_time_index(time)

  expect_s3_class(index, "gaussianprocesses_time_index")
  expect_equal(as.numeric(index), c(0, 1.5, 5))
  expect_equal(attr(index, "origin"), 10)
  expect_equal(attr(index, "unit"), "native")
  expect_equal(attr(index, "input_type"), "numeric")
})


test_that("Date and POSIXt time indices use explicit units", {
  dates <- as.Date("2026-01-01") + c(0, 2, 5)
  date_index <- gp_time_index(dates)

  expect_equal(as.numeric(date_index), c(0, 2, 5))
  expect_equal(attr(date_index, "unit"), "days")

  timestamps <- as.POSIXct(
    c(
      "2026-01-01 00:00:00",
      "2026-01-01 01:00:00",
      "2026-01-01 03:00:00"
    ),
    tz = "UTC"
  )
  hour_index <- gp_time_index(timestamps, unit = "hours")

  expect_equal(as.numeric(hour_index), c(0, 1, 3))
  expect_equal(attr(hour_index, "unit"), "hours")
})


test_that("time observations are never silently sorted", {
  expect_error(
    gp_time_index(c(0, 2, 1)),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    gp_time_index(c(0, 1, 1)),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("time-series kernel is explicitly trend plus periodic plus local", {
  kernel <- time_series_kernel(
    period = 12,
    trend_variance = 0.2,
    periodic_variance = 1.5,
    periodic_length_scale = 0.8,
    local_variance = 0.4,
    local_length_scale = 2
  )

  expect_equal(kernel$type, "sum")
  expect_length(kernel$children, 3L)
  expect_equal(kernel$children[[1]]$type, "linear")
  expect_equal(kernel$children[[2]]$type, "periodic")
  expect_equal(kernel$children[[3]]$type, "matern32")
})


test_that("time-aware fit stores reusable time metadata", {
  dates <- as.Date("2026-01-01") + 0:4

  model <- fit_time_series_gp(
    dates,
    c(0, 1, 0.5, 1.2, 0.9),
    kernel = rbf_kernel(length_scale = 2),
    noise_variance = 0.1
  )

  expect_equal(model$time$unit, "days")
  expect_equal(model$time$input_type, "Date")
  expect_equal(model$time$training_index, 0:4)
  expect_equal(model$time$training_time, dates)
})


test_that("forecast helper labels interpolation and extrapolation", {
  time <- 0:4
  y <- c(0, 1, 0.4, 1.1, 0.8)

  model <- fit_time_series_gp(
    time,
    y,
    kernel = rbf_kernel(length_scale = 1.2),
    noise_variance = 0.1
  )

  result <- forecast_gp(
    model,
    c(-1, 2.5, 5, 6)
  )

  expect_equal(
    result$region,
    c(
      "backward_extrapolation",
      "interpolation",
      "forecast",
      "forecast"
    )
  )
  expect_length(result$mean, 4L)
  expect_true(all(result$observation_variance >= 0))
})


test_that("forecast helper reuses Date origin and units", {
  dates <- as.Date("2026-01-01") + 0:4

  model <- fit_time_series_gp(
    dates,
    c(0, 1, 0.5, 1.2, 1),
    kernel = matern32_kernel(length_scale = 2),
    noise_variance = 0.1
  )

  future <- as.Date("2026-01-01") + c(2, 5, 7)
  result <- forecast_gp(model, future)

  expect_equal(result$time_index, c(2, 5, 7))
  expect_equal(
    result$region,
    c("interpolation", "forecast", "forecast")
  )
})


test_that("forecast helper supports observation-specific training noise", {
  time <- 0:7

  model <- fit_time_series_gp(
    time,
    sin(time / 2),
    kernel = rbf_kernel(length_scale = 2),
    noise_variance = seq(0.01, 0.08, by = 0.01)
  )

  # The prediction noise is unknown for such models and must be supplied.
  expect_error(
    forecast_gp(model, c(8, 9)),
    class = "gaussianprocesses_validation_error"
  )

  result <- forecast_gp(
    model,
    c(3.5, 9),
    observation_noise_variance = c(0.1, 0.2)
  )

  expect_equal(
    result$observation_variance - result$latent_variance,
    c(0.1, 0.2),
    tolerance = 1e-12
  )
  expect_equal(result$observation_noise_variance, c(0.1, 0.2))
  expect_equal(result$region, c("interpolation", "forecast"))
})


test_that("rolling-origin evaluation uses only past data", {
  time <- 1:8
  y <- c(0, 0.5, 1, 0.8, 0.3, 0, 0.4, 0.9)

  backtest <- rolling_origin_gp(
    time,
    y,
    kernel = rbf_kernel(
      variance = 1,
      length_scale = 1.5
    ),
    initial_window = 4,
    horizon = 2,
    step = 2,
    noise_variance = 0.1,
    interval_level = 0.9
  )

  expect_s3_class(backtest, "gaussianprocesses_backtest")
  expect_equal(
    unique(backtest$origin_index),
    c(4, 6)
  )
  expect_equal(
    backtest$target_index,
    c(5, 6, 7, 8)
  )
  expect_equal(
    backtest$horizon,
    c(1, 2, 1, 2)
  )
  expect_true(
    all(backtest$target_index > backtest$origin_index)
  )
})


test_that("forecast metrics combine point and uncertainty assessment", {
  time <- 1:9
  y <- c(0, 0.4, 0.8, 0.6, 0.2, 0.1, 0.5, 0.9, 0.7)

  backtest <- rolling_origin_gp(
    time,
    y,
    kernel = matern32_kernel(
      variance = 1,
      length_scale = 1.5
    ),
    initial_window = 5,
    horizon = 1,
    noise_variance = 0.1
  )

  metrics <- forecast_metrics(backtest)

  expect_equal(metrics$n_forecasts, nrow(backtest))
  expect_gte(metrics$rmse, 0)
  expect_gte(metrics$mae, 0)
  expect_gte(metrics$standardized_rmse, 0)
  expect_true(
    metrics$empirical_coverage >= 0 &&
      metrics$empirical_coverage <= 1
  )
  expect_gt(metrics$mean_interval_width, 0)
  expect_true(is.finite(metrics$mean_log_predictive_density))
})


test_that("time-series workflow validates type and split constraints", {
  model <- fit_time_series_gp(
    as.Date("2026-01-01") + 0:3,
    c(0, 1, 0, 1),
    kernel = rbf_kernel(),
    noise_variance = 0.1
  )

  expect_error(
    forecast_gp(model, c(4, 5)),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    rolling_origin_gp(
      1:5,
      c(0, 1, 0, 1, 0),
      kernel = rbf_kernel(),
      initial_window = 5
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    rolling_origin_gp(
      1:5,
      c(0, 1, 0, 1, 0),
      kernel = rbf_kernel(),
      initial_window = 4,
      horizon = 2
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    forecast_metrics(data.frame()),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("time forecast and backtest have concise print methods", {
  model <- fit_time_series_gp(
    1:5,
    c(0, 1, 0.5, 1, 0.7),
    kernel = rbf_kernel(),
    noise_variance = 0.1
  )

  forecast <- forecast_gp(model, c(3.5, 6))
  forecast_output <- capture.output(print(forecast))

  expect_true(
    any(grepl("time-series forecast", forecast_output, fixed = TRUE))
  )
  expect_true(
    any(grepl("interpolation points: 1", forecast_output, fixed = TRUE))
  )
  expect_true(
    any(grepl("forward forecast points: 1", forecast_output, fixed = TRUE))
  )

  backtest <- rolling_origin_gp(
    1:7,
    c(0, 1, 0, 1, 0, 1, 0),
    kernel = rbf_kernel(),
    initial_window = 4,
    horizon = 1,
    noise_variance = 0.1
  )
  backtest_output <- capture.output(print(backtest))

  expect_true(
    any(grepl("rolling-origin backtest", backtest_output, fixed = TRUE))
  )
})


test_that("time inputs and units are validated", {
  error <- "gaussianprocesses_validation_error"
  dates <- as.Date("2026-01-01") + 0:5
  times <- as.POSIXct("2026-01-01 00:00", tz = "UTC") + 3600 * 0:5

  expect_error(gp_time_index(list(1, 2)), class = error)
  expect_error(gp_time_index(matrix(1:4, 2)), class = error)
  expect_error(gp_time_index(dates, unit = "weeks"), class = error)
  expect_error(gp_time_index(dates, unit = c("days", "hours")), class = error)
  expect_error(gp_time_index(1:3, unit = "days"), "use unit = 'auto'", class = error)
  expect_error(gp_time_index(numeric(0)), "must not be empty", class = error)
  expect_error(gp_time_index(c(1, NA)), "must be finite", class = error)
  expect_error(gp_time_index(c(dates[1], NA)), "missing values", class = error)
  expect_error(gp_time_index(c(3, 1, 2)), class = error)
  expect_error(
    gaussianprocesses:::.time_unit_seconds("weeks"),
    "Unsupported time unit 'weeks'",
    class = error
  )

  expect_identical(attr(gp_time_index(times), "unit"), "seconds")
  expect_identical(attr(gp_time_index(dates), "unit"), "days")
})


test_that("forecast times must match the training time type", {
  error <- "gaussianprocesses_validation_error"
  y <- sin(0:5)
  numeric_model <- fit_time_series_gp(0:5, y, rbf_kernel(1, 2), noise_variance = 0.01)
  date_model <- fit_time_series_gp(
    as.Date("2026-01-01") + 0:5,
    y,
    rbf_kernel(1, 2),
    noise_variance = 0.01
  )
  time_model <- fit_time_series_gp(
    as.POSIXct("2026-01-01 00:00", tz = "UTC") + 3600 * 0:5,
    y,
    rbf_kernel(1, 2),
    unit = "hours",
    noise_variance = 0.01
  )

  expect_error(forecast_gp(numeric_model, "7"), class = error)
  expect_error(forecast_gp(numeric_model, c(6, Inf)), class = error)
  expect_error(forecast_gp(date_model, 6:7), class = error)
  expect_error(forecast_gp(time_model, 6:7), class = error)
  expect_error(forecast_gp(date_model, as.Date(NA)), class = error)
  expect_error(forecast_gp(numeric_model, numeric(0)), class = error)

  # A model fitted without time metadata cannot forecast.
  plain <- fit_gp(0:5, y, rbf_kernel(1, 2), noise_variance = 0.01)
  expect_error(forecast_gp(plain, 6), "no time metadata", class = error)

  expect_error(
    rolling_origin_gp(1:10, sin(1:9), rbf_kernel(), initial_window = 5),
    "same number of observations",
    class = error
  )
  expect_error(forecast_metrics(data.frame(observed = 1)), class = error)

  backtest <- rolling_origin_gp(1:12, sin(1:12), rbf_kernel(1, 3), 8, noise_variance = 0.01)
  broken <- backtest
  broken$observation_variance[1] <- 0
  expect_error(
    forecast_metrics(broken),
    "finite and strictly positive",
    class = "gaussianprocesses_numerical_error"
  )
})
