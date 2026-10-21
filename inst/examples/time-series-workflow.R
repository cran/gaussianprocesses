# Time-series Gaussian-process workflow
#
# This example keeps the model structure explicit:
# - a linear covariance component for nonstationary trend;
# - a periodic component;
# - a local Matérn component.
#
# The forecast helper labels interpolation and extrapolation separately.

time <- seq(0, 24, by = 1)
response <- 0.04 * time +
  sin(2 * pi * time / 12) +
  0.15 * cos(2 * pi * time / 3)

kernel <- time_series_kernel(
  period = 12,
  trend_variance = 0.05,
  periodic_variance = 1,
  periodic_length_scale = 0.8,
  local_variance = 0.2,
  local_length_scale = 2
)

model <- fit_time_series_gp(
  time,
  response,
  kernel = kernel,
  noise_variance = 0.05
)

# Includes both interpolation (inside the training range) and extrapolation.
requested_time <- c(6.5, 18.5, 25, 26, 27)

forecast <- forecast_gp(
  model,
  requested_time,
  interval_level = 0.95
)

data.frame(
  time = requested_time,
  region = forecast$region,
  mean = forecast$mean,
  lower = forecast$prediction_interval[, "lower"],
  upper = forecast$prediction_interval[, "upper"]
)

backtest <- rolling_origin_gp(
  time,
  response,
  kernel = kernel,
  initial_window = 12,
  horizon = 3,
  step = 3,
  noise_variance = 0.05
)

forecast_metrics(backtest)
