## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)
library(gaussianprocesses)

## -----------------------------------------------------------------------------
time <- as.Date("2026-01-01") + 0:9

index <- gp_time_index(time)

index
attributes(index)[c("origin", "unit", "input_type")]

## -----------------------------------------------------------------------------
kernel <- time_series_kernel(
  period = 12,
  trend_variance = 0.03,
  periodic_variance = 1,
  periodic_length_scale = 0.8,
  local_variance = 0.25,
  local_length_scale = 2
)

kernel

## -----------------------------------------------------------------------------
time_numeric <- 0:24

response <-
  0.04 * time_numeric +
  sin(2 * pi * time_numeric / 12) +
  0.15 * cos(2 * pi * time_numeric / 3)

model <- fit_time_series_gp(
  time_numeric,
  response,
  kernel = kernel,
  noise_variance = 0.05
)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
backtest <- rolling_origin_gp(
  time_numeric,
  response,
  kernel = kernel,
  initial_window = 12,
  horizon = 3,
  step = 3,
  noise_variance = 0.05,
  interval_level = 0.9
)

head(backtest)

## -----------------------------------------------------------------------------
forecast_metrics(backtest)

## -----------------------------------------------------------------------------
set.seed(28)
dates <- as.Date("2026-01-01") + 0:89
day <- as.numeric(gp_time_index(dates))
break_day <- as.numeric(as.Date("2026-03-15") - dates[1])

regime_signal <- function(day) {
  ifelse(
    day < break_day,
    sin(2 * pi * day / 30),
    2.5 + 0.4 * sin(2 * pi * day / 7)
  )
}
observed <- regime_signal(day) + rnorm(90, sd = 0.15)

## -----------------------------------------------------------------------------
regime <- function() {
  sum_kernel(
    rbf_kernel(variance = 1, length_scale = 60),
    matern52_kernel(variance = 0.5, length_scale = 5)
  )
}

single <- optimize_gp(
  day,
  observed,
  kernel = regime(),
  noise_variance = 0.05
)
changing <- optimize_gp(
  day,
  observed,
  kernel = changepoint_kernel(regime(), regime(), location = 45),
  noise_variance = 0.05,
  n_starts = 4
)

c(
  single = log_marginal_likelihood(single),
  changepoint = log_marginal_likelihood(changing)
)

## -----------------------------------------------------------------------------
dates[1] + changing$kernel$parameters$location

## -----------------------------------------------------------------------------
single_model <- fit_time_series_gp(
  dates,
  observed,
  kernel = single$kernel,
  noise_variance = single$noise_variance
)
changepoint_model <- fit_time_series_gp(
  dates,
  observed,
  kernel = changing$kernel,
  noise_variance = changing$noise_variance
)

future <- max(dates) + c(1, 7, 14, 21)
single_forecast <- forecast_gp(single_model, future)
changepoint_forecast <- forecast_gp(changepoint_model, future)

data.frame(
  date = future,
  truth = regime_signal(as.numeric(future - dates[1])),
  single = single_forecast$mean,
  single_sd = single_forecast$latent_sd,
  changepoint = changepoint_forecast$mean,
  changepoint_sd = changepoint_forecast$latent_sd
)

## ----break-forecasts, fig.alt = "Observed series from January to March 2026 with a jump on 15 March, followed by three weeks of forecasts. The single-kernel forecast (dashed) falls from the new level towards zero with a widening interval; the changepoint forecast (solid) stays near the new level of about 2.5 with a narrower interval."----
grid <- seq(min(dates), max(dates) + 21, by = 1)
single_path <- forecast_gp(single_model, grid)
changepoint_path <- forecast_gp(changepoint_model, grid)

plot(
  dates,
  observed,
  xlim = range(grid),
  ylim = range(
    single_path$latent_interval,
    changepoint_path$latent_interval,
    observed
  ),
  xlab = "date",
  ylab = "series"
)
abline(v = as.Date("2026-03-15"), col = "grey60", lty = 3)
lines(grid, changepoint_path$mean)
matlines(grid, changepoint_path$latent_interval, lty = 1, col = "grey40")
lines(grid, single_path$mean, lty = 2)
matlines(grid, single_path$latent_interval, lty = 2, col = "grey40")
legend(
  "topleft",
  legend = c("changepoint", "single kernel"),
  lty = c(1, 2),
  bty = "n"
)

## -----------------------------------------------------------------------------
set.seed(3)
long_time <- sort(runif(1500, 0, 300))
long_y <- sin(2 * pi * long_time / 50) + 0.3 * sin(long_time / 3) +
  rnorm(1500, sd = 0.2)
long_kernel <- sum_kernel(
  matern52_kernel(variance = 1, length_scale = 15),
  matern12_kernel(variance = 0.1, length_scale = 2)
)

exact_long <- fit_time_series_gp(
  long_time,
  long_y,
  kernel = long_kernel,
  noise_variance = 0.04
)
kalman_long <- fit_time_series_gp(
  long_time,
  long_y,
  kernel = long_kernel,
  noise_variance = 0.04,
  method = "state_space"
)

c(exact = logLik(exact_long), state_space = logLik(kalman_long))

## -----------------------------------------------------------------------------
requested_long <- c(-5, 150.25, 299.9, 305, 330)
exact_forecast <- forecast_gp(exact_long, requested_long)
kalman_forecast <- forecast_gp(kalman_long, requested_long)

data.frame(
  time = requested_long,
  region = kalman_forecast$region,
  exact_mean = exact_forecast$mean,
  state_space_mean = kalman_forecast$mean,
  exact_sd = exact_forecast$observation_sd,
  state_space_sd = kalman_forecast$observation_sd
)

c(
  mean = max(abs(kalman_forecast$mean - exact_forecast$mean)),
  sd = max(abs(kalman_forecast$observation_sd - exact_forecast$observation_sd))
)

## -----------------------------------------------------------------------------
tryCatch(
  fit_time_series_gp(
    time_numeric,
    response,
    kernel = kernel,
    noise_variance = 0.05,
    method = "state_space"
  ),
  gaussianprocesses_state_space_error = conditionMessage
)

## -----------------------------------------------------------------------------
first <- seq_len(300)
kalman_fit <- optimize_time_series_gp(
  long_time[first],
  long_y[first],
  kernel = matern52_kernel(length_scale = 10),
  noise_variance = 0.1,
  n_starts = 1,
  method = "state_space"
)
exact_fit <- optimize_time_series_gp(
  long_time[first],
  long_y[first],
  kernel = matern52_kernel(length_scale = 10),
  noise_variance = 0.1,
  n_starts = 1
)

rbind(
  exact = exact_fit$optimization$optimized_parameters,
  state_space = kalman_fit$optimization$optimized_parameters
)

