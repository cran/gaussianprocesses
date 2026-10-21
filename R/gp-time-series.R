# Time-series Gaussian-process workflows --------------------------------------

.time_unit_seconds <- function(unit) {
  switch(
    unit,
    seconds = 1,
    minutes = 60,
    hours = 3600,
    days = 86400,
    .gp_abort(
      sprintf("Unsupported time unit '%s'.", unit),
      "gaussianprocesses_validation_error"
    )
  )
}


.time_input_type <- function(time) {
  if (inherits(time, "Date")) {
    return("Date")
  }

  if (inherits(time, "POSIXt")) {
    return("POSIXt")
  }

  if (is.numeric(time) && is.null(dim(time))) {
    return("numeric")
  }

  .gp_abort(
    "'time' must be numeric, Date, or POSIXt.",
    "gaussianprocesses_validation_error"
  )
}


.resolve_time_unit <- function(type, unit) {
  allowed <- c("auto", "seconds", "minutes", "hours", "days")

  if (!is.character(unit) || length(unit) != 1L || !(unit %in% allowed)) {
    .gp_abort(
      paste0(
        "'unit' must be one of: ",
        paste(allowed, collapse = ", "),
        "."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (identical(type, "numeric")) {
    if (!identical(unit, "auto")) {
      .gp_abort(
        paste0(
          "Numeric time values already define their own unit; ",
          "use unit = 'auto'."
        ),
        "gaussianprocesses_validation_error"
      )
    }

    return("native")
  }

  if (!identical(unit, "auto")) {
    return(unit)
  }

  if (identical(type, "Date")) "days" else "seconds"
}


.time_to_elapsed <- function(time, origin, type, unit) {
  if (identical(type, "numeric")) {
    if (!is.numeric(time) || !is.null(dim(time))) {
      .gp_abort(
        "Future time values must be numeric to match the fitted model.",
        "gaussianprocesses_validation_error"
      )
    }

    if (any(!is.finite(time))) {
      .gp_abort(
        "Time values must be finite.",
        "gaussianprocesses_validation_error"
      )
    }

    return(as.numeric(time) - as.numeric(origin))
  }

  if (identical(type, "Date") && !inherits(time, "Date")) {
    .gp_abort(
      "Future time values must be Date objects to match the fitted model.",
      "gaussianprocesses_validation_error"
    )
  }

  if (identical(type, "POSIXt") && !inherits(time, "POSIXt")) {
    .gp_abort(
      "Future time values must be POSIXt objects to match the fitted model.",
      "gaussianprocesses_validation_error"
    )
  }

  if (any(is.na(time))) {
    .gp_abort(
      "Time values must not contain missing values.",
      "gaussianprocesses_validation_error"
    )
  }

  elapsed_seconds <- as.numeric(
    difftime(time, origin, units = "secs")
  )

  elapsed_seconds / .time_unit_seconds(unit)
}


.assert_strictly_increasing_time <- function(values, name = "time") {
  if (length(values) == 0L) {
    .gp_abort(
      sprintf("'%s' must not be empty.", name),
      "gaussianprocesses_validation_error"
    )
  }

  if (length(values) > 1L && any(diff(values) <= 0)) {
    .gp_abort(
      sprintf(
        "'%s' must be strictly increasing; observations are never sorted silently.",
        name
      ),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(values)
}


#' Build an explicit Gaussian-process time index
#'
#' Converts numeric, Date, or POSIXt values to elapsed numeric time while
#' retaining origin, unit, and input-class metadata. Input order is never
#' changed: non-increasing time values are rejected.
#'
#' @param time Numeric, Date, or POSIXt vector in strictly increasing order.
#' @param unit Time unit for Date or POSIXt inputs. `"auto"` uses days for
#'   Date values, seconds for POSIXt values, and the native scale for numeric
#'   inputs.
#'
#' @details
#' Time values are never sorted silently. Non-increasing input is rejected so
#' the temporal order supplied by the user remains explicit.
#'
#' @return A numeric vector of elapsed time with class
#'   `gaussianprocesses_time_index` and time metadata in attributes.
#' @examples
#' gp_time_index(as.Date("2026-01-01") + c(0, 1, 7))
#'
#' gp_time_index(
#'   as.POSIXct("2026-01-01 00:00", tz = "UTC") + c(0, 3600, 7200),
#'   unit = "hours"
#' )
#' @template stable
#' @export
gp_time_index <- function(time, unit = "auto") {
  type <- .time_input_type(time)
  resolved_unit <- .resolve_time_unit(type, unit)

  if (length(time) == 0L) {
    .gp_abort(
      "'time' must not be empty.",
      "gaussianprocesses_validation_error"
    )
  }

  if (identical(type, "numeric") && any(!is.finite(time))) {
    .gp_abort(
      "Numeric time values must be finite.",
      "gaussianprocesses_validation_error"
    )
  }

  if (!identical(type, "numeric") && any(is.na(time))) {
    .gp_abort(
      "Time values must not contain missing values.",
      "gaussianprocesses_validation_error"
    )
  }

  origin <- time[[1L]]
  values <- .time_to_elapsed(
    time,
    origin = origin,
    type = type,
    unit = resolved_unit
  )

  .assert_strictly_increasing_time(values)

  structure(
    values,
    origin = origin,
    unit = resolved_unit,
    input_type = type,
    class = c(
      "gaussianprocesses_time_index",
      "numeric"
    )
  )
}


#' Create a trend-periodic-local time-series kernel
#'
#' Constructs a transparent additive example consisting of a nonstationary
#' linear covariance component, a periodic covariance component, and a local
#' Matérn-3/2 component. The components are explicit so stationarity is never
#' implied for the full model.
#'
#' @param period Positive period in the same units as the time index.
#' @param trend_variance Positive variance of the linear trend component.
#' @param periodic_variance Positive variance of the periodic component.
#' @param periodic_length_scale Positive periodic length scale.
#' @param local_variance Positive variance of the local Matérn component.
#' @param local_length_scale Positive local Matérn length scale.
#' @param changepoint Optional location of a structural break, on the time
#'   index scale of [gp_time_index()]: elapsed time since the first training
#'   time, in its unit.
#' @param changepoint_steepness Positive steepness of the break.
#'
#' @details
#' With a `changepoint`, the kernel is a [changepoint_kernel()] between two
#' independent copies of the trend-periodic-local kernel, so every component
#' can change at the break; the copies start from the same parameter values
#' and are estimated separately. To place a break at a date, convert it with
#' the training times: `gp_time_index(c(time[1], date), unit)[2]`.
#'
#' @return A composite `gaussianprocesses_kernel`.
#' @examples
#' kernel <- time_series_kernel(period = 12)
#' kernel
#'
#' # A structural break 30 time units after the first observation.
#' names(kernel_parameters(
#'   time_series_kernel(period = 12, changepoint = 30),
#'   flatten = TRUE
#' ))
#' @templateVar experimental_parts The `changepoint` argument, which builds a [changepoint_kernel()], is experimental.
#' @template stable
#' @export
time_series_kernel <- function(
    period,
    trend_variance = 1,
    periodic_variance = 1,
    periodic_length_scale = 1,
    local_variance = 1,
    local_length_scale = 1,
    changepoint = NULL,
    changepoint_steepness = 1) {
  .assert_numeric_scalar(period, "period", lower = 0, lower_open = TRUE)

  kernel <- sum_kernel(
    linear_kernel(variance = trend_variance),
    periodic_kernel(
      variance = periodic_variance,
      length_scale = periodic_length_scale,
      period = period
    ),
    matern32_kernel(
      variance = local_variance,
      length_scale = local_length_scale
    )
  )

  if (is.null(changepoint)) {
    return(kernel)
  }

  changepoint_kernel(
    kernel,
    kernel,
    location = changepoint,
    steepness = changepoint_steepness
  )
}


#' Fit a GP to ordered time-series observations
#'
#' @param time Numeric, Date, or POSIXt vector in strictly increasing order.
#' @param y Numeric response vector.
#' @param kernel Gaussian-process kernel specification.
#' @param unit Time unit passed to [gp_time_index()].
#' @param method `"exact"` conditions a Gaussian prior on the observations
#'   with [fit_gp()]. `"state_space"` computes the same posterior and log
#'   marginal likelihood in \eqn{O(n)} time by Kalman filtering and
#'   smoothing; it needs a Matérn-1/2, 3/2, or 5/2 kernel, a sum of such
#'   kernels, or scaled versions of them (see State-space inference).
#' @param ... Additional arguments passed to [fit_gp()]. With
#'   `method = "state_space"` only `noise_variance` (positive) and `mean`
#'   (with fixed coefficients) apply.
#'
#' @return A fitted model containing explicit time-index metadata used by
#'   [forecast_gp()]: a `gaussianprocesses_model` for `method = "exact"` and
#'   a `gaussianprocesses_state_space_model` for `method = "state_space"`.
#'
#' @details
#' This function was called `fit_gp_time()` in version 0.1.0.
#'
#' @section State-space inference:
#' On one input, a Matérn process with \eqn{\nu = p + 1/2} is the first
#' component of the state \eqn{x = (f, f', \ldots, f^{(p)})} of the linear
#' stochastic differential equation \eqn{dx = F x\,dt + L\,dW}, with
#' \eqn{F} the companion matrix of \eqn{(s + \lambda)^{p+1}},
#' \eqn{\lambda = \sqrt{2\nu}/\ell}, and white noise \eqn{W} of spectral
#' density \eqn{q_c} (Hartikainen and Särkkä, 2010). The stationary state
#' covariance \eqn{P_\infty} solves
#' \eqn{F P_\infty + P_\infty F^\top + L q_c L^\top = 0}, and
#' \eqn{k(\tau) = H e^{F\tau} P_\infty H^\top} for \eqn{\tau \ge 0}
#' with \eqn{H = (1, 0, \ldots, 0)}. Sums of kernels stack their states, and
#' scaled kernels scale \eqn{P_\infty} and \eqn{q_c}. Between consecutive
#' times the state moves by \eqn{A = e^{F\Delta t}}, which has a closed
#' form because \eqn{F + \lambda I} is nilpotent, with process noise
#' \eqn{P_\infty - A P_\infty A^\top}.
#'
#' The Kalman filter (with Joseph-form updates) gives the log marginal
#' likelihood by the prediction-error decomposition, the Rauch-Tung-Striebel
#' smoother gives the posterior at the training times, and [forecast_gp()]
#' interleaves the requested times as unobserved steps. The results equal
#' those of the exact GP up to rounding, at \eqn{O(n)} cost instead of
#' \eqn{O(n^3)}, and no numerical jitter is used. [sample_gp_posterior()]
#' draws by forward filtering, backward sampling. Kernels without an exact
#' finite-dimensional state (RBF, rational quadratic, periodic, linear, and
#' products) raise an error of class `gaussianprocesses_state_space_error`.
#'
#' @references
#' Hartikainen, J., and Särkkä, S. (2010). Kalman filtering and smoothing
#' solutions to temporal Gaussian process regression models. *IEEE
#' International Workshop on Machine Learning for Signal Processing*,
#' 379--384.
#'
#' Särkkä, S., and Solin, A. (2019). *Applied Stochastic Differential
#' Equations*. Cambridge University Press.
#' @seealso [optimize_time_series_gp()] to estimate the hyperparameters.
#' @examples
#' time <- as.Date("2026-01-01") + 0:23
#' y <- sin(2 * pi * (0:23) / 12)
#'
#' model <- fit_time_series_gp(
#'   time,
#'   y,
#'   kernel = time_series_kernel(period = 12),
#'   noise_variance = 0.05
#' )
#' model$time[c("origin", "unit")]
#'
#' # The same posterior by Kalman smoothing, for a Matern kernel.
#' kernel <- matern32_kernel(length_scale = 3)
#' exact <- fit_time_series_gp(time, y, kernel, noise_variance = 0.05)
#' kalman <- fit_time_series_gp(time, y, kernel, noise_variance = 0.05,
#'   method = "state_space")
#' c(logLik(exact), logLik(kalman))
#' @templateVar experimental_parts `method = "state_space"` is experimental.
#' @template stable
#' @export
fit_time_series_gp <- function(
    time,
    y,
    kernel,
    unit = "auto",
    method = c("exact", "state_space"),
    ...) {
  method <- match.arg(method)
  index <- gp_time_index(time, unit = unit)

  model <- if (identical(method, "exact")) {
    fit_gp(as.numeric(index), y, kernel = kernel, ...)
  } else {
    .fit_state_space_gp(as.numeric(index), y, kernel, ...)
  }

  .with_time_metadata(model, time, index)
}


#' Optimize the hyperparameters of a time-series GP
#'
#' Estimates the kernel parameters and the noise of a GP on ordered
#' time-series observations by maximizing the log marginal likelihood, with
#' the optimizer of [optimize_gp()].
#'
#' @inheritParams fit_time_series_gp
#' @param kernel Initial kernel specification.
#' @param ... Additional arguments passed to [optimize_gp()]. With
#'   `method = "state_space"`, the arguments `noise_variance`, `mean`,
#'   `optimize_noise`, `n_starts`, `start_spread`, `lower`, `upper`, and
#'   `control` apply, with the meanings and defaults of [optimize_gp()].
#'
#' @details
#' With `method = "state_space"` each objective value is one Kalman filter
#' pass, in \eqn{O(n)} time ([fit_time_series_gp()] describes the method).
#' Gradients are finite-difference approximations by L-BFGS-B
#' (`gradient = "numerical"` in [optimize_gp()]); analytical sensitivity
#' recursions are not implemented yet.
#'
#' @return A fitted model as returned by [fit_time_series_gp()], with
#'   optimizer diagnostics in `$optimization`.
#' @examples
#' set.seed(1)
#' time <- sort(runif(200, 0, 50))
#' y <- sin(time / 3) + rnorm(200, sd = 0.2)
#'
#' model <- optimize_time_series_gp(time, y, matern52_kernel(),
#'   noise_variance = 0.1, n_starts = 1, method = "state_space")
#' model$optimization$optimized_parameters
#' @template experimental
#' @export
optimize_time_series_gp <- function(
    time,
    y,
    kernel,
    unit = "auto",
    method = c("exact", "state_space"),
    ...) {
  method <- match.arg(method)
  index <- gp_time_index(time, unit = unit)

  model <- if (identical(method, "exact")) {
    optimize_gp(as.numeric(index), y, kernel = kernel, ...)
  } else {
    .optimize_state_space_gp(as.numeric(index), y, kernel, ...)
  }

  .with_time_metadata(model, time, index)
}


.with_time_metadata <- function(model, time, index) {
  model$time <- list(
    origin = attr(index, "origin"),
    unit = attr(index, "unit"),
    input_type = attr(index, "input_type"),
    training_time = time,
    training_index = as.numeric(index)
  )

  model
}


.prepare_time_model <- function(model) {
  model <- if (inherits(model, "gaussianprocesses_state_space_model")) {
    .prepare_state_space_model(model)
  } else {
    .prepare_gp_model(model)
  }

  if (is.null(model$time)) {
    .gp_abort(
      paste0(
        "'model' has no time metadata. Fit the model with fit_time_series_gp() ",
        "before using time-series forecast helpers."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  model
}


.convert_future_time <- function(model, future_time) {
  model <- .prepare_time_model(model)

  values <- .time_to_elapsed(
    future_time,
    origin = model$time$origin,
    type = model$time$input_type,
    unit = model$time$unit
  )

  .assert_strictly_increasing_time(values, name = "future_time")
  values
}


#' Forecast or interpolate from a time-aware GP model
#'
#' Converts future time values using the training model's original time origin
#' and units, then evaluates the GP posterior: exactly, or by Kalman
#' smoothing for state-space models. Each requested point is labelled as
#' interpolation, backward extrapolation, or forward forecast.
#'
#' @param model A model fitted by [fit_time_series_gp()] or
#'   [optimize_time_series_gp()].
#' @param future_time Numeric, Date, or POSIXt values matching the training time
#'   type and supplied in strictly increasing order.
#' @param interval_level Probability level for posterior intervals.
#' @param include_covariance If `TRUE`, include full posterior covariance
#'   matrices.
#' @param variance_tolerance Relative tolerance for posterior variances.
#' @param observation_noise_variance Optional non-negative scalar or vector of
#'   observation-noise variances at `future_time`, passed to [predict_gp()].
#'   Required when the model was fitted with observation-specific training
#'   noise.
#'
#' @return An object of class `gaussianprocesses_time_forecast`.
#' @examples
#' time <- as.Date("2026-01-01") + 0:23
#' model <- fit_time_series_gp(
#'   time,
#'   sin(2 * pi * (0:23) / 12),
#'   kernel = time_series_kernel(period = 12),
#'   noise_variance = 0.05
#' )
#'
#' forecast <- forecast_gp(model, as.Date("2026-01-01") + c(10.5, 25, 30))
#' data.frame(
#'   region = forecast$region,
#'   mean = forecast$mean,
#'   forecast$prediction_interval
#' )
#' @template stable
#' @export
forecast_gp <- function(
    model,
    future_time,
    interval_level = 0.95,
    include_covariance = FALSE,
    variance_tolerance = sqrt(.Machine$double.eps),
    observation_noise_variance = NULL) {
  model <- .prepare_time_model(model)
  numeric_time <- .convert_future_time(model, future_time)

  predict_at <- if (inherits(model, "gaussianprocesses_state_space_model")) {
    .predict_state_space
  } else {
    predict_gp
  }
  prediction <- predict_at(
    model,
    numeric_time,
    interval_level = interval_level,
    include_covariance = include_covariance,
    variance_tolerance = variance_tolerance,
    observation_noise_variance = observation_noise_variance
  )

  training_range <- range(model$time$training_index)
  region <- ifelse(
    numeric_time < training_range[[1L]],
    "backward_extrapolation",
    ifelse(
      numeric_time > training_range[[2L]],
      "forecast",
      "interpolation"
    )
  )

  result <- unclass(prediction)
  result$time <- future_time
  result$time_index <- numeric_time
  result$region <- region

  structure(
    result,
    class = "gaussianprocesses_time_forecast"
  )
}


#' @export
print.gaussianprocesses_time_forecast <- function(x, ...) {
  cat("Gaussian-process time-series forecast\n")
  cat("  prediction points: ", length(x$mean), "\n", sep = "")
  cat(
    "  interpolation points: ",
    sum(x$region == "interpolation"),
    "\n",
    sep = ""
  )
  cat(
    "  forward forecast points: ",
    sum(x$region == "forecast"),
    "\n",
    sep = ""
  )
  cat(
    "  backward extrapolation points: ",
    sum(x$region == "backward_extrapolation"),
    "\n",
    sep = ""
  )
  invisible(x)
}


#' Rolling-origin evaluation for time-aware Gaussian processes
#'
#' Performs expanding-window rolling-origin evaluation with fixed kernel and
#' noise hyperparameters. Each split is fitted using only observations
#' available at that origin, preventing look-ahead leakage.
#'
#' @param time Numeric, Date, or POSIXt vector in strictly increasing order.
#' @param y Numeric response vector.
#' @param kernel Fixed Gaussian-process kernel specification.
#' @param initial_window Number of observations in the first training window.
#' @param horizon Number of future observations predicted at each origin.
#' @param step Number of observations between consecutive origins.
#' @param noise_variance Non-negative observation-noise variance: a scalar,
#'   or one value per observation. With one value per observation, each
#'   training window uses the values of its rows and each forecast uses the
#'   values of its target rows.
#' @param mean Gaussian-process mean specification.
#' @param unit Time unit passed to [gp_time_index()].
#' @param interval_level Probability level for prediction intervals.
#' @param ... Additional arguments passed to [fit_time_series_gp()], such
#'   as `method = "state_space"`.
#'
#' @details
#' Fixed hyperparameters prevent future observations from influencing earlier
#' origins through repeated whole-series hyperparameter tuning.
#' Hyperparameters can be selected separately using an appropriate training
#' period.
#'
#' @return A data frame with class `gaussianprocesses_backtest`.
#' @examples
#' time <- 0:23
#' y <- sin(2 * pi * time / 12)
#'
#' backtest <- rolling_origin_gp(
#'   time,
#'   y,
#'   kernel = time_series_kernel(period = 12),
#'   initial_window = 12,
#'   horizon = 3,
#'   step = 3,
#'   noise_variance = 0.05
#' )
#' head(backtest[, c("origin_index", "horizon", "observed", "mean", "covered")])
#' @template stable
#' @export
rolling_origin_gp <- function(
    time,
    y,
    kernel,
    initial_window,
    horizon = 1L,
    step = 1L,
    noise_variance = 1e-6,
    mean = zero_mean(),
    unit = "auto",
    interval_level = 0.95,
    ...) {
  index <- gp_time_index(time, unit = unit)
  .assert_numeric_vector(y, name = "y")

  if (length(index) != length(y)) {
    .gp_abort(
      "'time' and 'y' must contain the same number of observations.",
      "gaussianprocesses_validation_error"
    )
  }

  initial_window <- .assert_positive_integer(
    initial_window,
    "initial_window"
  )
  horizon <- .assert_positive_integer(horizon, "horizon")
  step <- .assert_positive_integer(step, "step")
  .validate_interval_level(interval_level)
  .assert_kernel(kernel)
  .assert_mean_function(mean)
  noise <- .resolve_training_noise_variance(
    noise_variance,
    n_observations = length(y)
  )
  per_observation <- length(noise_variance) > 1L

  n <- length(y)

  if (initial_window >= n) {
    .gp_abort(
      "'initial_window' must be smaller than the number of observations.",
      "gaussianprocesses_validation_error"
    )
  }

  last_origin <- n - horizon

  if (initial_window > last_origin) {
    .gp_abort(
      paste0(
        "'initial_window' and 'horizon' leave no complete ",
        "rolling-origin evaluation split."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  origins <- seq(
    from = initial_window,
    to = last_origin,
    by = step
  )

  blocks <- vector("list", length(origins))

  for (block_index in seq_along(origins)) {
    origin <- origins[[block_index]]
    training_rows <- seq_len(origin)
    target_rows <- origin + seq_len(horizon)

    model <- fit_time_series_gp(
      time[training_rows],
      y[training_rows],
      kernel = kernel,
      unit = unit,
      noise_variance = if (per_observation) {
        noise[training_rows]
      } else {
        noise_variance
      },
      mean = mean,
      ...
    )

    forecast <- forecast_gp(
      model,
      time[target_rows],
      interval_level = interval_level,
      observation_noise_variance = if (per_observation) {
        noise[target_rows]
      }
    )

    blocks[[block_index]] <- data.frame(
      origin_index = rep(origin, horizon),
      target_index = target_rows,
      horizon = seq_len(horizon),
      origin_time = rep(time[origin], horizon),
      target_time = time[target_rows],
      observed = y[target_rows],
      mean = forecast$mean,
      observation_variance = forecast$observation_variance,
      observation_sd = forecast$observation_sd,
      lower = forecast$prediction_interval[, "lower"],
      upper = forecast$prediction_interval[, "upper"],
      covered =
        y[target_rows] >= forecast$prediction_interval[, "lower"] &
        y[target_rows] <= forecast$prediction_interval[, "upper"],
      log_predictive_density = stats::dnorm(
        y[target_rows],
        mean = forecast$mean,
        sd = forecast$observation_sd,
        log = TRUE
      ),
      check.names = FALSE
    )
  }

  result <- do.call(rbind, blocks)
  rownames(result) <- NULL

  structure(
    result,
    interval_level = interval_level,
    unit = attr(index, "unit"),
    class = c(
      "gaussianprocesses_backtest",
      "data.frame"
    )
  )
}


#' Uncertainty-aware forecast metrics
#'
#' Summarizes a rolling-origin GP backtest with both point and probabilistic
#' forecast metrics.
#'
#' @param backtest A `gaussianprocesses_backtest` returned by
#'   [rolling_origin_gp()].
#'
#' @return A named list containing RMSE, MAE, standardized RMSE, empirical
#'   interval coverage, mean interval width, and mean log predictive density.
#' @examples
#' time <- 0:23
#' backtest <- rolling_origin_gp(
#'   time,
#'   sin(2 * pi * time / 12),
#'   kernel = time_series_kernel(period = 12),
#'   initial_window = 12,
#'   horizon = 3,
#'   step = 3,
#'   noise_variance = 0.05
#' )
#'
#' forecast_metrics(backtest)
#' @template stable
#' @export
forecast_metrics <- function(backtest) {
  if (!inherits(backtest, "gaussianprocesses_backtest")) {
    .gp_abort(
      "'backtest' must be returned by rolling_origin_gp().",
      "gaussianprocesses_validation_error"
    )
  }

  residual <- backtest$observed - backtest$mean
  variance <- backtest$observation_variance

  if (any(!is.finite(variance)) || any(variance <= 0)) {
    .gp_abort(
      paste0(
        "Forecast observation variances must be finite and strictly positive ",
        "for probabilistic metrics."
      ),
      "gaussianprocesses_numerical_error"
    )
  }

  list(
    n_forecasts = nrow(backtest),
    rmse = sqrt(mean(residual^2)),
    mae = mean(abs(residual)),
    standardized_rmse = sqrt(mean(residual^2 / variance)),
    empirical_coverage = mean(backtest$covered),
    interval_level = attr(backtest, "interval_level"),
    mean_interval_width = mean(backtest$upper - backtest$lower),
    mean_log_predictive_density = mean(
      backtest$log_predictive_density
    )
  )
}


#' @export
print.gaussianprocesses_backtest <- function(x, ...) {
  cat("Gaussian-process rolling-origin backtest\n")
  cat("  forecasts: ", nrow(x), "\n", sep = "")
  cat(
    "  origins: ",
    length(unique(x$origin_index)),
    "\n",
    sep = ""
  )
  cat(
    "  maximum horizon: ",
    max(x$horizon),
    "\n",
    sep = ""
  )
  invisible(x)
}
