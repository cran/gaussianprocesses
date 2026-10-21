# Heteroscedastic Gaussian-process approximation -----------------------------

.clip_noise_variance <- function(
    variance,
    noise_floor,
    noise_ceiling) {
  clipped <- pmax(variance, noise_floor)

  if (is.finite(noise_ceiling)) {
    clipped <- pmin(clipped, noise_ceiling)
  }

  clipped
}


.log_squared_residual_target <- function(
    residual,
    noise_floor) {
  # For e ~ N(0, sigma^2):
  # E[log(e^2)] = log(sigma^2) + digamma(1/2) + log(2).
  correction <- digamma(0.5) + log(2)

  log(pmax(residual^2, noise_floor)) - correction
}


# How the noise GP's targets are formed from the leave-one-out residuals
# e_i = y_i - mu_{-i} of the current mean model. Under the model
# e_i ~ N(0, v_i + sigma_i^2), where v_i is the leave-one-out variance of the
# latent function, so log(e_i^2) estimates log(v_i + sigma_i^2) rather than
# log(sigma_i^2):
#   "residual"      log(e_i^2) - E[log chi^2_1]; biased upwards by
#                   log(1 + v_i / sigma_i^2). Used before this was studied.
#   "standardized"  the same minus log((v_i + s_i^2) / s_i^2), with s_i^2 the
#                   current noise; its fixed point is unbiased.
#   "subtracted"    log(max(e_i^2 - v_i, floor)) - E[log chi^2_1].
# inst/benchmarks/heteroscedastic-noise.R compares them; the option exists
# for that study and is not part of the documented API.
.heteroscedastic_noise_targets <- c("standardized", "residual", "subtracted")


.heteroscedastic_noise_target <- function() {
  rule <- getOption(
    "gaussianprocesses.heteroscedastic_target",
    .heteroscedastic_noise_targets[[1L]]
  )

  if (!is.character(rule) || length(rule) != 1L ||
      !(rule %in% .heteroscedastic_noise_targets)) {
    .gp_abort(
      sprintf(
        "Option 'gaussianprocesses.heteroscedastic_target' must be one of %s.",
        paste0("\"", .heteroscedastic_noise_targets, "\"", collapse = ", ")
      ),
      "gaussianprocesses_validation_error"
    )
  }

  rule
}


# `loo` comes from loo_gp() on a mean model fitted with `training_noise`, so
# loo$variance = v_i + s_i^2.
.noise_targets <- function(loo, training_noise, noise_floor, rule) {
  latent_variance <- pmax(loo$variance - training_noise, 0)

  switch(
    rule,
    residual = .log_squared_residual_target(loo$residual, noise_floor),
    standardized = .log_squared_residual_target(loo$residual, noise_floor) -
      log(loo$variance / training_noise),
    subtracted = log(
      pmax(loo$residual^2 - latent_variance, noise_floor)
    ) - (digamma(0.5) + log(2))
  )
}


#' Fit an approximate heteroscedastic Gaussian process
#'
#' Fits two coupled Gaussian processes: one for the latent regression function
#' and one for the logarithm of input-dependent observation-noise variance.
#' Noise targets are constructed from exact leave-one-out residuals of the
#' current mean GP, standardized by their leave-one-out variance so that the
#' latent function's uncertainty is not counted as noise, and corrected for
#' the known mean of a log chi-squared variable with one degree of freedom.
#'
#' This is a deterministic plug-in approximation inspired by most-likely
#' heteroscedastic GP methods. It is not full joint Bayesian inference over the
#' latent function and noise process.
#'
#' @param x Numeric vector or matrix of training inputs.
#' @param y Numeric response vector.
#' @param kernel Kernel specification for the latent regression function.
#' @param noise_kernel Kernel specification for log observation-noise variance.
#' @param mean Mean specification for the latent regression function, with
#'   fixed coefficients.
#' @param initial_noise_variance Optional positive scalar used for the initial
#'   homoscedastic mean fit. If `NULL`, a small fraction of response variance
#'   is used, bounded below by `noise_floor`.
#' @param noise_floor Strictly positive lower bound for estimated noise
#'   variances.
#' @param noise_ceiling Optional upper bound for estimated noise variances.
#' @param damping Value in \eqn{(0, 1]} controlling how strongly each noise
#'   update replaces the previous estimate.
#' @param max_iterations Maximum number of mean/noise refinement iterations.
#'   In a simulation study every fit that converged did so within 50
#'   iterations, and most within 15.
#' @param convergence_tolerance Convergence tolerance on the maximum absolute
#'   log change in training noise variance.
#' @param initial_jitter Non-negative jitter tried first, relative to the
#'   covariance scale (the mean of its diagonal).
#' @param fallback_jitter Positive relative jitter tried after the first
#'   failed factorization.
#' @param jitter_multiplier Multiplicative jitter escalation factor.
#' @param max_attempts Maximum Cholesky attempts.
#' @param symmetry_tolerance Relative tolerance for checking that
#'   covariance matrices are symmetric.
#'
#' @details
#' For a Gaussian residual \eqn{e \sim N(0, \sigma^2)}, \eqn{\log(e^2)} has
#' additive log-chi-squared noise with known mean \eqn{\psi(1/2) + \log 2} and
#' variance \eqn{\psi_1(1/2)}. The implementation uses these known moments
#' when constructing and fitting the log-noise GP.
#'
#' A leave-one-out residual \eqn{e_i} has variance \eqn{v_i + \sigma_i^2},
#' where \eqn{v_i} is the leave-one-out variance of the latent function, so
#' \eqn{\log(e_i^2)} alone overestimates the noise. The target is therefore
#' \eqn{\log(e_i^2) - \psi(1/2) - \log 2 - \log((v_i + s_i^2) / s_i^2)}, with
#' \eqn{s_i^2} the current noise estimate, which is unbiased at a fixed point
#' of the iteration. The research note
#' `system.file("notes", "heteroscedastic-gp.md", package =
#' "gaussianprocesses")` derives the correction and reports a simulation
#' study of its bias and convergence.
#'
#' The approach alternates plug-in estimates rather than integrating jointly
#' over the latent function and noise process. If the estimate has not
#' converged after `max_iterations` iterations, the function warns with class
#' `gaussianprocesses_convergence_warning` and returns the last iterate with
#' `converged = FALSE`.
#'
#' @return An object of class
#'   `gaussianprocesses_heteroscedastic_model`.
#' @seealso [gaussianprocesses_model] for the structure, versioning, and
#'   persistence of fitted models.
#' @examples
#' simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 60)
#'
#' model <- fit_heteroscedastic_gp(
#'   simulation$x,
#'   simulation$observed,
#'   kernel = matern52_kernel(length_scale = 0.2),
#'   noise_kernel = rbf_kernel(length_scale = 0.4),
#'   max_iterations = 30
#' )
#' model
#'
#' # The estimated noise variance rises with x, as in the simulation.
#' rows <- c(1, 20, 40, 60)
#' cbind(
#'   x = simulation$x[rows],
#'   true = simulation$noise_variance[rows],
#'   estimated = model$training_noise_variance[rows]
#' )
#' @template experimental
#' @export
fit_heteroscedastic_gp <- function(
    x,
    y,
    kernel,
    noise_kernel,
    mean = zero_mean(),
    initial_noise_variance = NULL,
    noise_floor = 1e-6,
    noise_ceiling = Inf,
    damping = 0.5,
    max_iterations = 50L,
    convergence_tolerance = 1e-3,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps)) {
  data <- .validate_gp_training_data(x, y)
  .assert_kernel(kernel)
  .assert_kernel(noise_kernel, name = "noise_kernel")
  .assert_kernel_dimensions(kernel, ncol(data$x))
  .assert_kernel_dimensions(noise_kernel, ncol(data$x))
  .assert_mean_function(mean)
  .assert_fixed_mean(mean, "fit_heteroscedastic_gp()")

  if (nrow(data$x) < 3L) {
    .gp_abort(
      "Heteroscedastic fitting requires at least three observations.",
      "gaussianprocesses_validation_error"
    )
  }

  .assert_numeric_scalar(
    noise_floor,
    "noise_floor",
    lower = 0,
    lower_open = TRUE
  )

  .assert_numeric_scalar(
    noise_ceiling,
    "noise_ceiling",
    lower = noise_floor,
    lower_open = TRUE,
    finite = FALSE
  )

  .assert_numeric_scalar(
    damping,
    "damping",
    lower = 0,
    lower_open = TRUE
  )

  if (damping > 1) {
    .gp_abort(
      "'damping' must be less than or equal to 1.",
      "gaussianprocesses_validation_error"
    )
  }

  max_iterations <- .assert_positive_integer(
    max_iterations,
    "max_iterations"
  )
  .assert_numeric_scalar(
    convergence_tolerance,
    "convergence_tolerance",
    lower = 0,
    lower_open = TRUE
  )

  if (is.null(initial_noise_variance)) {
    response_variance <- stats::var(data$y)

    initial_noise_variance <- if (
      is.finite(response_variance) &&
        response_variance > 0
    ) {
      max(0.05 * response_variance, noise_floor)
    } else {
      noise_floor
    }
  } else {
    .assert_numeric_scalar(
      initial_noise_variance,
      "initial_noise_variance",
      lower = 0,
      lower_open = TRUE
    )
  }

  initial_noise_variance <- min(
    max(initial_noise_variance, noise_floor),
    noise_ceiling
  )

  training_noise <- rep(
    initial_noise_variance,
    nrow(data$x)
  )

  mean_model <- fit_gp(
    data$x,
    data$y,
    kernel = kernel,
    noise_variance = training_noise,
    mean = mean,
    initial_jitter = initial_jitter,
    fallback_jitter = fallback_jitter,
    jitter_multiplier = jitter_multiplier,
    max_attempts = max_attempts,
    symmetry_tolerance = symmetry_tolerance
  )

  history <- vector("list", max_iterations)
  noise_model <- NULL
  converged <- FALSE
  iterations_completed <- 0L

  # If squared Gaussian residuals are used as noisy observations of sigma^2,
  # log residual-squared has additive log-chi-squared noise. Its variance is
  # trigamma(1/2), which we use as the observation variance for the noise GP.
  log_residual_observation_variance <- trigamma(0.5)

  target_rule <- .heteroscedastic_noise_target()

  for (iteration in seq_len(max_iterations)) {
    loo <- loo_gp(mean_model)
    target <- .noise_targets(
      loo,
      training_noise = training_noise,
      noise_floor = noise_floor,
      rule = target_rule
    )

    noise_mean <- constant_mean(mean(target))

    noise_model <- fit_gp(
      data$x,
      target,
      kernel = noise_kernel,
      noise_variance = log_residual_observation_variance,
      mean = noise_mean,
      initial_jitter = initial_jitter,
      fallback_jitter = fallback_jitter,
      jitter_multiplier = jitter_multiplier,
      max_attempts = max_attempts,
      symmetry_tolerance = symmetry_tolerance
    )

    noise_prediction <- predict_gp(
      noise_model,
      data$x,
      observation_noise_variance = 0
    )

    proposed_noise <- .clip_noise_variance(
      exp(noise_prediction$mean),
      noise_floor = noise_floor,
      noise_ceiling = noise_ceiling
    )

    updated_noise <- (1 - damping) * training_noise +
      damping * proposed_noise

    updated_noise <- .clip_noise_variance(
      updated_noise,
      noise_floor = noise_floor,
      noise_ceiling = noise_ceiling
    )

    log_change <- max(
      abs(log(updated_noise) - log(training_noise))
    )

    history[[iteration]] <- list(
      iteration = iteration,
      max_log_noise_change = log_change,
      min_noise_variance = min(updated_noise),
      median_noise_variance = stats::median(updated_noise),
      max_noise_variance = max(updated_noise),
      mean_loo_standardized_residual =
        mean(loo$standardized_residual),
      sd_loo_standardized_residual =
        stats::sd(loo$standardized_residual)
    )

    training_noise <- updated_noise

    mean_model <- fit_gp(
      data$x,
      data$y,
      kernel = kernel,
      noise_variance = training_noise,
      mean = mean,
      initial_jitter = initial_jitter,
      fallback_jitter = fallback_jitter,
      jitter_multiplier = jitter_multiplier,
      max_attempts = max_attempts,
      symmetry_tolerance = symmetry_tolerance
    )

    iterations_completed <- iteration

    if (log_change <= convergence_tolerance) {
      converged <- TRUE
      break
    }
  }

  history <- history[seq_len(iterations_completed)]

  if (!converged) {
    warning(
      warningCondition(
        sprintf(
          paste0(
            "The heteroscedastic noise estimate did not converge in %d ",
            "iterations: the last largest change in log noise variance ",
            "was %.3g, above convergence_tolerance = %.3g. Increase ",
            "max_iterations."
          ),
          iterations_completed,
          history[[iterations_completed]]$max_log_noise_change,
          convergence_tolerance
        ),
        class = "gaussianprocesses_convergence_warning"
      )
    )
  }

  structure(
    list(
      x = data$x,
      y = data$y,
      mean_model = mean_model,
      noise_model = noise_model,
      kernel = kernel,
      noise_kernel = noise_kernel,
      training_noise_variance = training_noise,
      noise_floor = noise_floor,
      noise_ceiling = noise_ceiling,
      damping = damping,
      tolerance = convergence_tolerance,
      converged = converged,
      iterations = iterations_completed,
      history = history,
      log_residual_observation_variance =
        log_residual_observation_variance,
      approximation = "iterative_log_residual_plugin",
      noise_target = target_rule,
      schema_version = .model_schema_version,
      package_version = .package_version()
    ),
    class = "gaussianprocesses_heteroscedastic_model"
  )
}


#' Predict from an approximate heteroscedastic GP
#'
#' The latent-function posterior and the input-dependent noise process are
#' predicted separately. Observation variance is the sum of latent posterior
#' variance and the plug-in estimated noise variance. Uncertainty about the
#' log-noise process is returned separately and is not folded into the
#' observation variance.
#'
#' @param model A fitted heteroscedastic GP model.
#' @param x Numeric vector or matrix of prediction inputs.
#' @param interval_level Probability level for returned intervals.
#' @param include_covariance If `TRUE`, include the full latent and
#'   observation covariance matrices.
#' @param variance_tolerance Relative tolerance for latent posterior variances.
#'
#' @return An object of class
#'   `gaussianprocesses_heteroscedastic_prediction`. It has the fields that
#'   every prediction in the package shares (see [predict_gp()]), with the
#'   plug-in noise estimate as `observation_noise_variance`, and in addition
#'   `noise_log_variance`, its posterior variance
#'   `noise_log_variance_uncertainty`, and `noise_variance_interval` from the
#'   log-noise GP.
#' @examples
#' simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 60)
#' model <- fit_heteroscedastic_gp(
#'   simulation$x,
#'   simulation$observed,
#'   kernel = matern52_kernel(length_scale = 0.2),
#'   noise_kernel = rbf_kernel(length_scale = 0.4),
#'   max_iterations = 30
#' )
#'
#' prediction <- predict_heteroscedastic_gp(model, c(0.1, 0.9))
#'
#' # Observation variance = latent variance + estimated noise variance.
#' cbind(
#'   latent = prediction$latent_variance,
#'   noise = prediction$observation_noise_variance,
#'   observation = prediction$observation_variance
#' )
#' @template experimental
#' @export
predict_heteroscedastic_gp <- function(
    model,
    x,
    interval_level = 0.95,
    include_covariance = FALSE,
    variance_tolerance = sqrt(.Machine$double.eps)) {
  model <- .prepare_heteroscedastic_model(model)
  .validate_interval_level(interval_level)

  mean_prediction <- predict_gp(
    model$mean_model,
    x,
    interval_level = interval_level,
    include_covariance = include_covariance,
    variance_tolerance = variance_tolerance,
    observation_noise_variance = 0
  )

  noise_prediction <- predict_gp(
    model$noise_model,
    x,
    interval_level = interval_level,
    include_covariance = FALSE,
    variance_tolerance = variance_tolerance,
    observation_noise_variance = 0
  )

  estimated_noise_variance <- .clip_noise_variance(
    exp(noise_prediction$mean),
    noise_floor = model$noise_floor,
    noise_ceiling = model$noise_ceiling
  )

  noise_variance_interval <- exp(
    noise_prediction$latent_interval
  )
  noise_variance_interval[, "lower"] <- .clip_noise_variance(
    noise_variance_interval[, "lower"],
    model$noise_floor,
    model$noise_ceiling
  )
  noise_variance_interval[, "upper"] <- .clip_noise_variance(
    noise_variance_interval[, "upper"],
    model$noise_floor,
    model$noise_ceiling
  )

  observation_variance <-
    mean_prediction$latent_variance +
    estimated_noise_variance
  observation_sd <- sqrt(observation_variance)
  quantile <- stats::qnorm((1 + interval_level) / 2)

  prediction_interval <- cbind(
    lower = mean_prediction$mean -
      quantile * observation_sd,
    upper = mean_prediction$mean +
      quantile * observation_sd
  )

  result <- list(
    x = mean_prediction$x,
    mean = mean_prediction$mean,
    latent_variance = mean_prediction$latent_variance,
    latent_sd = mean_prediction$latent_sd,
    latent_interval = mean_prediction$latent_interval,
    observation_noise_variance = estimated_noise_variance,
    noise_log_variance = noise_prediction$mean,
    noise_log_variance_uncertainty =
      noise_prediction$latent_variance,
    noise_variance_interval = noise_variance_interval,
    observation_variance = observation_variance,
    observation_sd = observation_sd,
    prediction_interval = prediction_interval,
    interval_level = interval_level
  )

  if (include_covariance) {
    result$latent_covariance <-
      mean_prediction$latent_covariance
    result$observation_covariance <- mean_prediction$latent_covariance +
      diag(estimated_noise_variance, length(estimated_noise_variance))
  }

  structure(
    result,
    class = "gaussianprocesses_heteroscedastic_prediction"
  )
}


#' @export
print.gaussianprocesses_heteroscedastic_model <- function(x, ...) {
  x <- .prepare_heteroscedastic_model(x, "x")
  cat("Approximate heteroscedastic Gaussian-process model\n")
  cat("  observations: ", nrow(x$x), "\n", sep = "")
  cat(
    "  iterations: ",
    x$iterations,
    "\n",
    sep = ""
  )
  cat(
    "  converged: ",
    if (x$converged) "yes" else "no",
    "\n",
    sep = ""
  )
  cat(
    "  training noise variance range: [",
    format(
      signif(min(x$training_noise_variance), 6),
      trim = TRUE
    ),
    ", ",
    format(
      signif(max(x$training_noise_variance), 6),
      trim = TRUE
    ),
    "]\n",
    sep = ""
  )
  cat(
    "  approximation: ",
    x$approximation,
    "\n",
    sep = ""
  )
  invisible(x)
}


#' @export
print.gaussianprocesses_heteroscedastic_prediction <- function(
    x,
    ...) {
  cat("Heteroscedastic Gaussian-process prediction\n")
  cat("  prediction points: ", length(x$mean), "\n", sep = "")
  cat(
    "  estimated noise variance range: [",
    format(signif(min(x$observation_noise_variance), 6), trim = TRUE),
    ", ",
    format(signif(max(x$observation_noise_variance), 6), trim = TRUE),
    "]\n",
    sep = ""
  )
  invisible(x)
}
