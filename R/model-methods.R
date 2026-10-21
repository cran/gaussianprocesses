# Methods for base R generics -------------------------------------------------

#' Base R methods for fitted Gaussian-process models
#'
#' Methods for [stats::predict()], [stats::logLik()], [stats::nobs()],
#' [stats::fitted()], [stats::residuals()], and [summary()] on exact
#' (`gaussianprocesses_model`), sparse (`gaussianprocesses_sparse_model`),
#' heteroscedastic (`gaussianprocesses_heteroscedastic_model`), latent
#' (`gaussianprocesses_latent_model`), and state-space
#' (`gaussianprocesses_state_space_model`) models.
#'
#' @param object A fitted model.
#' @param newdata Prediction inputs. If missing, the training inputs (or,
#'   for time-series models, the training times) are used.
#' @param x A summary object.
#' @param ... Further arguments passed to the package's prediction function:
#'   [predict_gp()], [forecast_gp()], [predict_sparse_gp()],
#'   [predict_heteroscedastic_gp()], or [predict_latent_gp()].
#'
#' @details
#' `predict()` returns exactly what the package's prediction function for the
#' model returns. For models from [fit_time_series_gp()], `newdata` are time
#' values and `predict()` calls [forecast_gp()]; this includes every
#' state-space model with time metadata.
#'
#' `logLik()` returns the log marginal likelihood
#' \eqn{\log p(y \mid X, \hat\theta)}: [log_marginal_likelihood()] for exact
#' and state-space models and, for sparse models, the FITC approximate
#' marginal likelihood or the VFE lower bound. This is a
#' type-II likelihood, with the latent function integrated out and the
#' hyperparameters \eqn{\hat\theta} plugged in. Its `df` attribute counts the
#' hyperparameters estimated from the data: the parameters optimized by
#' [optimize_gp()] or [optimize_sparse_gp()], and 0 for models whose
#' hyperparameters were supplied,
#' plus the mean coefficients when they are estimated. For restricted
#' likelihoods (REML estimates and vague coefficient priors) the `nobs`
#' attribute is \eqn{n - p}, as for restricted likelihoods elsewhere in R,
#' and such likelihoods are comparable only between models with the same
#' mean.
#' [stats::AIC()] and [stats::BIC()] therefore work, but comparing kernels
#' by them is a heuristic: it is not an ordinary maximum-likelihood
#' comparison, and the penalty does not account for the flexibility of the
#' kernel. Heteroscedastic models have no `logLik()` method, because their
#' noise function is a plug-in estimate rather than a likelihood maximizer.
#' For latent models it is the Laplace approximation to the marginal
#' likelihood, and `df` counts the parameters optimized by
#' [optimize_latent_gp()] or the mean coefficients estimated by
#' [fit_latent_gp()]. For state-space models `df` counts the parameters
#' optimized by [optimize_time_series_gp()].
#'
#' `fitted()` returns the posterior mean of the latent function at the
#' training inputs, and `residuals()` returns the observed responses minus
#' those fitted values. For exact models the fitted values are
#' \eqn{y - (\Sigma_\varepsilon + \delta I)\alpha}, which equals
#' \eqn{m(X) + K\alpha} without forming \eqn{K}; \eqn{\delta} is the
#' numerical jitter. For state-space models they are the Kalman-smoother
#' means. For latent models `fitted()` returns the posterior mode
#' \eqn{\hat f} of the latent function, and there is no `residuals()`
#' method, because the responses are not on the scale of \eqn{f}.
#'
#' `summary()` collects the training size and dimension, kernel, mean, noise,
#' log marginal likelihood, optimization status, and numerical status, and
#' prints them compactly.
#'
#' @return `predict()`: a prediction object. `logLik()`: an object of class
#'   `logLik`. `nobs()`: the number of training observations. `fitted()`
#'   and `residuals()`: numeric vectors with one value per training
#'   observation. `summary()`: an object of class
#'   `summary.gaussianprocesses_model`.
#'
#' @examples
#' x <- seq(-2, 2, length.out = 20)
#' y <- sin(2 * x) + 0.1 * cos(7 * x)
#' model <- optimize_gp(
#'   x,
#'   y,
#'   kernel = rbf_kernel(),
#'   noise_variance = 0.05,
#'   n_starts = 1
#' )
#'
#' summary(model)
#' predict(model, c(-1, 0.5))$mean
#' logLik(model)
#' AIC(model)
#' head(residuals(model))
#' @template stable
#' @name gp_model_methods
#' @importFrom stats predict logLik nobs fitted residuals
NULL


# predict ---------------------------------------------------------------------

#' @rdname gp_model_methods
#' @export
predict.gaussianprocesses_model <- function(object, newdata, ...) {
  object <- .prepare_gp_model(object, "object")

  if (!is.null(object$time)) {
    if (missing(newdata)) {
      newdata <- object$time$training_time
    }

    return(forecast_gp(object, newdata, ...))
  }

  if (missing(newdata)) {
    newdata <- object$x
  }

  predict_gp(object, newdata, ...)
}


#' @rdname gp_model_methods
#' @export
predict.gaussianprocesses_sparse_model <- function(object, newdata, ...) {
  object <- .prepare_sparse_gp_model(object, "object")

  if (missing(newdata)) {
    newdata <- object$x
  }

  predict_sparse_gp(object, newdata, ...)
}


#' @rdname gp_model_methods
#' @export
predict.gaussianprocesses_state_space_model <- function(object, newdata, ...) {
  object <- .prepare_state_space_model(object, "object")

  if (!is.null(object$time)) {
    if (missing(newdata)) {
      newdata <- object$time$training_time
    }

    return(forecast_gp(object, newdata, ...))
  }

  if (missing(newdata)) {
    newdata <- object$x
  }

  .predict_state_space(object, newdata, ...)
}


#' @rdname gp_model_methods
#' @export
predict.gaussianprocesses_heteroscedastic_model <- function(
    object,
    newdata,
    ...) {
  object <- .prepare_heteroscedastic_model(object, "object")

  if (missing(newdata)) {
    newdata <- object$x
  }

  predict_heteroscedastic_gp(object, newdata, ...)
}


#' @rdname gp_model_methods
#' @export
predict.gaussianprocesses_latent_model <- function(object, newdata, ...) {
  object <- .prepare_latent_model(object, "object")

  if (missing(newdata)) {
    newdata <- object$x
  }

  predict_latent_gp(object, newdata, ...)
}


# logLik and nobs -------------------------------------------------------------

# Hyperparameters estimated from the data: those optimize_gp() optimized.
.estimated_hyperparameter_count <- function(model) {
  if (is.null(model$optimization)) {
    return(0L)
  }

  length(model$optimization$optimized_parameters)
}


.new_log_lik <- function(value, df, n_observations) {
  structure(
    value,
    df = df,
    nobs = n_observations,
    class = "logLik"
  )
}


#' @rdname gp_model_methods
#' @export
logLik.gaussianprocesses_model <- function(object, ...) {
  object <- .prepare_gp_model(object, "object")

  n_coefficients <- length(object$mean_fit$coefficients)
  estimated <- identical(object$mean_fit$treatment, "estimated")
  restricted <- identical(.likelihood_type(object), "restricted")

  .new_log_lik(
    log_marginal_likelihood(object),
    df = .estimated_hyperparameter_count(object) +
      if (estimated) n_coefficients else 0L,
    n_observations = object$n_observations -
      if (restricted) n_coefficients else 0L
  )
}


#' @rdname gp_model_methods
#' @export
logLik.gaussianprocesses_sparse_model <- function(object, ...) {
  object <- .prepare_sparse_gp_model(object, "object")

  .new_log_lik(
    object$log_marginal_likelihood,
    df = .estimated_hyperparameter_count(object),
    n_observations = object$n_observations
  )
}


#' @rdname gp_model_methods
#' @export
logLik.gaussianprocesses_state_space_model <- function(object, ...) {
  object <- .prepare_state_space_model(object, "object")

  .new_log_lik(
    object$log_marginal_likelihood,
    df = .estimated_hyperparameter_count(object),
    n_observations = object$n_observations
  )
}


#' @rdname gp_model_methods
#' @export
logLik.gaussianprocesses_heteroscedastic_model <- function(object, ...) {
  .gp_abort(
    paste0(
      "logLik() is not defined for heteroscedastic models: the noise ",
      "function is a plug-in estimate, not a maximizer of a likelihood. ",
      "Use the leave-one-out or held-out predictive density instead."
    ),
    "gaussianprocesses_validation_error"
  )
}


#' @rdname gp_model_methods
#' @export
logLik.gaussianprocesses_latent_model <- function(object, ...) {
  object <- .prepare_latent_model(object, "object")

  # Coefficients estimated by fit_latent_gp(); optimize_latent_gp() counts
  # them among its optimized parameters.
  estimated_coefficients <- if (is.null(object$optimization)) {
    length(object$mean_fit$coefficients)
  } else {
    0L
  }

  .new_log_lik(
    object$log_marginal_likelihood,
    df = .estimated_hyperparameter_count(object) + estimated_coefficients,
    n_observations = object$n_observations
  )
}


#' @rdname gp_model_methods
#' @export
nobs.gaussianprocesses_model <- function(object, ...) {
  .prepare_gp_model(object, "object")$n_observations
}


#' @rdname gp_model_methods
#' @export
nobs.gaussianprocesses_sparse_model <- function(object, ...) {
  .prepare_sparse_gp_model(object, "object")$n_observations
}


#' @rdname gp_model_methods
#' @export
nobs.gaussianprocesses_state_space_model <- function(object, ...) {
  .prepare_state_space_model(object, "object")$n_observations
}


#' @rdname gp_model_methods
#' @export
nobs.gaussianprocesses_heteroscedastic_model <- function(object, ...) {
  nrow(.prepare_heteroscedastic_model(object, "object")$x)
}


#' @rdname gp_model_methods
#' @export
nobs.gaussianprocesses_latent_model <- function(object, ...) {
  .prepare_latent_model(object, "object")$n_observations
}


# fitted and residuals --------------------------------------------------------

# Posterior latent means at the training inputs. Since
# (K + Sigma + delta I) alpha = y - m(X), the mean m(X) + K alpha equals
# y - (Sigma + delta I) alpha.
.exact_fitted_values <- function(model) {
  model$y - (model$training_noise_variance + model$jitter) * model$alpha
}


#' @rdname gp_model_methods
#' @export
fitted.gaussianprocesses_model <- function(object, ...) {
  .exact_fitted_values(.prepare_gp_model(object, "object"))
}


#' @rdname gp_model_methods
#' @export
fitted.gaussianprocesses_sparse_model <- function(object, ...) {
  object <- .prepare_sparse_gp_model(object, "object")
  predict_sparse_gp(object, object$x)$mean
}


#' @rdname gp_model_methods
#' @export
fitted.gaussianprocesses_state_space_model <- function(object, ...) {
  .prepare_state_space_model(object, "object")$fitted_mean
}


#' @rdname gp_model_methods
#' @export
fitted.gaussianprocesses_heteroscedastic_model <- function(object, ...) {
  object <- .prepare_heteroscedastic_model(object, "object")
  .exact_fitted_values(object$mean_model)
}


#' @rdname gp_model_methods
#' @export
fitted.gaussianprocesses_latent_model <- function(object, ...) {
  .prepare_latent_model(object, "object")$mode
}


#' @rdname gp_model_methods
#' @export
residuals.gaussianprocesses_model <- function(object, ...) {
  object <- .prepare_gp_model(object, "object")
  object$y - .exact_fitted_values(object)
}


#' @rdname gp_model_methods
#' @export
residuals.gaussianprocesses_sparse_model <- function(object, ...) {
  object <- .prepare_sparse_gp_model(object, "object")
  object$y - fitted(object)
}


#' @rdname gp_model_methods
#' @export
residuals.gaussianprocesses_state_space_model <- function(object, ...) {
  object <- .prepare_state_space_model(object, "object")
  object$y - object$fitted_mean
}


#' @rdname gp_model_methods
#' @export
residuals.gaussianprocesses_heteroscedastic_model <- function(object, ...) {
  object <- .prepare_heteroscedastic_model(object, "object")
  object$y - fitted(object)
}


# summary ---------------------------------------------------------------------

.format_noise_summary <- function(model) {
  if (identical(model$noise_structure, "homoscedastic")) {
    return(paste0("homoscedastic, variance ", .format_noise_variance(model)))
  }

  paste0("known per observation, variances ", .format_noise_variance(model))
}


.new_model_summary <- function(
    model_type,
    x,
    kernel,
    mean,
    optimization,
    log_marginal_likelihood,
    noise,
    numerical,
    details = list()) {
  structure(
    list(
      model_type = model_type,
      n_observations = nrow(x),
      n_features = ncol(x),
      kernel = .format_kernel(kernel),
      mean = .format_mean(mean),
      noise = noise,
      log_marginal_likelihood = log_marginal_likelihood,
      optimization = optimization[
        intersect(
          c(
            "converged", "method", "gradient", "log_marginal_likelihood",
            "n_starts", "optima"
          ),
          names(optimization)
        )
      ],
      numerical = numerical,
      details = details
    ),
    class = "summary.gaussianprocesses_model"
  )
}


#' @rdname gp_model_methods
#' @export
summary.gaussianprocesses_model <- function(object, ...) {
  object <- .prepare_gp_model(object, "object")

  .new_model_summary(
    "Exact Gaussian-process regression",
    x = object$x,
    kernel = object$kernel,
    mean = object$mean,
    optimization = object$optimization,
    log_marginal_likelihood = log_marginal_likelihood(object),
    noise = .format_noise_summary(object),
    numerical = .format_jitter(object$jitter, .model_covariance_scale(object)),
    details = if (!is.null(object$time)) {
      list(time_unit = object$time$unit)
    } else {
      list()
    }
  )
}


#' @rdname gp_model_methods
#' @export
summary.gaussianprocesses_sparse_model <- function(object, ...) {
  object <- .prepare_sparse_gp_model(object, "object")

  .new_model_summary(
    paste0("Sparse Gaussian-process regression (", object$approximation, ")"),
    x = object$x,
    kernel = object$kernel,
    mean = object$mean,
    optimization = object$optimization,
    log_marginal_likelihood = object$log_marginal_likelihood,
    noise = .format_noise_summary(object),
    numerical = paste(
      format(signif(object$kuu_jitter, 3), trim = TRUE),
      "in the inducing-point covariance"
    ),
    details = list(
      inducing_points = object$n_inducing,
      trace_term = signif(object$trace_term, 4)
    )
  )
}


#' @rdname gp_model_methods
#' @export
summary.gaussianprocesses_state_space_model <- function(object, ...) {
  object <- .prepare_state_space_model(object, "object")

  .new_model_summary(
    "State-space Gaussian-process regression (Kalman smoothing)",
    x = object$x,
    kernel = object$kernel,
    mean = object$mean,
    optimization = object$optimization,
    log_marginal_likelihood = object$log_marginal_likelihood,
    noise = .format_noise_summary(object),
    numerical = "none (no covariance matrix is factorized)",
    details = c(
      list(state_dimension = object$state_dimension),
      if (!is.null(object$time)) list(time_unit = object$time$unit)
    )
  )
}


#' @rdname gp_model_methods
#' @export
summary.gaussianprocesses_heteroscedastic_model <- function(object, ...) {
  object <- .prepare_heteroscedastic_model(object, "object")
  noise <- object$training_noise_variance

  .new_model_summary(
    "Heteroscedastic Gaussian-process regression (plug-in)",
    x = object$x,
    kernel = object$kernel,
    mean = object$mean_model$mean,
    optimization = NULL,
    log_marginal_likelihood = NULL,
    noise = paste0(
      "input-dependent, training range [",
      format(signif(min(noise), 3), trim = TRUE),
      ", ",
      format(signif(max(noise), 3), trim = TRUE),
      "]"
    ),
    numerical = .format_jitter(
      object$mean_model$jitter,
      .model_covariance_scale(object$mean_model)
    ),
    details = list(
      noise_iterations = object$iterations,
      noise_converged = object$converged
    )
  )
}


#' @rdname gp_model_methods
#' @export
summary.gaussianprocesses_latent_model <- function(object, ...) {
  object <- .prepare_latent_model(object, "object")

  .new_model_summary(
    "Latent Gaussian-process (Laplace)",
    x = object$x,
    kernel = object$kernel,
    mean = object$mean,
    optimization = object$optimization,
    log_marginal_likelihood = object$log_marginal_likelihood,
    noise = NULL,
    numerical = format(signif(object$jitter, 3), trim = TRUE),
    details = list(
      likelihood = utils::capture.output(print(object$likelihood)),
      newton_iterations = object$laplace$iterations,
      newton_converged = object$laplace$converged
    )
  )
}


#' @rdname gp_model_methods
#' @method print summary.gaussianprocesses_model
#' @export
print.summary.gaussianprocesses_model <- function(x, ...) {
  cat(x$model_type, "model\n")
  cat(
    "  observations: ", x$n_observations,
    ", input dimensions: ", x$n_features, "\n",
    sep = ""
  )
  cat("  kernel: ", paste(trimws(x$kernel), collapse = " "), "\n", sep = "")
  cat("  mean: ", x$mean, "\n", sep = "")

  if (!is.null(x$noise)) {
    cat("  noise: ", x$noise, "\n", sep = "")
  }

  if (!is.null(x$log_marginal_likelihood)) {
    cat(
      "  log marginal likelihood: ",
      format(signif(x$log_marginal_likelihood, 7), trim = TRUE),
      "\n",
      sep = ""
    )
  }

  if (length(x$optimization) > 0L) {
    cat(
      "  optimization: ",
      if (isTRUE(x$optimization$converged)) "converged" else "not converged",
      " (", x$optimization$method, ", ", x$optimization$gradient,
      " gradient)\n",
      sep = ""
    )

    optima <- x$optimization$optima

    if (!is.null(optima)) {
      cat(
        "  starts: ", x$optimization$n_starts, ", reaching ",
        nrow(optima), " distinct ",
        if (nrow(optima) == 1L) "optimum" else "optima",
        sep = ""
      )

      if (nrow(optima) > 1L) {
        cat(
          " (log marginal likelihood gaps ",
          paste(
            vapply(signif(optima$gap, 3), format, character(1), trim = TRUE),
            collapse = ", "
          ),
          ")",
          sep = ""
        )
      }

      cat("\n")
    }
  }

  for (name in names(x$details)) {
    cat(
      "  ", gsub("_", " ", name), ": ", format(x$details[[name]]), "\n",
      sep = ""
    )
  }

  cat("  numerical jitter: ", x$numerical, "\n", sep = "")
  invisible(x)
}
