# Posterior prediction -------------------------------------------------------

.validate_prediction_inputs <- function(model, x) {
  x_matrix <- .as_input_matrix(x, name = "x")

  if (ncol(x_matrix) != model$n_features) {
    .gp_abort(
      sprintf(
        paste0(
          "'x' must have %d input dimension(s) to match the fitted model. ",
          "Received %d."
        ),
        model$n_features,
        ncol(x_matrix)
      ),
      "gaussianprocesses_validation_error"
    )
  }

  x_matrix
}


.sanitize_posterior_variance <- function(
    variance,
    reference_scale,
    tolerance = .default_tolerance()) {
  .assert_numeric_vector(variance, name = "variance")
  .assert_numeric_scalar(reference_scale, "reference_scale", lower = 0)
  .assert_numeric_scalar(tolerance, "tolerance", lower = 0)

  scale <- max(1, reference_scale)
  lower_limit <- -tolerance * scale

  if (any(variance < lower_limit)) {
    .gp_abort(
      paste0(
        "Posterior variance is materially negative. ",
        "This indicates numerical instability or an inconsistent covariance."
      ),
      "gaussianprocesses_numerical_error"
    )
  }

  pmax(variance, 0)
}


.validate_interval_level <- function(interval_level) {
  .assert_numeric_scalar(
    interval_level,
    "interval_level",
    lower = 0,
    lower_open = TRUE
  )

  if (interval_level >= 1) {
    .gp_abort(
      "'interval_level' must be strictly less than 1.",
      "gaussianprocesses_validation_error"
    )
  }

  invisible(interval_level)
}


.resolve_prediction_noise_variance <- function(
    model,
    observation_noise_variance,
    n_predictions) {
  if (is.null(observation_noise_variance)) {
    if (length(model$noise_variance) == 1L) {
      return(rep(model$noise_variance, n_predictions))
    }

    .gp_abort(
      paste0(
        "The fitted model used observation-specific training noise. ",
        "Supply 'observation_noise_variance' for the prediction inputs, ",
        "or use a model that predicts input-dependent noise."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (!is.numeric(observation_noise_variance) ||
      !is.null(dim(observation_noise_variance)) ||
      length(observation_noise_variance) == 0L ||
      any(!is.finite(observation_noise_variance)) ||
      any(observation_noise_variance < 0)) {
    .gp_abort(
      paste0(
        "'observation_noise_variance' must contain finite non-negative ",
        "numeric values."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (!(length(observation_noise_variance) %in%
      c(1L, n_predictions))) {
    .gp_abort(
      sprintf(
        paste0(
          "'observation_noise_variance' must contain either one value or ",
          "one value per prediction input. Received %d value(s) for %d ",
          "prediction point(s)."
        ),
        length(observation_noise_variance),
        n_predictions
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (length(observation_noise_variance) == 1L) {
    return(
      rep(as.numeric(observation_noise_variance), n_predictions)
    )
  }

  as.numeric(observation_noise_variance)
}


# Marginal prediction processes the prediction inputs in blocks so that the
# n x b cross-covariance block holds at most this many elements (8 MiB of
# doubles). Peak memory then depends on the training size n, not on the
# number of prediction inputs. Kernel evaluation creates several temporaries
# of this size, and for n = 500 larger blocks were no faster.
.prediction_block_elements <- 2^20


.prediction_block_size <- function(
    block_size,
    n_training,
    n_predictions) {
  if (is.null(block_size)) {
    block_size <- floor(.prediction_block_elements / n_training)
  } else {
    block_size <- .assert_positive_integer(block_size, "block_size")
  }

  as.integer(max(1, min(block_size, n_predictions)))
}


# Posterior mean offsets k(X, x*)' alpha and explained variances
# sum_i V_ij^2 with V = R^{-T} k(X, x*), one block of prediction inputs at a
# time. The prior covariance between prediction inputs is never formed. For
# marginalized mean coefficients, `coefficient_variance` is the diagonal of
# R_* A^{-1} R_*' (R/mean-coefficients.R), and zero otherwise.
.blocked_marginal_posterior <- function(model, prediction_x, block_size) {
  n_predictions <- nrow(prediction_x)
  mean_offset <- numeric(n_predictions)
  explained_variance <- numeric(n_predictions)
  coefficient_variance <- numeric(n_predictions)
  marginalized <- .marginalized_mean(model)

  if (marginalized) {
    whitened_basis <- .whitened_basis(model)
  }

  for (start in seq.int(1L, n_predictions, by = block_size)) {
    rows <- start:min(start + block_size - 1L, n_predictions)
    cross_covariance <- .model_cross_covariance(
      model,
      prediction_x[rows, , drop = FALSE]
    )

    mean_offset[rows] <- drop(crossprod(cross_covariance, model$alpha))

    # Solves R' V = k(X, x*) without forming the transposed factor.
    projected <- backsolve(
      model$cholesky,
      cross_covariance,
      upper.tri = TRUE,
      transpose = TRUE
    )
    explained_variance[rows] <- colSums(projected^2)

    if (marginalized) {
      uncertainty <- .coefficient_uncertainty(
        model,
        projected,
        .mean_basis(model$mean, prediction_x[rows, , drop = FALSE]),
        whitened_basis
      )
      coefficient_variance[rows] <- colSums(uncertainty^2)
    }
  }

  list(
    mean_offset = mean_offset,
    explained_variance = explained_variance,
    coefficient_variance = coefficient_variance
  )
}


#' Predict from an exact Gaussian-process regression model
#'
#' Computes the posterior mean and uncertainty at new input locations. Latent
#' function uncertainty and future-observation uncertainty are reported
#' separately. The stored Cholesky factor from model fitting is reused; no
#' covariance-matrix inverse is formed.
#'
#' @param model A fitted `gaussianprocesses_model`, or a state-space model
#'   from [fit_time_series_gp()], which is predicted by Kalman smoothing at
#'   inputs on its time-index scale (`block_size` does not apply; memory is
#'   linear in the number of inputs).
#' @param x Numeric vector or matrix of prediction inputs. Matrix rows are
#'   observations and columns must match the fitted input dimension.
#' @param interval_level Probability level for latent credible intervals and
#'   noisy-observation prediction intervals.
#' @param include_covariance If `TRUE`, return the full latent and observation
#'   posterior covariance matrices. Otherwise only marginal variances are
#'   returned.
#' @param variance_tolerance Relative tolerance used when deciding whether a
#'   small negative posterior variance is floating-point error or a material
#'   numerical failure.
#' @param observation_noise_variance Optional non-negative scalar or vector
#'   giving observation-noise variance at the prediction inputs. It defaults to
#'   the fitted scalar noise variance for homoscedastic models. Models fitted
#'   with observation-specific training noise require this argument.
#' @param block_size Number of prediction inputs processed together when
#'   `include_covariance = FALSE`. `NULL` (the default) chooses the largest
#'   block whose training-by-block cross-covariance holds at most \eqn{2^{20}}
#'   values (8 MiB). Results do not depend on the block size. It is ignored
#'   when `include_covariance = TRUE`.
#'
#' @details
#' With \eqn{n} training and \eqn{m} prediction inputs, \eqn{C = R^\top R}
#' the factorized training covariance, and
#' \eqn{V = R^{-\top} K(X, X_*)}, the latent posterior variance at
#' prediction input \eqn{j} is
#' \deqn{\mathrm{Var}[f(x_{*j}) \mid y] = k(x_{*j}, x_{*j}) -
#' \sum_i V_{ij}^2.}
#' This is exact. The default marginal prediction evaluates only the
#' cross-covariance \eqn{K(X, X_*)}, in blocks of `block_size` prediction
#' inputs, and the prior variances \eqn{k(x_{*j}, x_{*j})} with
#' [kernel_diagonal()]. It takes \eqn{O(n^2 m)} time and \eqn{O(n b)} memory
#' for blocks of \eqn{b} inputs, beyond the \eqn{O(n^2)} fitted factor.
#'
#' `include_covariance = TRUE` also forms the \eqn{m \times m} prior
#' covariance and the latent and observation posterior covariances, which
#' takes \eqn{O(n^2 m + n m^2)} time and \eqn{O(m^2)} memory. Each
#' \eqn{m \times m} matrix takes \eqn{8 m^2} bytes, 512 MB for 8000 inputs,
#' and several are held at once. It is opt-in because pointwise means,
#' variances, and intervals do not need it; joint quantities such as
#' function draws ([sample_gp_posterior()]) do.
#'
#' A parametric mean contributes \eqn{h(x_*)^\top \beta} with its fixed,
#' estimated, or posterior-mean coefficients. Estimated coefficients are
#' plugged in. Marginalized coefficients add their uncertainty,
#' \eqn{R_* (B^{-1} + H^\top C^{-1} H)^{-1} R_*^\top} with
#' \eqn{R_* = H_* - V^\top R^{-\top} H}, to the latent variances and
#' covariances ([mean_coefficients]).
#'
#' @return An object of class `gaussianprocesses_prediction`. Its fields
#'   are shared by every prediction in the package: [predict_sparse_gp()],
#'   [predict_heteroscedastic_gp()], and [forecast_gp()] return the same
#'   fields with the same meaning, and may add their own.
#'   \describe{
#'     \item{`x`}{The prediction inputs as a matrix.}
#'     \item{`mean`}{The posterior mean of the latent function.}
#'     \item{`latent_variance`, `latent_sd`}{The posterior variance and
#'     standard deviation of the latent function.}
#'     \item{`observation_noise_variance`}{The noise variance of a new
#'     observation.}
#'     \item{`observation_variance`, `observation_sd`}{The variance and
#'     standard deviation of a new noisy observation:
#'     `latent_variance + observation_noise_variance`.}
#'     \item{`interval_level`, `latent_interval`,
#'     `prediction_interval`}{The interval level and the intervals for the
#'     latent function and for a new observation.}
#'     \item{`latent_covariance`, `observation_covariance`}{The full
#'     covariance matrices, only with `include_covariance = TRUE`.}
#'   }
#' @examples
#' x <- seq(-2, 2, length.out = 15)
#' model <- fit_gp(
#'   x,
#'   sin(2 * x),
#'   kernel = rbf_kernel(length_scale = 0.6),
#'   noise_variance = 0.01
#' )
#'
#' prediction <- predict_gp(model, c(0, 3))
#'
#' # Latent-function and noisy-observation uncertainty are reported separately;
#' # both grow away from the data.
#' prediction$latent_sd
#' prediction$observation_sd
#' prediction$prediction_interval
#' @template stable
#' @export
predict_gp <- function(
    model,
    x,
    interval_level = 0.95,
    include_covariance = FALSE,
    variance_tolerance = sqrt(.Machine$double.eps),
    observation_noise_variance = NULL,
    block_size = NULL) {
  if (inherits(model, "gaussianprocesses_state_space_model")) {
    return(.predict_state_space(
      .prepare_state_space_model(model),
      x,
      interval_level = interval_level,
      include_covariance = include_covariance,
      variance_tolerance = variance_tolerance,
      observation_noise_variance = observation_noise_variance
    ))
  }

  model <- .prepare_gp_model(model)
  prediction_x <- .validate_prediction_inputs(model, x)
  .validate_interval_level(interval_level)
  .assert_numeric_scalar(
    variance_tolerance,
    "variance_tolerance",
    lower = 0
  )

  if (!is.logical(include_covariance) ||
      length(include_covariance) != 1L ||
      is.na(include_covariance)) {
    .gp_abort(
      "'include_covariance' must be TRUE or FALSE.",
      "gaussianprocesses_validation_error"
    )
  }

  prediction_noise_variance <- .resolve_prediction_noise_variance(
    model,
    observation_noise_variance,
    n_predictions = nrow(prediction_x)
  )

  block_size <- .prediction_block_size(
    block_size,
    n_training = model$n_observations,
    n_predictions = nrow(prediction_x)
  )

  prior_mean <- .model_mean(model, prediction_x)

  if (include_covariance) {
    cross_covariance <- .model_cross_covariance(model, prediction_x)
    prior_covariance <- evaluate_kernel(
      model$kernel,
      prediction_x
    )

    posterior_mean <- prior_mean +
      drop(crossprod(cross_covariance, model$alpha))

    projected <- forwardsolve(
      t(model$cholesky),
      cross_covariance,
      upper.tri = FALSE,
      transpose = FALSE
    )

    latent_covariance <- prior_covariance - crossprod(projected)

    if (.marginalized_mean(model)) {
      latent_covariance <- latent_covariance + crossprod(
        .coefficient_uncertainty(
          model,
          projected,
          .mean_basis(model$mean, prediction_x),
          .whitened_basis(model)
        )
      )
    }

    latent_covariance <- 0.5 * (
      latent_covariance + t(latent_covariance)
    )

    prior_reference <- max(abs(diag(prior_covariance)))
    latent_variance <- .sanitize_posterior_variance(
      diag(latent_covariance),
      reference_scale = prior_reference,
      tolerance = variance_tolerance
    )

    diag(latent_covariance) <- latent_variance
  } else {
    marginal <- .blocked_marginal_posterior(model, prediction_x, block_size)
    prior_variance <- kernel_diagonal(model$kernel, prediction_x)

    posterior_mean <- prior_mean + marginal$mean_offset
    latent_variance <- .sanitize_posterior_variance(
      prior_variance - marginal$explained_variance +
        marginal$coefficient_variance,
      reference_scale = max(abs(prior_variance)),
      tolerance = variance_tolerance
    )
  }

  observation_variance <- latent_variance +
    prediction_noise_variance
  latent_sd <- sqrt(latent_variance)
  observation_sd <- sqrt(observation_variance)

  quantile <- stats::qnorm((1 + interval_level) / 2)

  latent_interval <- cbind(
    lower = posterior_mean - quantile * latent_sd,
    upper = posterior_mean + quantile * latent_sd
  )
  prediction_interval <- cbind(
    lower = posterior_mean - quantile * observation_sd,
    upper = posterior_mean + quantile * observation_sd
  )

  result <- list(
    x = prediction_x,
    mean = posterior_mean,
    latent_variance = latent_variance,
    observation_noise_variance = prediction_noise_variance,
    observation_variance = observation_variance,
    latent_sd = latent_sd,
    observation_sd = observation_sd,
    interval_level = interval_level,
    latent_interval = latent_interval,
    prediction_interval = prediction_interval
  )

  if (include_covariance) {
    result$latent_covariance <- latent_covariance
    result$observation_covariance <- latent_covariance +
      diag(prediction_noise_variance, nrow(prediction_x))
  }

  structure(
    result,
    class = "gaussianprocesses_prediction"
  )
}


#' @export
print.gaussianprocesses_prediction <- function(x, ...) {
  cat("Gaussian-process posterior prediction\n")
  cat("  prediction points: ", length(x$mean), "\n", sep = "")
  cat(
    "  interval level: ",
    format(100 * x$interval_level, trim = TRUE),
    "%\n",
    sep = ""
  )
  cat(
    "  mean range: [",
    format(signif(min(x$mean), 6), trim = TRUE),
    ", ",
    format(signif(max(x$mean), 6), trim = TRUE),
    "]\n",
    sep = ""
  )
  invisible(x)
}
