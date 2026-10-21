# Prediction and sampling for latent Gaussian-process models -----------------

# Latent posterior mean and covariance terms at prediction inputs (Rasmussen
# and Williams, 2006, Algorithm 3.2): with k_* = k(X, x_*),
#   E[f_*] = m(x_*) + k_*' a,
#   Var[f_*] = k(x_*, x_*) - v' v,   v = R^{-T} W^(1/2) k_*,
# where R'R = B. Marginal variances are computed in blocks of prediction
# inputs, as in predict_gp().
.latent_posterior <- function(model, prediction_x, include_covariance,
                              block_size) {
  root <- sqrt(model$curvature)
  prior_mean <- .latent_mean_values(model, prediction_x)

  whiten <- function(cross_covariance) {
    backsolve(
      model$b_cholesky,
      root * cross_covariance,
      upper.tri = TRUE,
      transpose = TRUE
    )
  }

  if (include_covariance) {
    cross_covariance <- .latent_prior_covariance(
      model$kernel,
      model$mean,
      model$x,
      prediction_x
    )
    whitened <- whiten(cross_covariance)
    prior_covariance <- .latent_prior_covariance(
      model$kernel,
      model$mean,
      prediction_x
    )

    return(list(
      mean = prior_mean + drop(crossprod(cross_covariance, model$alpha)),
      prior_variance = diag(prior_covariance),
      covariance = prior_covariance - crossprod(whitened)
    ))
  }

  n_predictions <- nrow(prediction_x)
  mean_offset <- numeric(n_predictions)
  explained_variance <- numeric(n_predictions)

  for (start in seq.int(1L, n_predictions, by = block_size)) {
    rows <- start:min(start + block_size - 1L, n_predictions)
    cross_covariance <- .latent_prior_covariance(
      model$kernel,
      model$mean,
      model$x,
      prediction_x[rows, , drop = FALSE]
    )
    mean_offset[rows] <- drop(crossprod(cross_covariance, model$alpha))
    explained_variance[rows] <- colSums(whiten(cross_covariance)^2)
  }

  prior_variance <- .latent_prior_variance(
    model$kernel,
    model$mean,
    prediction_x
  )

  list(
    mean = prior_mean + mean_offset,
    prior_variance = prior_variance,
    variance = prior_variance - explained_variance
  )
}


#' Predict from a latent Gaussian-process model
#'
#' Computes the approximate posterior of the latent function at new inputs
#' and the predictive distribution of new responses through the likelihood:
#' class probabilities for Bernoulli likelihoods and predictive means and
#' variances for every likelihood.
#'
#' @param model A fitted `gaussianprocesses_latent_model` from
#'   [fit_latent_gp()] or [optimize_latent_gp()].
#' @param x Numeric vector or matrix of prediction inputs.
#' @param interval_level Probability level of the latent credible intervals.
#' @param include_covariance If `TRUE`, also return the full latent
#'   posterior covariance matrix.
#' @param variance_tolerance Relative tolerance used when deciding whether a
#'   small negative posterior variance is floating-point error or a material
#'   numerical failure.
#' @param exposure Optional positive exposure of the new responses for
#'   Poisson likelihoods: one value or one per prediction input. It defaults
#'   to 1.
#' @param quadrature_order Number of Gauss-Hermite nodes for predictive
#'   integrals without a closed form ([likelihood_predictive()]).
#' @param block_size Number of prediction inputs processed together when
#'   `include_covariance = FALSE`, as in [predict_gp()].
#'
#' @details
#' Under the Laplace approximation the latent values at new inputs are
#' Gaussian (Rasmussen and Williams, 2006, Algorithm 3.2), with
#' \deqn{\mathrm{E}[f_*] = m(x_*) + k_*^\top \nabla \log p(y \mid \hat f),
#' \qquad \mathrm{Var}[f_*] = k(x_*, x_*) - v^\top v,}
#' where \eqn{k_* = k(X, x_*)}, \eqn{v = R^{-\top} W^{1/2} k_*}, and
#' \eqn{R^\top R = B}. The predictive distribution of a new response
#' integrates the likelihood against this Gaussian with
#' [likelihood_predictive()].
#'
#' @section Class probabilities:
#' For a Bernoulli likelihood, `probability` is \eqn{\Pr(y_* = 1 \mid y)}
#' for the response coded 1: the second level of a factor response
#' (`model$response_levels`), or `TRUE`. It averages the link over the
#' latent posterior,
#' \deqn{\Pr(y_* = 1 \mid y) = \int \pi(f_*) N(f_* \mid \mu_*, \sigma_*^2)
#' df_*,}
#' which for the probit link is \eqn{\Phi(\mu_* / \sqrt{1 + \sigma_*^2})}
#' exactly and for the logistic link is computed by Gauss-Hermite
#' quadrature. It is not \eqn{\pi(\mu_*)}: the plug-in probability ignores
#' the latent uncertainty and is too confident, by up to 0.04 at
#' \eqn{\sigma_*^2 = 1} and 0.11 at \eqn{\sigma_*^2 = 4}. MacKay's
#' approximation \eqn{\sigma(\kappa \mu_*)} with
#' \eqn{\kappa = (1 + \pi \sigma_*^2 / 8)^{-1/2}} is closer but still off by
#' up to 0.005 at \eqn{\sigma_*^2 = 1} and 0.01 at \eqn{\sigma_*^2 = 4}, so
#' it is not used.
#'
#' `probability_interval` describes the uncertainty of the class probability
#' \eqn{\pi(f_*)} itself: because the link is increasing, the
#' `interval_level` quantiles of \eqn{\pi(f_*)} are the link applied to the
#' latent interval, exactly under the Gaussian approximation. Joint
#' uncertainty, for example of a difference of probabilities, comes from
#' the `probability` draws of [sample_gp_posterior()].
#'
#' Multi-class (softmax) classification is not supported.
#'
#' @section Counts:
#' For a Poisson likelihood the rate \eqn{\lambda_* = E_* e^{f_*}} is
#' lognormal: its mean is \eqn{E_* e^{\mu_* + \sigma_*^2 / 2}}
#' (`response_mean`, also the mean of the count) and its quantiles are
#' \eqn{E_* e^{\mu_* \pm z \sigma_*}} exactly, because the exponential is
#' increasing (`rate_interval`). A new count has the Poisson-lognormal
#' predictive distribution; `prediction_interval` gives its equal-tailed
#' quantiles at `interval_level`, computed by quadrature
#' ([gp_count_scores()]). Counts are discrete, so these intervals cover at
#' least `interval_level`: their coverage is conservative.
#'
#' @return An object of class `gaussianprocesses_latent_prediction` with the
#'   shared latent fields of [predict_gp()], `x`, `mean`, `latent_variance`,
#'   `latent_sd`, `interval_level`, and `latent_interval` (and
#'   `latent_covariance` with `include_covariance = TRUE`), and the
#'   predictive distribution of a new response: `response_mean`,
#'   `response_variance`, `probability` \eqn{\Pr(y_* = 1)} and
#'   `probability_interval` for Bernoulli likelihoods, `rate_interval`,
#'   `prediction_interval`, and `exposure` for Poisson likelihoods, and
#'   `prediction_interval` for Gaussian likelihoods. A non-Gaussian response
#'   has no observation-noise variance, so the observation fields of
#'   [predict_gp()] are replaced by these; `prediction_interval` has its
#'   meaning there, an interval for a new response.
#' @seealso [fit_latent_gp()], [sample_gp_posterior()].
#' @examples
#' x <- seq(-3, 3, length.out = 30)
#' y <- as.numeric(sin(x) + 0.3 * cos(5 * x) > 0)
#' model <- fit_latent_gp(x, y, rbf_kernel(variance = 4), bernoulli_likelihood())
#'
#' prediction <- predict_latent_gp(model, c(-2, 0, 2, 6))
#'
#' # Far from the data the latent variance grows and the probability moves
#' # towards 1/2, unlike the logistic function of the latent mean.
#' cbind(prediction$probability, plogis(prediction$mean))
#' @template experimental
#' @export
predict_latent_gp <- function(
    model,
    x,
    interval_level = 0.95,
    include_covariance = FALSE,
    variance_tolerance = sqrt(.Machine$double.eps),
    exposure = NULL,
    quadrature_order = 40L,
    block_size = NULL) {
  model <- .prepare_latent_model(model)
  prediction_x <- .validate_prediction_inputs(model, x)
  .validate_interval_level(interval_level)
  .assert_numeric_scalar(variance_tolerance, "variance_tolerance", lower = 0)

  if (!is.logical(include_covariance) ||
      length(include_covariance) != 1L ||
      is.na(include_covariance)) {
    .gp_abort(
      "'include_covariance' must be TRUE or FALSE.",
      "gaussianprocesses_validation_error"
    )
  }

  n_predictions <- nrow(prediction_x)
  resolved_exposure <- .likelihood_exposure(
    model$likelihood,
    exposure,
    n_predictions
  )
  block_size <- .prediction_block_size(
    block_size,
    n_training = model$n_observations,
    n_predictions = n_predictions
  )
  posterior <- .latent_posterior(
    model,
    prediction_x,
    include_covariance,
    block_size
  )

  if (include_covariance) {
    latent_covariance <- 0.5 * (posterior$covariance + t(posterior$covariance))
    posterior$variance <- diag(latent_covariance)
  }

  latent_variance <- .sanitize_posterior_variance(
    posterior$variance,
    reference_scale = max(abs(posterior$prior_variance)),
    tolerance = variance_tolerance
  )
  latent_sd <- sqrt(latent_variance)
  quantile <- stats::qnorm((1 + interval_level) / 2)
  response <- likelihood_predictive(
    model$likelihood,
    posterior$mean,
    latent_variance,
    exposure = if (identical(model$likelihood$family, "poisson")) {
      resolved_exposure
    },
    quadrature_order = quadrature_order
  )

  latent_interval <- cbind(
    lower = posterior$mean - quantile * latent_sd,
    upper = posterior$mean + quantile * latent_sd
  )
  family <- model$likelihood$family
  bernoulli <- identical(family, "bernoulli")
  poisson <- identical(family, "poisson")
  entry <- .likelihood_entry(model$likelihood)

  result <- list(
    x = prediction_x,
    mean = posterior$mean,
    latent_variance = latent_variance,
    latent_sd = latent_sd,
    interval_level = interval_level,
    latent_interval = latent_interval,
    response_mean = response$mean,
    response_variance = response$variance,
    probability = response$probability,
    # The link is increasing, so it maps latent quantiles to quantiles of
    # the class probability.
    probability_interval = if (bernoulli) {
      entry$conditional_mean(latent_interval, list(), 1)
    },
    rate_interval = if (poisson) {
      entry$conditional_mean(latent_interval, list(), resolved_exposure)
    },
    prediction_interval = if (poisson) {
      .count_prediction_interval(
        interval_level,
        posterior$mean,
        latent_variance,
        resolved_exposure,
        quadrature_order
      )
    } else if (identical(family, "gaussian")) {
      response_sd <- sqrt(response$variance)
      cbind(
        lower = response$mean - quantile * response_sd,
        upper = response$mean + quantile * response_sd
      )
    },
    exposure = if (poisson) resolved_exposure
  )

  if (include_covariance) {
    diag(latent_covariance) <- latent_variance
    result$latent_covariance <- latent_covariance
  }

  structure(
    result[!vapply(result, is.null, logical(1))],
    class = "gaussianprocesses_latent_prediction"
  )
}


#' @export
print.gaussianprocesses_latent_prediction <- function(x, ...) {
  cat("Latent Gaussian-process prediction (Laplace approximation)\n")
  cat("  prediction points: ", length(x$mean), "\n", sep = "")
  cat(
    "  latent mean range: [",
    format(signif(min(x$mean), 6), trim = TRUE),
    ", ",
    format(signif(max(x$mean), 6), trim = TRUE),
    "]\n",
    sep = ""
  )

  if (!is.null(x$probability)) {
    cat(
      "  probability range: [",
      format(signif(min(x$probability), 4), trim = TRUE),
      ", ",
      format(signif(max(x$probability), 4), trim = TRUE),
      "]\n",
      sep = ""
    )
  } else {
    cat(
      "  response mean range: [",
      format(signif(min(x$response_mean), 6), trim = TRUE),
      ", ",
      format(signif(max(x$response_mean), 6), trim = TRUE),
      "]\n",
      sep = ""
    )
  }

  invisible(x)
}


#' @template stable
#' @rdname sample_gp_posterior
#' @param exposure For latent Poisson models, optional exposure of the
#'   observation draws, as in [predict_latent_gp()].
#' @export
sample_gp_posterior.gaussianprocesses_latent_model <- function(
    model,
    x,
    n_draws = 1L,
    type = c("latent", "observation"),
    seed = NULL,
    exposure = NULL,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps),
    ...) {
  model <- .prepare_latent_model(model)
  type <- match.arg(type)
  n_draws <- .assert_draw_count(n_draws)
  .validate_seed(seed)

  prediction <- predict_latent_gp(
    model,
    x,
    include_covariance = TRUE,
    exposure = exposure
  )
  decomposition <- .cholesky_for_sampling(
    prediction$latent_covariance,
    initial_jitter,
    fallback_jitter,
    jitter_multiplier,
    max_attempts,
    symmetry_tolerance
  )
  entry <- .likelihood_entry(model$likelihood)
  draw_exposure <- if (is.null(prediction$exposure)) {
    rep(1, length(prediction$mean))
  } else {
    prediction$exposure
  }

  draws <- .with_local_seed(
    seed,
    function() {
      latent <- .gaussian_draws(
        prediction$mean,
        decomposition$factor,
        n_draws
      )
      observation <- if (identical(type, "observation")) {
        entry$simulate(latent, model$likelihood$parameters, draw_exposure)
      }

      list(latent = latent, observation = observation)
    }
  )

  samples <- .new_gp_samples(
    x = prediction$x,
    distribution = "posterior",
    model = "Laplace approximation",
    latent = draws$latent,
    observation = draws$observation,
    mean = prediction$mean,
    latent_covariance = prediction$latent_covariance,
    observation_noise_variance = NULL,
    decomposition = decomposition,
    seed = seed
  )

  if (identical(model$likelihood$family, "bernoulli")) {
    samples$probability <- entry$conditional_mean(draws$latent, list(), 1)
  }

  if (identical(model$likelihood$family, "poisson")) {
    samples$rate <- entry$conditional_mean(draws$latent, list(), draw_exposure)
  }

  samples
}
