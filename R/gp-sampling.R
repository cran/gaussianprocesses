# Prior and posterior function draws ------------------------------------------
#
# A draw from N(mu, Sigma) is mu + R' z with Sigma = R'R and z ~ N(0, I).
# R comes from the package's stable Cholesky routine, so nearly singular
# covariances follow the usual jitter policy and the jitter is reported. No
# covariance matrix is inverted.

.assert_draw_count <- function(n_draws) {
  .assert_positive_integer(n_draws, "n_draws")
}


.cholesky_for_sampling <- function(
    covariance,
    initial_jitter,
    fallback_jitter,
    jitter_multiplier,
    max_attempts,
    tolerance) {
  .stable_cholesky(
    covariance,
    initial_jitter = initial_jitter,
    fallback_jitter = fallback_jitter,
    jitter_multiplier = jitter_multiplier,
    max_attempts = max_attempts,
    tolerance = tolerance
  )
}


# One draw per column: mean + R' Z with Z an n x n_draws standard normal
# matrix. Uses the same arithmetic as simulate_gp_data(), so a single draw
# with the same seed reproduces its latent values.
.gaussian_draws <- function(mean, factor, n_draws) {
  standard <- matrix(
    stats::rnorm(length(mean) * n_draws),
    nrow = length(mean)
  )

  mean + t(factor) %*% standard
}


.new_gp_samples <- function(
    x,
    distribution,
    model,
    latent,
    observation,
    mean,
    latent_covariance,
    observation_noise_variance,
    decomposition,
    seed) {
  structure(
    list(
      x = x,
      distribution = distribution,
      model = model,
      n_draws = ncol(latent),
      latent = latent,
      observation = observation,
      mean = mean,
      latent_covariance = latent_covariance,
      observation_noise_variance = observation_noise_variance,
      jitter = decomposition$jitter,
      relative_jitter = decomposition$relative_jitter,
      covariance_scale = decomposition$covariance_scale,
      cholesky_attempts = decomposition$attempts,
      seed = seed
    ),
    class = "gaussianprocesses_samples"
  )
}


#' Draw functions from a Gaussian-process prior
#'
#' Draws latent function values \eqn{f(x) \sim \mathcal{N}(m(x), K(x, x))}
#' at the supplied inputs. Use [simulate_gp_data()] for a single prior draw
#' with added observation noise.
#'
#' @param x Numeric vector or matrix of inputs. Matrix rows are input points.
#' @param kernel A Gaussian-process kernel specification, including ARD and
#'   composite kernels.
#' @param mean Gaussian-process mean specification. Its coefficients must
#'   be fixed or have a Gaussian [coefficient_prior()], whose uncertainty the
#'   draws include.
#' @param n_draws Number of independent function draws.
#' @param seed Optional non-negative integer seed. The caller's random state
#'   is restored afterwards.
#' @param initial_jitter Non-negative jitter tried first, relative to the
#'   covariance scale (the mean of its diagonal).
#' @param fallback_jitter Positive relative jitter tried after the first
#'   failed factorization.
#' @param jitter_multiplier Multiplicative jitter escalation factor.
#' @param max_attempts Maximum Cholesky attempts.
#' @param symmetry_tolerance Relative tolerance for checking that the
#'   covariance matrix is symmetric.
#'
#' @details
#' Each draw is \eqn{m(x) + R^\top z} with \eqn{K(x, x) = R^\top R} and
#' \eqn{z \sim \mathcal{N}(0, I)}. Prior covariances on closely spaced inputs
#' are often numerically singular; the Cholesky factor then uses the package's
#' jitter policy, and the jitter is reported in the result. A single draw
#' with the same seed equals the latent values of [simulate_gp_data()].
#'
#' @return An object of class `gaussianprocesses_samples`. Its `latent`
#'   element is a matrix with one row per input point and one column per draw.
#'   It also contains the prior `mean`, the `latent_covariance` that was
#'   sampled, the numerical `jitter`, the `relative_jitter` and
#'   `covariance_scale` it was computed from (jitter is relative jitter times
#'   the mean diagonal of the covariance), and the `cholesky_attempts`.
#'
#' @examples
#' x <- seq(0, 5, length.out = 100)
#' draws <- sample_gp_prior(x, rbf_kernel(length_scale = 0.8), n_draws = 3, seed = 1)
#' draws
#'
#' matplot(x, draws$latent, type = "l", lty = 1, ylab = "f(x)")
#' @template stable
#' @export
sample_gp_prior <- function(
    x,
    kernel,
    mean = zero_mean(),
    n_draws = 1L,
    seed = NULL,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps)) {
  x_matrix <- .as_input_matrix(x, name = "x")
  .assert_kernel(kernel)
  .assert_kernel_dimensions(kernel, ncol(x_matrix))
  .assert_mean_function(mean)
  n_draws <- .assert_draw_count(n_draws)
  .validate_seed(seed)

  prior_mean <- evaluate_mean(mean, x_matrix)
  covariance <- evaluate_kernel(kernel, x_matrix)
  coefficient_covariance <- .mean_prior_covariance(mean, x_matrix)

  # Coefficients with a Gaussian prior N(b, B) add h(x)' B h(x').
  if (!is.null(coefficient_covariance)) {
    covariance <- covariance + coefficient_covariance
  }
  decomposition <- .cholesky_for_sampling(
    covariance,
    initial_jitter,
    fallback_jitter,
    jitter_multiplier,
    max_attempts,
    symmetry_tolerance
  )

  latent <- .with_local_seed(
    seed,
    function() .gaussian_draws(prior_mean, decomposition$factor, n_draws)
  )

  .new_gp_samples(
    x = x_matrix,
    distribution = "prior",
    model = "prior",
    latent = latent,
    observation = NULL,
    mean = prior_mean,
    latent_covariance = covariance,
    observation_noise_variance = NULL,
    decomposition = decomposition,
    seed = seed
  )
}


# Latent draws and, optionally, noisy observation draws y = f + e with
# e ~ N(0, diag(noise)), all within one seeded scope.
.posterior_draws <- function(
    prediction_x,
    distribution_model,
    mean,
    latent_covariance,
    observation_noise_variance,
    type,
    n_draws,
    seed,
    initial_jitter,
    fallback_jitter,
    jitter_multiplier,
    max_attempts,
    tolerance) {
  decomposition <- .cholesky_for_sampling(
    latent_covariance,
    initial_jitter,
    fallback_jitter,
    jitter_multiplier,
    max_attempts,
    tolerance
  )

  draws <- .with_local_seed(
    seed,
    function() {
      latent <- .gaussian_draws(mean, decomposition$factor, n_draws)
      observation <- if (identical(type, "observation")) {
        latent + sqrt(observation_noise_variance) *
          matrix(stats::rnorm(length(mean) * n_draws), nrow = length(mean))
      }

      list(latent = latent, observation = observation)
    }
  )

  .new_gp_samples(
    x = prediction_x,
    distribution = "posterior",
    model = distribution_model,
    latent = draws$latent,
    observation = draws$observation,
    mean = mean,
    latent_covariance = latent_covariance,
    observation_noise_variance = if (identical(type, "observation")) {
      observation_noise_variance
    },
    decomposition = decomposition,
    seed = seed
  )
}


#' Draw functions from a Gaussian-process posterior
#'
#' Draws latent function values from the posterior of a fitted model and,
#' optionally, noisy observations from the posterior predictive distribution.
#'
#' @param model A fitted `gaussianprocesses_model` (from [fit_gp()],
#'   [optimize_gp()], or [fit_time_series_gp()]), a
#'   `gaussianprocesses_heteroscedastic_model` (from
#'   [fit_heteroscedastic_gp()]), a `gaussianprocesses_latent_model` (from
#'   [fit_latent_gp()] or [optimize_latent_gp()]), or a
#'   `gaussianprocesses_sparse_model` (from [fit_sparse_gp()] or
#'   [optimize_sparse_gp()]).
#' @param x Numeric vector or matrix of prediction inputs. For models fitted
#'   with [fit_time_series_gp()], use the numeric time index.
#' @param n_draws Number of independent draws.
#' @param type `"latent"` draws the latent function \eqn{f(x)}. `"observation"`
#'   also draws noisy observations \eqn{y = f(x) + \varepsilon}; both are
#'   returned, with the observations generated from the same latent draws.
#' @param observation_noise_variance For exact and sparse models, optional
#'   non-negative scalar or vector of noise variances at `x`, used only when
#'   `type = "observation"`. It defaults to the fitted noise variance of
#'   homoscedastic models and is required for models fitted with
#'   observation-specific noise.
#' @param seed Optional non-negative integer seed. The caller's random state
#'   is restored afterwards.
#' @param initial_jitter Non-negative jitter tried first, relative to the
#'   covariance scale (the mean of its diagonal).
#' @param fallback_jitter Positive relative jitter tried after the first
#'   failed factorization.
#' @param jitter_multiplier Multiplicative jitter escalation factor.
#' @param max_attempts Maximum Cholesky attempts.
#' @param symmetry_tolerance Relative tolerance for checking that the
#'   covariance matrix is symmetric.
#' @param ... Further arguments passed to methods.
#'
#' @details
#' The posterior mean \eqn{\mu} and latent covariance \eqn{\Sigma} come from
#' [predict_gp()] (or [predict_heteroscedastic_gp()]), which uses triangular
#' solves with the fitted Cholesky factor; no covariance matrix is inverted.
#' Latent draws are \eqn{\mu + R^\top z} with \eqn{\Sigma = R^\top R} and
#' \eqn{z \sim \mathcal{N}(0, I)}. Posterior covariances are often nearly
#' singular, for example at closely spaced or repeated inputs; the factor then
#' uses the package's jitter policy, and the jitter is reported.
#'
#' Observation draws add independent noise,
#' \eqn{y = f + \varepsilon} with \eqn{\varepsilon \sim
#' \mathcal{N}(0, \mathrm{diag}(\sigma^2_*))}, so their covariance is
#' \eqn{\Sigma + \mathrm{diag}(\sigma^2_*)}.
#'
#' For heteroscedastic models, the latent draws come from the mean GP and the
#' noise variances \eqn{\sigma^2_*} are the plug-in estimates of
#' [predict_heteroscedastic_gp()]. The draws are therefore conditional on the
#' estimated noise function: uncertainty in the noise process is not
#' propagated, consistent with the approximation being a plug-in method rather
#' than full Bayesian inference.
#'
#' For latent models, the latent draws come from the Gaussian (Laplace)
#' approximation of [predict_latent_gp()], and observation draws are drawn
#' from the likelihood given each latent draw: Bernoulli, Poisson (with
#' `exposure`), or Gaussian responses. Their `observation_noise_variance` is
#' `NULL`. For Bernoulli likelihoods the samples also contain `probability`,
#' the class probability \eqn{\pi(f)} of each latent draw, whose mean over
#' draws is the predictive probability of [predict_latent_gp()]; for Poisson
#' likelihoods they contain `rate`, the rate \eqn{E e^f} of each latent
#' draw.
#'
#' For sparse models (FITC or VFE, [fit_sparse_gp()]), the draws use the
#' sparse posterior covariance
#' \eqn{K_{**} - Q_{**} + K_{*u} A^{-1} K_{u*}} of [predict_sparse_gp()] at
#' the \eqn{k} sampling inputs. It costs \eqn{O(k m^2 + k^2 m)} time for
#' \eqn{m} inducing points, and its factorization \eqn{O(k^3)} time and
#' \eqn{O(k^2)} memory, as for exact models. The number of inducing points
#' enters only through \eqn{m \times k} cross-covariances, so the limit is
#' the number of sampling inputs, a few thousand, rather than \eqn{m}.
#' Observation draws need `observation_noise_variance` for models fitted with
#' observation-specific noise.
#'
#' @return An object of class `gaussianprocesses_samples` with one row per
#'   input point and one column per draw in `latent` and, for
#'   `type = "observation"`, in `observation`. It also contains the posterior
#'   `mean`, the `latent_covariance` that was sampled, the
#'   `observation_noise_variance` used for observation draws, the numerical
#'   `jitter` with its `relative_jitter` and `covariance_scale`, and the
#'   `cholesky_attempts`.
#'
#' @examples
#' x <- c(-2, -1, 0, 1.5)
#' y <- sin(2 * x)
#' model <- fit_gp(x, y, rbf_kernel(length_scale = 0.7), noise_variance = 0.01)
#'
#' x_new <- seq(-3, 3, length.out = 100)
#' draws <- sample_gp_posterior(model, x_new, n_draws = 5, seed = 1)
#' draws
#'
#' matplot(x_new, draws$latent, type = "l", lty = 1, ylab = "f(x)")
#' points(x, y, pch = 19)
#'
#' # Noisy observation draws are generated from the same latent draws.
#' noisy <- sample_gp_posterior(
#'   model,
#'   x_new,
#'   n_draws = 1,
#'   type = "observation",
#'   seed = 2
#' )
#' head(cbind(latent = noisy$latent[, 1], observation = noisy$observation[, 1]))
#' @export
sample_gp_posterior <- function(model, x, n_draws = 1L, ...) {
  UseMethod("sample_gp_posterior")
}


#' @rdname sample_gp_posterior
#' @export
sample_gp_posterior.default <- function(model, x, n_draws = 1L, ...) {
  .gp_abort(
    paste0(
      "'model' must be a fitted exact, heteroscedastic, or latent ",
      "Gaussian-process model."
    ),
    "gaussianprocesses_validation_error"
  )
}


#' @rdname sample_gp_posterior
#' @export
sample_gp_posterior.gaussianprocesses_model <- function(
    model,
    x,
    n_draws = 1L,
    type = c("latent", "observation"),
    observation_noise_variance = NULL,
    seed = NULL,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps),
    ...) {
  model <- .prepare_gp_model(model)
  type <- match.arg(type)
  n_draws <- .assert_draw_count(n_draws)
  .validate_seed(seed)

  if (identical(type, "latent") && !is.null(observation_noise_variance)) {
    .gp_abort(
      paste0(
        "'observation_noise_variance' is used only for observation draws; ",
        "use type = \"observation\"."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  # Latent draws need no noise variance, even for models fitted with
  # observation-specific noise.
  prediction <- predict_gp(
    model,
    x,
    include_covariance = TRUE,
    observation_noise_variance = if (identical(type, "latent")) {
      0
    } else {
      observation_noise_variance
    }
  )

  .posterior_draws(
    prediction_x = prediction$x,
    distribution_model = "exact",
    mean = prediction$mean,
    latent_covariance = prediction$latent_covariance,
    observation_noise_variance = prediction$observation_noise_variance,
    type = type,
    n_draws = n_draws,
    seed = seed,
    initial_jitter = initial_jitter,
    fallback_jitter = fallback_jitter,
    jitter_multiplier = jitter_multiplier,
    max_attempts = max_attempts,
    tolerance = symmetry_tolerance
  )
}


#' @rdname sample_gp_posterior
#' @export
sample_gp_posterior.gaussianprocesses_heteroscedastic_model <- function(
    model,
    x,
    n_draws = 1L,
    type = c("latent", "observation"),
    seed = NULL,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps),
    ...) {
  model <- .prepare_heteroscedastic_model(model)
  type <- match.arg(type)
  n_draws <- .assert_draw_count(n_draws)
  .validate_seed(seed)

  prediction <- predict_heteroscedastic_gp(
    model,
    x,
    include_covariance = TRUE
  )

  .posterior_draws(
    prediction_x = prediction$x,
    distribution_model = "heteroscedastic plug-in",
    mean = prediction$mean,
    latent_covariance = prediction$latent_covariance,
    observation_noise_variance = prediction$observation_noise_variance,
    type = type,
    n_draws = n_draws,
    seed = seed,
    initial_jitter = initial_jitter,
    fallback_jitter = fallback_jitter,
    jitter_multiplier = jitter_multiplier,
    max_attempts = max_attempts,
    tolerance = symmetry_tolerance
  )
}


#' @export
print.gaussianprocesses_samples <- function(x, ...) {
  cat("Gaussian-process samples\n")
  cat(
    "  distribution: ",
    x$distribution,
    if (!identical(x$model, x$distribution)) paste0(" (", x$model, ")"),
    "\n",
    sep = ""
  )
  cat(
    "  draws: ",
    x$n_draws,
    " at ",
    nrow(x$x),
    " input point(s)\n",
    sep = ""
  )
  cat(
    "  contents: latent function draws",
    if (!is.null(x$observation) && is.null(x$observation_noise_variance)) {
      " and observation draws"
    } else if (!is.null(x$observation)) {
      " and noisy observation draws"
    },
    "\n",
    sep = ""
  )
  cat(
    "  numerical jitter: ",
    .format_jitter(x$jitter, x$covariance_scale),
    "\n",
    sep = ""
  )
  invisible(x)
}
