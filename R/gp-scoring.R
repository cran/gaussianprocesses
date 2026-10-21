# Proper scoring rules -------------------------------------------------------
#
# Scores of Gaussian predictive distributions N(mu, sigma^2) for observations
# y, with z = (y - mu) / sigma:
#   log score       log phi(z) - log sigma (higher is better);
#   CRPS            sigma (z (2 Phi(z) - 1) + 2 phi(z) - 1 / sqrt(pi))
#                   (Gneiting and Raftery, 2007; lower is better);
#   interval score  (u - l) + (2 / a) (l - y) 1{y < l}
#                   + (2 / a) (y - u) 1{y > u} for the central interval
#                   [l, u] of level 1 - a (lower is better);
#   PIT             Phi(z), uniform for calibrated forecasts.
# SMSE is the mean squared error over the variance of the observations, and
# MSLL the mean of -log p(y) minus that under a reference Gaussian, usually
# the training mean and variance (Rasmussen and Williams, 2006, section 2.5).


.crps_gaussian <- function(observed, mean, sd) {
  z <- (observed - mean) / sd

  sd * (z * (2 * stats::pnorm(z) - 1) + 2 * stats::dnorm(z) - 1 / sqrt(pi))
}


.interval_score <- function(observed, lower, upper, alpha) {
  (upper - lower) +
    (2 / alpha) * pmax(lower - observed, 0) +
    (2 / alpha) * pmax(observed - upper, 0)
}


.validate_interval_levels <- function(interval_levels) {
  if (!is.numeric(interval_levels) || length(interval_levels) == 0L ||
      any(!is.finite(interval_levels)) || any(interval_levels <= 0) ||
      any(interval_levels >= 1) || anyDuplicated(interval_levels)) {
    .gp_abort(
      paste0(
        "'interval_levels' must be distinct probabilities strictly between ",
        "0 and 1."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  sort(interval_levels)
}


#' Proper scoring rules for Gaussian predictions
#'
#' Scores Gaussian predictive distributions against observations with proper
#' scoring rules, interval scores and coverage across a grid of levels, and
#' the probability integral transform. `gp_holdout_scores()` scores a fitted
#' model on held-out data and `gp_loo_scores()` on its exact leave-one-out
#' predictions.
#'
#' @param observed Numeric vector of observations.
#' @param mean,sd Predictive means and standard deviations, one per
#'   observation; the standard deviations must be positive.
#' @param interval_levels Levels of the central prediction intervals whose
#'   coverage and interval score are reported: the calibration curve.
#' @param reference_mean,reference_variance Mean and variance of the
#'   reference Gaussian for the MSLL, usually those of the training
#'   responses. Without them the MSLL is `NA`.
#'
#' @details
#' With \eqn{z_i = (y_i - \mu_i) / \sigma_i}, the scores are
#' \describe{
#'   \item{log score}{\eqn{\log \varphi(z_i) - \log \sigma_i}, the log
#'   predictive density. Higher is better.}
#'   \item{CRPS}{\eqn{\sigma_i [z_i (2 \Phi(z_i) - 1) + 2 \varphi(z_i) -
#'   1 / \sqrt{\pi}]}, the continuous ranked probability score of
#'   Gneiting and Raftery (2007), in closed form. Lower is better; it is in
#'   the units of \eqn{y}.}
#'   \item{interval score}{\eqn{(u - l) + \frac{2}{a}(l - y)
#'   \mathbf{1}\{y < l\} + \frac{2}{a}(y - u) \mathbf{1}\{y > u\}} for the
#'   central interval \eqn{[l, u]} of level \eqn{1 - a}. Lower is better.}
#'   \item{PIT}{\eqn{\Phi(z_i)}, uniform on \eqn{[0, 1]} when the forecasts
#'   are calibrated.}
#'   \item{SMSE}{The mean squared error divided by the variance of the
#'   observations (Rasmussen and Williams, 2006, section 2.5).}
#'   \item{MSLL}{The mean of \eqn{-\log p(y_i)} minus the same under the
#'   reference Gaussian; negative values beat that reference.}
#' }
#' The log score, CRPS, and interval score are proper: in expectation they
#' are optimized by the true predictive distribution.
#'
#' `gp_holdout_scores()` predicts with [stats::predict()], so it accepts
#' every model class, and scores the observation predictive distribution,
#' including noise. Its MSLL reference is the mean and variance of the
#' training responses. `gp_loo_scores()` scores the exact leave-one-out
#' predictions of [loo_gp()] against the same reference.
#'
#' For latent models ([fit_latent_gp()]), `gp_holdout_scores()` scores the
#' predictive distribution of new responses: with a Gaussian likelihood as
#' above, with a Bernoulli likelihood by [gp_classification_scores()], whose
#' log loss uses the predictive log probabilities of
#' [likelihood_predictive()], and with a Poisson likelihood by
#' [gp_count_scores()], with the exposure of the predictions.
#'
#' @return An object of class `gaussianprocesses_scores` with
#' \describe{
#'   \item{`summary`}{Named numeric vector: `n`, `log_score` (mean log
#'   predictive density), `crps` (mean), `smse`, `msll`, `pit_mean`, and
#'   `pit_variance` (1/12 for uniform PIT values).}
#'   \item{`calibration`}{Data frame with one row per interval level: the
#'   `level`, empirical `coverage`, mean `width`, and mean
#'   `interval_score`.}
#'   \item{`pointwise`}{Data frame with one row per observation:
#'   `observed`, `mean`, `sd`, `standardized_error`, `pit`, `log_score`, and
#'   `crps`.}
#' }
#' @references Gneiting, T. and Raftery, A. E. (2007). Strictly proper
#'   scoring rules, prediction, and estimation. *Journal of the American
#'   Statistical Association*, 102(477), 359--378.
#'
#'   Rasmussen, C. E. and Williams, C. K. I. (2006). *Gaussian Processes for
#'   Machine Learning*. MIT Press. Section 2.5.
#' @seealso [gp_calibration_diagnostics()], [forecast_metrics()].
#' @examples
#' x <- seq(0, 6, length.out = 40)
#' y <- sin(x) + rnorm(40, sd = 0.1)
#' training <- seq(1, 40, by = 2)
#'
#' model <- fit_gp(x[training], y[training], rbf_kernel(),
#'   noise_variance = 0.01)
#' scores <- gp_holdout_scores(model, x[-training], y[-training])
#' scores
#' scores$calibration
#'
#' gp_loo_scores(model)$summary
#' @template stable
#' @export
gp_scores <- function(
    observed,
    mean,
    sd,
    interval_levels = c(0.5, 0.8, 0.9, 0.95),
    reference_mean = NULL,
    reference_variance = NULL) {
  .assert_numeric_vector(observed, "observed")
  .assert_numeric_vector(mean, "mean")
  .assert_numeric_vector(sd, "sd")

  if (length(mean) != length(observed) || length(sd) != length(observed)) {
    .gp_abort(
      "'observed', 'mean', and 'sd' must have the same length.",
      "gaussianprocesses_validation_error"
    )
  }

  if (any(sd <= 0)) {
    .gp_abort(
      "Predictive standard deviations must be positive.",
      "gaussianprocesses_validation_error"
    )
  }

  interval_levels <- .validate_interval_levels(interval_levels)
  standardized <- (observed - mean) / sd
  log_score <- stats::dnorm(observed, mean, sd, log = TRUE)

  msll <- NA_real_

  if (!is.null(reference_mean) || !is.null(reference_variance)) {
    .assert_numeric_scalar(reference_mean, "reference_mean")
    .assert_numeric_scalar(
      reference_variance,
      "reference_variance",
      lower = 0,
      lower_open = TRUE
    )
    msll <- mean(
      -log_score +
        stats::dnorm(observed, reference_mean, sqrt(reference_variance),
                     log = TRUE)
    )
  }

  calibration <- do.call(
    rbind,
    lapply(
      interval_levels,
      function(level) {
        quantile <- stats::qnorm((1 + level) / 2)
        lower <- mean - quantile * sd
        upper <- mean + quantile * sd

        data.frame(
          level = level,
          coverage = mean(observed >= lower & observed <= upper),
          width = mean(upper - lower),
          interval_score = mean(
            .interval_score(observed, lower, upper, 1 - level)
          )
        )
      }
    )
  )

  pit <- stats::pnorm(standardized)

  structure(
    list(
      summary = c(
        n = length(observed),
        log_score = mean(log_score),
        crps = mean(.crps_gaussian(observed, mean, sd)),
        smse = if (length(observed) > 1L) {
          mean((observed - mean)^2) / stats::var(observed)
        } else {
          NA_real_
        },
        msll = msll,
        pit_mean = mean(pit),
        pit_variance = if (length(pit) > 1L) stats::var(pit) else NA_real_
      ),
      calibration = calibration,
      pointwise = data.frame(
        observed = observed,
        mean = mean,
        sd = sd,
        standardized_error = standardized,
        pit = pit,
        log_score = log_score,
        crps = .crps_gaussian(observed, mean, sd)
      )
    ),
    class = "gaussianprocesses_scores"
  )
}


# Training responses that observe function values, for the MSLL reference.
.reference_responses <- function(model) {
  if (is.null(model$derivative)) {
    return(model$y)
  }

  model$y[model$derivative == 0L]
}


#' @rdname gp_scores
#' @param model A fitted model of any class.
#' @param x Held-out inputs, as accepted by the model's [stats::predict()]
#'   method: time values for time-series models.
#' @param y Held-out observations, one per row of `x`.
#' @param n_bins Number of probability bins for the reliability of
#'   classification scores ([gp_classification_scores()]).
#' @param seed Optional seed for the randomized PIT values of count scores
#'   ([gp_count_scores()]).
#' @param ... Further arguments to [stats::predict()], such as
#'   `observation_noise_variance` for models fitted with
#'   observation-specific noise.
#' @export
gp_holdout_scores <- function(
    model,
    x,
    y,
    interval_levels = c(0.5, 0.8, 0.9, 0.95),
    n_bins = 10L,
    seed = NULL,
    ...) {
  if (inherits(model, "gaussianprocesses_latent_model")) {
    return(.latent_holdout_scores(model, x, y, interval_levels, n_bins, seed,
                                  ...))
  }

  prediction <- stats::predict(model, x, ...)
  reference <- .reference_responses(model)

  gp_scores(
    y,
    prediction$mean,
    prediction$observation_sd,
    interval_levels = interval_levels,
    reference_mean = mean(reference),
    reference_variance = stats::var(reference)
  )
}


# Held-out scores of the predictive distribution of a latent model.
.latent_holdout_scores <- function(model, x, y, interval_levels, n_bins, seed,
                                   ...) {
  model <- .prepare_latent_model(model)
  prediction <- predict_latent_gp(model, x, ...)
  family <- model$likelihood$family

  if (identical(family, "bernoulli")) {
    observed <- .binary_response(y, "y")
    .assert_same_length(observed, prediction$mean)

    return(gp_classification_scores(
      observed,
      prediction$probability,
      n_bins = n_bins,
      log_density = likelihood_predictive(
        model$likelihood,
        prediction$mean,
        prediction$latent_variance,
        observed = observed
      )$log_density
    ))
  }

  if (identical(family, "poisson")) {
    observed <- .count_response(y, "y")
    .assert_same_length(observed, prediction$mean)

    return(gp_count_scores(
      observed,
      prediction$mean,
      prediction$latent_variance,
      exposure = prediction$exposure,
      interval_levels = interval_levels,
      seed = seed
    ))
  }

  gp_scores(
    y,
    prediction$response_mean,
    sqrt(prediction$response_variance),
    interval_levels = interval_levels,
    reference_mean = mean(model$y),
    reference_variance = stats::var(model$y)
  )
}


.assert_same_length <- function(observed, predicted) {
  if (length(observed) != length(predicted)) {
    .gp_abort(
      "'y' must have one value per row of 'x'.",
      "gaussianprocesses_validation_error"
    )
  }

  invisible(observed)
}


#' @rdname gp_scores
#' @export
gp_loo_scores <- function(model, interval_levels = c(0.5, 0.8, 0.9, 0.95)) {
  model <- .prepare_gp_model(model)
  loo <- loo_gp(model)
  reference <- .reference_responses(model)

  gp_scores(
    model$y,
    loo$mean,
    loo$sd,
    interval_levels = interval_levels,
    reference_mean = mean(reference),
    reference_variance = stats::var(reference)
  )
}


#' @export
print.gaussianprocesses_scores <- function(x, ...) {
  summary <- x$summary
  number <- function(value) format(signif(value, 4), trim = TRUE)

  cat("Predictive scores for ", summary[["n"]], " observations\n", sep = "")
  cat("  log score (mean): ", number(summary[["log_score"]]), "\n", sep = "")
  cat("  CRPS (mean): ", number(summary[["crps"]]), "\n", sep = "")
  cat("  SMSE: ", number(summary[["smse"]]), "\n", sep = "")
  cat("  MSLL: ", number(summary[["msll"]]), "\n", sep = "")
  cat(
    "  PIT mean and variance: ", number(summary[["pit_mean"]]), ", ",
    number(summary[["pit_variance"]]), " (0.5 and 0.0833 if calibrated)\n",
    sep = ""
  )
  cat(
    "  coverage: ",
    paste0(
      format(100 * x$calibration$level, trim = TRUE), "% ",
      format(round(100 * x$calibration$coverage, 1), trim = TRUE), "%",
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
  invisible(x)
}
