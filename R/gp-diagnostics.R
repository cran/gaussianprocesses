# Gaussian-process diagnostics -----------------------------------------------

.precision_diagonal_from_cholesky <- function(factor) {
  .assert_square_matrix(factor, name = "factor")

  # Exact LOO needs diag(C^{-1}) itself. With C = R'R this is the row-wise
  # sum of squares of R^{-1}, obtained by triangular back-substitution; the
  # dense covariance is never inverted.
  inverse_factor <- backsolve(
    factor,
    diag(nrow(factor)),
    upper.tri = TRUE,
    transpose = FALSE
  )

  rowSums(inverse_factor^2)
}


#' Exact leave-one-out Gaussian-process diagnostics
#'
#' Computes exact leave-one-out (LOO) predictive summaries from a fitted exact
#' Gaussian-process model without refitting the model once per observation.
#' The calculation uses the fitted precision diagonal implied by the stored
#' Cholesky factor.
#'
#' @param model A fitted `gaussianprocesses_model`.
#' @param interval_level Probability level for LOO prediction intervals.
#'
#' @details
#' When numerical jitter was required during fitting, the LOO calculation
#' necessarily reflects the stabilized covariance used in that factorization;
#' the amount of jitter is reported separately by [gp_numerical_diagnostics()].
#'
#' Estimated mean coefficients stay at their full-data estimate, as in
#' [predict_gp()]. Marginalized coefficients are integrated out as in the
#' posterior: the precision diagonal is that of
#' \eqn{P_A = C^{-1} - C^{-1} H (B^{-1} + H^\top C^{-1} H)^{-1} H^\top C^{-1}},
#' so with a vague prior each prediction re-estimates the coefficients
#' without the left-out observation.
#'
#' @return A data frame with observed values, LOO predictive means and
#'   variances, residuals, standardized residuals, probability integral
#'   transform (PIT) values, interval bounds, coverage indicators, and
#'   pointwise log predictive densities.
#' @examples
#' x <- seq(-2, 2, length.out = 20)
#' model <- fit_gp(
#'   x,
#'   sin(2 * x),
#'   kernel = rbf_kernel(length_scale = 0.6),
#'   noise_variance = 0.01
#' )
#'
#' loo <- loo_gp(model)
#' head(loo[, c("observed", "mean", "sd", "pit", "covered")])
#' @template stable
#' @export
loo_gp <- function(model, interval_level = 0.95) {
  model <- .prepare_gp_model(model)
  .validate_interval_level(interval_level)

  precision_diagonal <- .precision_diagonal_from_cholesky(
    model$cholesky
  )

  if (.marginalized_mean(model)) {
    precision_diagonal <- precision_diagonal -
      rowSums(.coefficient_projection(model)^2)
  }

  if (any(!is.finite(precision_diagonal)) ||
      any(precision_diagonal <= 0)) {
    .gp_abort(
      paste0(
        "The fitted precision diagonal contains invalid values. ",
        "LOO diagnostics cannot be computed reliably."
      ),
      "gaussianprocesses_numerical_error"
    )
  }

  loo_variance <- 1 / precision_diagonal
  loo_sd <- sqrt(loo_variance)
  loo_mean <- model$y - model$alpha / precision_diagonal
  residual <- model$y - loo_mean
  standardized_residual <- residual / loo_sd
  pit <- stats::pnorm(
    model$y,
    mean = loo_mean,
    sd = loo_sd
  )
  quantile <- stats::qnorm((1 + interval_level) / 2)
  lower <- loo_mean - quantile * loo_sd
  upper <- loo_mean + quantile * loo_sd
  covered <- model$y >= lower & model$y <= upper
  log_predictive_density <- stats::dnorm(
    model$y,
    mean = loo_mean,
    sd = loo_sd,
    log = TRUE
  )

  data.frame(
    observation = seq_len(model$n_observations),
    observed = model$y,
    mean = loo_mean,
    variance = loo_variance,
    sd = loo_sd,
    residual = residual,
    standardized_residual = standardized_residual,
    pit = pit,
    lower = lower,
    upper = upper,
    covered = covered,
    log_predictive_density = log_predictive_density,
    check.names = FALSE
  )
}


#' Summarize Gaussian-process uncertainty calibration
#'
#' Summarizes exact LOO predictive diagnostics. The returned statistics are
#' descriptive diagnostics rather than a pass/fail calibration test.
#'
#' @details
#' PIT values should be assessed as a distribution rather than individually.
#' Under a correctly specified continuous predictive distribution they are
#' uniform in repeated sampling.
#'
#' @inheritParams loo_gp
#'
#' @return A named list of calibration summaries.
#' @examples
#' x <- seq(-2, 2, length.out = 20)
#' model <- fit_gp(
#'   x,
#'   sin(2 * x),
#'   kernel = rbf_kernel(length_scale = 0.6),
#'   noise_variance = 0.01
#' )
#'
#' gp_calibration_diagnostics(model, interval_level = 0.9)
#' @template experimental
#' @export
gp_calibration_diagnostics <- function(
    model,
    interval_level = 0.95) {
  model <- .prepare_gp_model(model)

  .calibration_summary(
    loo_gp(model, interval_level = interval_level),
    interval_level = interval_level
  )
}


# Calibration summaries of a leave-one-out table from loo_gp(), so that
# gp_diagnostics() can reuse one leave-one-out computation.
.calibration_summary <- function(loo, interval_level) {
  list(
    interval_level = interval_level,
    empirical_coverage = mean(loo$covered),
    coverage_error = mean(loo$covered) - interval_level,
    mean_standardized_residual = mean(
      loo$standardized_residual
    ),
    sd_standardized_residual = stats::sd(
      loo$standardized_residual
    ),
    mean_absolute_standardized_residual = mean(
      abs(loo$standardized_residual)
    ),
    pit_mean = mean(loo$pit),
    pit_variance = stats::var(loo$pit),
    mean_log_predictive_density = mean(
      loo$log_predictive_density
    )
  )
}


#' Summarize numerical conditioning of a fitted GP
#'
#' Reports numerical quantities separately from statistical fit diagnostics.
#' In particular, observation-noise variance and numerical jitter remain
#' distinct.
#'
#' @param model A fitted `gaussianprocesses_model`.
#' @param condition_tolerance Positive reciprocal-condition threshold used only
#'   to flag a potentially ill-conditioned stabilized covariance matrix.
#'
#' @return A named list of numerical diagnostics, including:
#' \describe{
#'   \item{`observation_noise_variance`}{The statistical noise variance.}
#'   \item{`numerical_jitter`}{The diagonal jitter \eqn{\delta} that the
#'   factorization of the training covariance \eqn{C} needed.}
#'   \item{`covariance_scale`}{The mean of the diagonal of \eqn{C}.}
#'   \item{`relative_jitter`}{`numerical_jitter / covariance_scale`.}
#'   \item{`cholesky_attempts`}{Factorization attempts, 1 if no jitter was
#'   needed.}
#'   \item{`reciprocal_condition_number`, `condition_number`}{The 1-norm
#'   condition estimate of \eqn{C + \delta I} and its reciprocal.}
#'   \item{`potentially_ill_conditioned`}{Whether the reciprocal condition
#'   number is below `condition_tolerance`.}
#'   \item{`observation_types`}{A data frame with one row per observation
#'   type (`"value"`, or `"derivative x1"`, ... for derivative
#'   observations): the number of observations, the mean diagonal of their
#'   block of \eqn{C}, and the reciprocal condition number of that block
#'   with the jitter. Derivative blocks scale with \eqn{1 / \ell^2}, so
#'   their conditioning can differ from that of the values.}
#' }
#' @examples
#' # Duplicated inputs without observation noise need numerical jitter, which is
#' # reported separately from the (zero) noise variance.
#' model <- fit_gp(c(0, 0, 1), c(1, 1, 0), kernel = rbf_kernel(), noise_variance = 0)
#'
#' diagnostics <- gp_numerical_diagnostics(model)
#' diagnostics[c(
#'   "observation_noise_variance",
#'   "numerical_jitter",
#'   "relative_jitter",
#'   "cholesky_attempts"
#' )]
#' @template experimental
#' @export
gp_numerical_diagnostics <- function(
    model,
    condition_tolerance = sqrt(.Machine$double.eps)) {
  model <- .prepare_gp_model(model)
  .assert_numeric_scalar(
    condition_tolerance,
    "condition_tolerance",
    lower = 0,
    lower_open = TRUE
  )

  observation_covariance <- .model_observation_covariance(model)
  stabilized_covariance <- observation_covariance +
    diag(model$jitter, model$n_observations)

  reciprocal_condition <- rcond(stabilized_covariance)
  condition_number <- if (
    is.finite(reciprocal_condition) &&
      reciprocal_condition > 0
  ) {
    1 / reciprocal_condition
  } else {
    Inf
  }

  covariance_scale <- .covariance_scale(observation_covariance)

  list(
    observation_noise_variance = model$noise_variance,
    numerical_jitter = model$jitter,
    relative_jitter = model$jitter / covariance_scale,
    covariance_scale = covariance_scale,
    cholesky_attempts = model$cholesky_attempts,
    used_numerical_jitter = model$jitter > 0,
    reciprocal_condition_number = reciprocal_condition,
    condition_number = condition_number,
    condition_tolerance = condition_tolerance,
    potentially_ill_conditioned =
      !is.finite(reciprocal_condition) ||
      reciprocal_condition < condition_tolerance,
    min_cholesky_diagonal = min(diag(model$cholesky)),
    max_cholesky_diagonal = max(diag(model$cholesky)),
    observation_types = .conditioning_by_type(model, stabilized_covariance)
  )
}


.conditioning_by_type <- function(model, covariance) {
  types <- if (is.null(model$derivative)) {
    integer(model$n_observations)
  } else {
    model$derivative
  }
  present <- sort(unique(types))

  data.frame(
    type = .observation_type_labels(present),
    observations = vapply(present, function(type) sum(types == type), 1L),
    covariance_scale = vapply(
      present,
      function(type) mean(diag(covariance)[types == type]),
      numeric(1)
    ),
    reciprocal_condition_number = vapply(
      present,
      function(type) {
        rows <- types == type
        rcond(covariance[rows, rows, drop = FALSE])
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}


#' Comprehensive diagnostics for an exact Gaussian-process model
#'
#' Combines exact leave-one-out predictive diagnostics, uncertainty-calibration
#' summaries, and numerical-conditioning diagnostics. The components are kept
#' separate so predictive misfit is not automatically attributed to numerical
#' instability, and numerical instability is not mistaken for model misfit.
#'
#' @param model A fitted `gaussianprocesses_model`.
#' @param interval_level Probability level for LOO prediction intervals.
#' @param condition_tolerance Positive reciprocal-condition threshold used by
#'   the numerical diagnostic summary.
#'
#' @return An object of class `gaussianprocesses_diagnostics`.
#' @examples
#' x <- seq(-2, 2, length.out = 20)
#' model <- fit_gp(
#'   x,
#'   sin(2 * x),
#'   kernel = rbf_kernel(length_scale = 0.6),
#'   noise_variance = 0.01
#' )
#'
#' gp_diagnostics(model)
#' @templateVar experimental_parts Its `numerical` component comes from the experimental [gp_numerical_diagnostics()].
#' @template stable
#' @export
gp_diagnostics <- function(
    model,
    interval_level = 0.95,
    condition_tolerance = sqrt(.Machine$double.eps)) {
  model <- .prepare_gp_model(model)

  loo <- loo_gp(
    model,
    interval_level = interval_level
  )
  calibration <- .calibration_summary(loo, interval_level = interval_level)
  numerical <- gp_numerical_diagnostics(
    model,
    condition_tolerance = condition_tolerance
  )

  structure(
    list(
      loo = loo,
      calibration = calibration,
      numerical = numerical
    ),
    class = "gaussianprocesses_diagnostics"
  )
}


#' @export
print.gaussianprocesses_diagnostics <- function(x, ...) {
  cat("Gaussian-process diagnostics\n")
  cat(
    "  LOO empirical coverage: ",
    format(
      signif(x$calibration$empirical_coverage, 4),
      trim = TRUE
    ),
    " (nominal ",
    format(
      signif(x$calibration$interval_level, 4),
      trim = TRUE
    ),
    ")\n",
    sep = ""
  )
  cat(
    "  mean standardized residual: ",
    format(
      signif(
        x$calibration$mean_standardized_residual,
        5
      ),
      trim = TRUE
    ),
    "\n",
    sep = ""
  )
  cat(
    "  SD standardized residual: ",
    format(
      signif(
        x$calibration$sd_standardized_residual,
        5
      ),
      trim = TRUE
    ),
    "\n",
    sep = ""
  )
  cat(
    "  numerical jitter: ",
    .format_jitter(
      x$numerical$numerical_jitter,
      x$numerical$covariance_scale
    ),
    "\n",
    sep = ""
  )
  cat(
    "  reciprocal condition number: ",
    format(
      signif(
        x$numerical$reciprocal_condition_number,
        5
      ),
      trim = TRUE
    ),
    "\n",
    sep = ""
  )
  invisible(x)
}
