# Exact Gaussian-process regression ------------------------------------------

.cholesky_solve <- function(factor, rhs) {
  .assert_square_matrix(factor, name = "factor")
  .assert_numeric_vector(rhs, name = "rhs")

  if (nrow(factor) != length(rhs)) {
    .gp_abort(
      "'factor' and 'rhs' have incompatible dimensions.",
      "gaussianprocesses_validation_error"
    )
  }

  intermediate <- forwardsolve(
    t(factor),
    rhs,
    upper.tri = FALSE,
    transpose = FALSE
  )

  backsolve(
    factor,
    intermediate,
    upper.tri = TRUE,
    transpose = FALSE
  )
}


.validate_gp_training_data <- function(x, y) {
  x_matrix <- .as_input_matrix(x, name = "x")
  .assert_numeric_vector(y, name = "y")

  if (nrow(x_matrix) != length(y)) {
    .gp_abort(
      sprintf(
        paste0(
          "'x' and 'y' must contain the same number of observations. ",
          "Received %d rows in 'x' and %d values in 'y'."
        ),
        nrow(x_matrix),
        length(y)
      ),
      "gaussianprocesses_validation_error"
    )
  }

  list(x = x_matrix, y = as.numeric(y))
}


.resolve_training_noise_variance <- function(
    noise_variance,
    n_observations) {
  if (!is.numeric(noise_variance) ||
      !is.null(dim(noise_variance)) ||
      length(noise_variance) == 0L ||
      any(!is.finite(noise_variance)) ||
      any(noise_variance < 0)) {
    .gp_abort(
      paste0(
        "'noise_variance' must contain finite non-negative numeric values."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (!(length(noise_variance) %in% c(1L, n_observations))) {
    .gp_abort(
      sprintf(
        paste0(
          "'noise_variance' must contain either one value or one value per ",
          "observation. Received %d value(s) for %d observation(s)."
        ),
        length(noise_variance),
        n_observations
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (length(noise_variance) == 1L) {
    return(rep(as.numeric(noise_variance), n_observations))
  }

  as.numeric(noise_variance)
}


#' Fit an exact Gaussian-process regression model
#'
#' Fits a Gaussian-process regression model by conditioning a Gaussian prior on
#' observed responses. The implementation uses a Cholesky factorization and
#' triangular solves; it never forms an explicit matrix inverse.
#'
#' @param x Numeric vector or matrix of training inputs. Matrix rows are
#'   observations and columns are input dimensions.
#' @param y Numeric response vector with one value per training observation.
#' @param kernel A Gaussian-process kernel specification.
#' @param noise_variance Non-negative scalar observation-noise variance or a
#'   non-negative vector with one known variance per training observation.
#' @param mean Gaussian-process mean specification: [zero_mean()],
#'   [constant_mean()], or a parametric mean such as [linear_mean()] whose
#'   coefficients are fixed, estimated, or marginalized
#'   ([mean_coefficients]).
#' @param initial_jitter Non-negative jitter tried first, relative to the
#'   covariance scale (see Numerical jitter).
#' @param fallback_jitter Positive relative jitter tried after the first
#'   failed factorization.
#' @param jitter_multiplier Multiplicative increase applied to jitter after each
#'   failed factorization.
#' @param max_attempts Maximum number of Cholesky attempts.
#' @param symmetry_tolerance Relative tolerance for checking that the
#'   covariance matrix is symmetric.
#' @param derivative Optional observation types, one per observation: 0
#'   if `y[i]` observes the function value \eqn{f(x_i)}, and \eqn{d} if it
#'   observes the derivative \eqn{\partial f(x_i) / \partial x_d}. `NULL`
#'   (the default) means every observation is a function value.
#'
#' @section Derivative observations:
#' Differentiation is linear, so derivatives of a Gaussian process are
#' jointly Gaussian with its values, and observing them is exact
#' conditioning. With values at \eqn{X} and derivatives at \eqn{X'},
#' \deqn{\begin{pmatrix} y \\ y' \end{pmatrix} \sim N\left(
#' \begin{pmatrix} m(X) \\ \partial_d m(X') \end{pmatrix},
#' \begin{pmatrix} K(X, X) & \partial_{x'_d} K(X, X') \\
#' \partial_{x_d} K(X', X) & \partial_{x_d} \partial_{x'_e} K(X', X')
#' \end{pmatrix} + \Sigma \right),}
#' where each derivative row uses its own direction \eqn{d}. The blocks come
#' from the kernel's input derivatives ([kernel_input_gradient()]), so the
#' kernel must be mean-square differentiable; otherwise an error of class
#' `gaussianprocesses_smoothness_error` is raised. Parametric means
#' contribute their input derivatives. `noise_variance` may differ between
#' value and derivative observations.
#'
#' Conditioning, prediction, [log_marginal_likelihood()], its gradient, and
#' [optimize_gp()] work unchanged with this covariance;
#' [predict_gp()] and [predict_gradient_gp()] predict values and
#' derivatives at new inputs. Derivative blocks scale with
#' \eqn{1 / \ell^2}, so [gp_numerical_diagnostics()] reports the
#' conditioning of each observation type separately.
#'
#' @section Numerical jitter:
#' If \eqn{C = K + \Sigma_\varepsilon} does not factorize, the factorization
#' is retried with \eqn{C + \delta I}, where \eqn{\delta = \epsilon s} and
#' \eqn{s} is the covariance scale, the mean of the diagonal of \eqn{C}. The
#' relative jitter \eqn{\epsilon} starts at `initial_jitter` (if positive),
#' then takes the values `fallback_jitter`, `fallback_jitter *
#' jitter_multiplier`, and so on, for at most `max_attempts` attempts.
#' Because \eqn{\delta} is proportional to \eqn{C}, results do not depend
#' on the units of `y`: multiplying `y` by \eqn{c} and every variance
#' parameter by \eqn{c^2} multiplies the jitter by \eqn{c^2}. The fitted
#' model reports the absolute jitter \eqn{\delta} in `jitter`, and
#' [gp_numerical_diagnostics()] also reports \eqn{\epsilon} and \eqn{s}.
#'
#' @details
#' A scalar `noise_variance` gives the usual homoscedastic model. A vector
#' adds a known diagonal observation-noise variance for each training point.
#' Observation noise and numerical jitter are always stored separately.
#'
#' @return An object of class `gaussianprocesses_model`.
#' @seealso [gaussianprocesses_model] for the structure, versioning, and
#'   persistence of fitted models.
#' @examples
#' x <- seq(-2, 2, length.out = 15)
#' y <- sin(2 * x)
#'
#' model <- fit_gp(x, y, kernel = rbf_kernel(length_scale = 0.6), noise_variance = 0.01)
#' model
#'
#' # Known, observation-specific noise variances.
#' fit_gp(
#'   x,
#'   y,
#'   kernel = rbf_kernel(length_scale = 0.6),
#'   noise_variance = seq(0.01, 0.1, length.out = 15)
#' )
#' @templateVar experimental_parts The `derivative` argument (derivative observations) is experimental.
#' @template stable
#' @export
fit_gp <- function(
    x,
    y,
    kernel,
    noise_variance = 1e-6,
    mean = zero_mean(),
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps),
    derivative = NULL) {
  data <- .validate_gp_training_data(x, y)
  .assert_kernel(kernel)
  .assert_kernel_dimensions(kernel, ncol(data$x))
  .assert_mean_function(mean)
  mean <- .resolve_mean(mean, data$x)
  derivative <- .validate_derivative_observations(derivative, data$x)

  if (!is.null(derivative)) {
    .assert_differentiable_kernel(kernel)
  }

  training_noise_variance <- .resolve_training_noise_variance(
    noise_variance,
    n_observations = nrow(data$x)
  )

  # Only the factor of C = K + Sigma is kept; K and C are recomputed on
  # demand (see ?gaussianprocesses_model).
  observation_covariance <- .observation_covariance(
    if (is.null(derivative)) {
      evaluate_kernel(kernel, data$x)
    } else {
      .typed_covariance(kernel, data$x, derivative, data$x, derivative)
    },
    training_noise_variance
  )

  decomposition <- .stable_cholesky(
    observation_covariance,
    initial_jitter = initial_jitter,
    fallback_jitter = fallback_jitter,
    jitter_multiplier = jitter_multiplier,
    max_attempts = max_attempts,
    tolerance = symmetry_tolerance
  )

  # Estimated or marginalized coefficients need the factor of C
  # (R/mean-coefficients.R).
  fitted_mean <- .fit_mean(
    mean,
    data$x,
    data$y,
    decomposition$factor,
    derivative
  )
  prior_mean <- fitted_mean$fitted_mean
  centered_y <- data$y - prior_mean
  alpha <- .cholesky_solve(decomposition$factor, centered_y)

  structure(
    list(
      x = data$x,
      y = data$y,
      derivative = derivative,
      centered_y = centered_y,
      kernel = kernel,
      mean = mean,
      mean_fit = fitted_mean$mean_fit,
      prior_mean = prior_mean,
      noise_variance = as.numeric(noise_variance),
      training_noise_variance = training_noise_variance,
      noise_structure = if (length(noise_variance) == 1L) {
        "homoscedastic"
      } else {
        "known_heteroscedastic"
      },
      cholesky = decomposition$factor,
      alpha = alpha,
      jitter = decomposition$jitter,
      cholesky_attempts = decomposition$attempts,
      numerical_policy = list(
        initial_jitter = initial_jitter,
        fallback_jitter = fallback_jitter,
        jitter_multiplier = jitter_multiplier,
        max_attempts = max_attempts,
        symmetry_tolerance = symmetry_tolerance
      ),
      n_observations = nrow(data$x),
      n_features = ncol(data$x),
      schema_version = .model_schema_version,
      package_version = .package_version()
    ),
    class = "gaussianprocesses_model"
  )
}


.format_noise_variance <- function(model) {
  if (length(model$noise_variance) == 1L) {
    return(
      format(signif(model$noise_variance, 6), trim = TRUE)
    )
  }

  paste0(
    "[",
    format(
      signif(min(model$training_noise_variance), 6),
      trim = TRUE
    ),
    ", ",
    format(
      signif(max(model$training_noise_variance), 6),
      trim = TRUE
    ),
    "]"
  )
}


#' @export
print.gaussianprocesses_model <- function(x, ...) {
  x <- .prepare_gp_model(x, "x")
  cat("Exact Gaussian-process regression model\n")
  cat("  observations: ", x$n_observations, sep = "")

  if (!is.null(x$derivative)) {
    n_derivatives <- sum(x$derivative > 0L)
    cat(
      " (",
      x$n_observations - n_derivatives,
      " values, ",
      n_derivatives,
      " derivatives)",
      sep = ""
    )
  }

  cat("\n")
  cat("  input dimensions: ", x$n_features, "\n", sep = "")
  cat("  mean: ", .format_mean(x$mean), "\n", sep = "")

  if (!is.null(x$mean_fit)) {
    coefficients <- x$mean_fit$coefficients
    cat(
      "  mean coefficients: ",
      paste0(
        names(coefficients),
        "=",
        vapply(signif(coefficients, 6), format, character(1), trim = TRUE),
        collapse = ", "
      ),
      "\n",
      sep = ""
    )
  }

  cat(
    "  noise structure: ",
    x$noise_structure,
    "\n",
    sep = ""
  )
  cat(
    if (length(x$noise_variance) == 1L) {
      "  noise variance: "
    } else {
      "  training noise variance range: "
    },
    .format_noise_variance(x),
    "\n",
    sep = ""
  )
  cat(
    "  numerical jitter: ",
    .format_jitter(x$jitter, .model_covariance_scale(x)),
    "\n",
    sep = ""
  )

  if (!is.null(x$optimization)) {
    cat(
      "  optimization: ",
      if (isTRUE(x$optimization$converged)) "converged" else "not converged",
      " (",
      x$optimization$method,
      if (!is.null(x$optimization$gradient)) {
        paste0(", ", x$optimization$gradient, " gradient")
      },
      ")\n",
      sep = ""
    )

    if (!is.null(x$optimization$function_evaluations)) {
      cat(
        "  evaluations: ",
        x$optimization$function_evaluations,
        " objective, ",
        x$optimization$gradient_evaluations,
        " gradient, ",
        x$optimization$model_fits,
        " fits\n",
        sep = ""
      )
    }

    cat(
      "  log marginal likelihood: ",
      format(
        signif(x$optimization$log_marginal_likelihood, 7),
        trim = TRUE
      ),
      "\n",
      sep = ""
    )
  }

  cat("  kernel:\n")
  cat(paste0("    ", .format_kernel(x$kernel)), sep = "\n")
  cat("\n")
  invisible(x)
}
