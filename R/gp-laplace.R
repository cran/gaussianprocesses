# Laplace approximation for latent Gaussian-process models -------------------
#
# A latent model has f ~ N(m, K) and responses y_i ~ p(y_i | f_i). Newton's
# method finds the mode f_hat of
#   Psi(f) = log p(y | f) - (f - m)' K^{-1} (f - m) / 2
# in the parameterization f = m + K a (Rasmussen and Williams, 2006,
# Algorithm 3.1). With W = -diag(d^2 log p / df^2), each step factorizes
#   B = I + W^(1/2) K W^(1/2),
# whose eigenvalues are at least 1, with the package's Cholesky policy, and
# moves towards
#   a_new = (I + W K)^{-1} c,   c = W (f - m) + d log p / df,
# which is the Newton step f_new - m = (K^{-1} + W)^{-1} c. A step-halving
# line search on Psi guards against overshooting. The Gaussian approximation
# q(f | y) = N(f_hat, (K^{-1} + W)^{-1}) gives
#   log q(y) = log p(y | f_hat) - a' (f_hat - m) / 2 - sum(log diag(R)),
# with R'R = B, so that log |B| = 2 sum(log diag(R)).
#
# Mean coefficients with a Gaussian prior N(b, B_beta) are folded into the
# prior: m = H b and K + H B_beta H'. Estimated coefficients maximize log q
# (type-II maximum likelihood, like the hyperparameters); vague priors are not
# supported.


.prepare_latent_model <- function(model, name = "model") {
  .prepare_model(model, "gaussianprocesses_latent_model", name)
}


.assert_latent_mean <- function(mean) {
  treatment <- .mean_treatment(mean)

  if (treatment$treatment %in% c("none", "fixed") ||
      (identical(treatment$treatment, "estimated") &&
       identical(treatment$method, "ml")) ||
      (identical(treatment$treatment, "marginalized") &&
       !is.null(treatment$covariance))) {
    return(invisible(mean))
  }

  .gp_abort(
    paste0(
      "Latent models support means with fixed coefficients, coefficients ",
      "estimated by maximum likelihood (estimate_coefficients(\"ml\")), or ",
      "a Gaussian coefficient_prior() with a covariance. Restricted ",
      "likelihoods and vague priors need an exact model (fit_gp())."
    ),
    "gaussianprocesses_validation_error"
  )
}


.estimated_mean <- function(mean) {
  identical(.mean_treatment(mean)$treatment, "estimated")
}


# The mean with its coefficients fixed at `values`.
.with_coefficient_values <- function(mean, values) {
  mean$coefficients <- list(treatment = "fixed", values = unname(values))
  mean
}


# Prior mean of a latent model at x, with estimated coefficients from
# `mean_fit`.
.latent_mean_values <- function(model, x) {
  if (is.null(model$mean_fit)) {
    return(evaluate_mean(model$mean, x))
  }

  drop(.mean_basis(model$mean, x) %*% model$mean_fit$coefficients)
}


# Starting values for estimated coefficients: the generalized linear model
# without the latent Gaussian process, or zeros where it fails.
.latent_coefficient_start <- function(basis, y, likelihood, exposure) {
  family <- switch(
    paste(likelihood$family, likelihood$link),
    "gaussian identity" = stats::gaussian(),
    "bernoulli logit" = stats::binomial(),
    "bernoulli probit" = stats::binomial("probit"),
    "poisson log" = stats::poisson()
  )
  fit <- tryCatch(
    suppressWarnings(stats::glm.fit(
      basis,
      y,
      family = family,
      offset = if (identical(likelihood$family, "poisson")) log(exposure)
    )),
    error = function(error) NULL
  )
  start <- if (!is.null(fit) && isTRUE(fit$converged) &&
               all(is.finite(fit$coefficients))) {
    fit$coefficients
  } else {
    rep(0, ncol(basis))
  }

  stats::setNames(as.numeric(start), colnames(basis))
}


# Coefficients that maximize log q at fixed hyperparameters, by BFGS with
# the analytical gradient; each fit warm-starts from the previous mode.
.estimate_latent_coefficients <- function(
    x,
    y,
    response_levels,
    kernel,
    likelihood,
    mean,
    exposure,
    policy,
    max_iterations,
    convergence_tolerance) {
  basis <- .mean_basis(mean, x)
  covariance <- .latent_prior_covariance(kernel, mean, x)
  last_coefficients <- NULL
  last_model <- NULL
  warm_alpha <- NULL

  fit_at <- function(coefficients) {
    if (!identical(coefficients, last_coefficients)) {
      last_model <<- tryCatch(
        .fit_latent_model(
          x, y, response_levels, kernel, likelihood, mean, exposure, policy,
          max_iterations, convergence_tolerance,
          initial_alpha = warm_alpha,
          coefficients = coefficients
        ),
        error = function(error) NULL
      )
      last_coefficients <<- coefficients

      if (!is.null(last_model)) {
        warm_alpha <<- last_model$alpha
      }
    }

    last_model
  }

  result <- stats::optim(
    .latent_coefficient_start(basis, y, likelihood, exposure),
    function(coefficients) {
      model <- fit_at(coefficients)

      if (is.null(model) || !is.finite(model$log_marginal_likelihood)) {
        return(1e100)
      }

      -model$log_marginal_likelihood
    },
    function(coefficients) {
      model <- fit_at(coefficients)

      if (is.null(model)) {
        return(rep(0, length(coefficients)))
      }

      -unname(.laplace_gradient(
        model,
        include_kernel = FALSE,
        include_likelihood = FALSE,
        include_mean = TRUE,
        covariance = covariance
      ))
    },
    method = "BFGS",
    control = list(reltol = 1e-14, maxit = 500)
  )

  if (result$value >= 1e99) {
    .gp_abort(
      "The mean coefficients could not be estimated: no fit was finite.",
      "gaussianprocesses_numerical_error"
    )
  }

  list(
    coefficients = stats::setNames(result$par, colnames(basis)),
    optimization = list(
      converged = identical(result$convergence, 0L),
      iterations = unname(result$counts[["gradient"]]),
      message = if (is.null(result$message)) "" else result$message
    )
  )
}


.assert_log_concave <- function(entry, likelihood) {
  if (!isTRUE(entry$log_concave)) {
    .gp_abort(
      sprintf(
        paste0(
          "The Laplace approximation needs a log-concave likelihood; ",
          "%s(link=%s) is not."
        ),
        entry$label,
        likelihood$link
      ),
      c(
        "gaussianprocesses_log_concavity_error",
        "gaussianprocesses_validation_error"
      )
    )
  }

  invisible(entry)
}


# h(x)' B_beta h(y) for mean coefficients with a Gaussian prior; NULL for
# other means.
.coefficient_prior_term <- function(mean, x, y = x) {
  treatment <- .mean_treatment(mean)

  if (!identical(treatment$treatment, "marginalized")) {
    return(NULL)
  }

  basis <- .mean_basis(mean, x)
  factor <- .coefficient_prior_factor(treatment, ncol(basis))

  tcrossprod(basis %*% t(factor), .mean_basis(mean, y) %*% t(factor))
}


# Prior covariance of the latent function between x and y.
.latent_prior_covariance <- function(kernel, mean, x, y = NULL) {
  covariance <- evaluate_kernel(kernel, x, y)
  coefficient <- .coefficient_prior_term(mean, x, if (is.null(y)) x else y)

  if (!is.null(coefficient)) {
    covariance <- covariance + coefficient
  }

  covariance
}


.latent_prior_variance <- function(kernel, mean, x) {
  variance <- kernel_diagonal(kernel, x)
  treatment <- .mean_treatment(mean)

  if (identical(treatment$treatment, "marginalized")) {
    basis <- .mean_basis(mean, x)
    factor <- .coefficient_prior_factor(treatment, ncol(basis))
    variance <- variance + rowSums((basis %*% t(factor))^2)
  }

  variance
}


# Mode of Psi by Newton's method. `initial_alpha`, from an earlier fit, is a
# warm start used when it is better than starting from the prior mean.
# Iterations stop when the Newton decrement
#   lambda^2 = g' (K^{-1} + W)^{-1} g,   g = d Psi / df = d log p / df - a,
# twice the increase of Psi that the next Newton step predicts, is at most
# `convergence_tolerance`; lambda is about the distance to the mode in
# posterior standard deviations.
.laplace_mode <- function(
    covariance,
    prior_mean,
    y,
    entry,
    parameters,
    exposure,
    policy,
    max_iterations,
    convergence_tolerance,
    initial_alpha = NULL) {
  n_observations <- length(y)
  identity <- diag(n_observations)

  objective <- function(alpha, latent) {
    sum(entry$log_density(y, latent, parameters, exposure)) -
      0.5 * sum(alpha * (latent - prior_mean))
  }

  # a_new = (I + W K)^{-1} c = W^(1/2) B^{-1} W^(-1/2) c. The usual form
  # c - W^(1/2) B^{-1} W^(1/2) K c subtracts two terms of the size of
  # c = W (f - m) + d log p / df, which for large W (counts in the
  # thousands) is far larger than a, and loses all accuracy; this form has
  # only the error of solving with B, relative to a. W is floored at the
  # smallest positive double, so that W^(-1/2) c is defined where W
  # underflows; that changes B by less than rounding.
  newton <- function(alpha, latent) {
    derivatives <- entry$derivatives(y, latent, parameters, exposure)
    curvature <- pmax(-derivatives$second, .Machine$double.xmin)
    root <- sqrt(curvature)
    decomposition <- .stable_cholesky(
      outer(root, root) * covariance + identity,
      initial_jitter = policy$initial_jitter,
      fallback_jitter = policy$fallback_jitter,
      jitter_multiplier = policy$jitter_multiplier,
      max_attempts = policy$max_attempts,
      tolerance = policy$symmetry_tolerance
    )
    target <- curvature * (latent - prior_mean) + derivatives$first

    list(
      derivatives = derivatives,
      curvature = curvature,
      decomposition = decomposition,
      alpha = root * .cholesky_solve(decomposition$factor, target / root)
    )
  }

  alpha <- numeric(n_observations)
  latent <- prior_mean
  value <- objective(alpha, latent)

  if (!is.null(initial_alpha)) {
    warm_latent <- prior_mean + drop(covariance %*% initial_alpha)
    warm_value <- objective(initial_alpha, warm_latent)

    if (is.finite(warm_value) && (!is.finite(value) || warm_value > value)) {
      alpha <- initial_alpha
      latent <- warm_latent
      value <- warm_value
    }
  }

  if (!is.finite(value)) {
    .gp_abort(
      "The latent log density is not finite at the starting point.",
      "gaussianprocesses_numerical_error"
    )
  }

  iterations <- 0L
  step_halvings <- 0L
  objective_change <- NA_real_
  converged <- FALSE
  message <- ""

  repeat {
    step <- newton(alpha, latent)
    direction <- step$alpha - alpha
    latent_direction <- drop(covariance %*% direction)
    decrement <- sum((step$derivatives$first - alpha) * latent_direction)

    # lambda^2 >= 0 in exact arithmetic; a negative value beyond the
    # tolerance means the step was lost to rounding and is not convergence.
    if (abs(decrement) <= convergence_tolerance) {
      converged <- TRUE
      message <- "converged"
      break
    }

    if (iterations >= max_iterations) {
      message <- "iteration limit reached"
      break
    }

    # A step is rejected only if Psi falls by more than rounding: near the
    # mode the predicted increase lambda^2 / 2 can be smaller than the
    # rounding error of Psi, and comparing exactly would then reject the
    # Newton steps at random and stall the iterations.
    fraction <- 1
    accepted <- FALSE
    rounding <- 1e-12 * (1 + abs(value))

    for (halving in 0:40) {
      candidate_alpha <- alpha + fraction * direction
      candidate_latent <- latent + fraction * latent_direction
      candidate_value <- objective(candidate_alpha, candidate_latent)

      if (is.finite(candidate_value) && candidate_value >= value - rounding) {
        accepted <- TRUE
        break
      }

      fraction <- fraction / 2
    }

    if (!accepted) {
      message <- "the line search found no increase of the objective"
      break
    }

    iterations <- iterations + 1L
    step_halvings <- step_halvings + halving
    objective_change <- candidate_value - value
    alpha <- candidate_alpha
    latent <- candidate_latent
    value <- candidate_value
  }

  list(
    mode = latent,
    alpha = alpha,
    curvature = step$curvature,
    decomposition = step$decomposition,
    objective = value,
    log_marginal_likelihood = value -
      sum(log(diag(step$decomposition$factor))),
    diagnostics = list(
      converged = converged,
      iterations = iterations,
      decrement = decrement,
      gradient_norm = max(abs(step$derivatives$first - alpha)),
      objective_change = objective_change,
      step_halvings = step_halvings,
      message = message
    )
  )
}


# A latent model fitted at given hyperparameters; validation is left to the
# callers.
.fit_latent_model <- function(
    x,
    y,
    response_levels,
    kernel,
    likelihood,
    mean,
    exposure,
    policy,
    max_iterations,
    convergence_tolerance,
    initial_alpha = NULL,
    coefficients = NULL) {
  entry <- .likelihood_entry(likelihood)
  mean_fit <- NULL
  fixed_mean <- mean

  # Estimated coefficients: given (by the optimizer), or estimated here.
  if (.estimated_mean(mean)) {
    estimation <- NULL

    if (is.null(coefficients)) {
      estimation <- .estimate_latent_coefficients(
        x, y, response_levels, kernel, likelihood, mean, exposure, policy,
        max_iterations, convergence_tolerance
      )
      coefficients <- estimation$coefficients
    }

    names(coefficients) <- colnames(.mean_basis(mean, x[1L, , drop = FALSE]))
    fixed_mean <- .with_coefficient_values(mean, coefficients)
    mean_fit <- list(
      treatment = "estimated",
      method = "ml",
      coefficients = coefficients,
      optimization = estimation$optimization
    )
  }

  prior_mean <- evaluate_mean(fixed_mean, x)
  fit <- .laplace_mode(
    .latent_prior_covariance(kernel, fixed_mean, x),
    prior_mean,
    y,
    entry,
    likelihood$parameters,
    exposure,
    policy,
    max_iterations,
    convergence_tolerance,
    initial_alpha
  )

  structure(
    list(
      x = x,
      y = y,
      response_levels = response_levels,
      exposure = if (identical(likelihood$family, "poisson")) exposure,
      kernel = kernel,
      mean = mean,
      mean_fit = mean_fit,
      likelihood = likelihood,
      method = "laplace",
      prior_mean = prior_mean,
      mode = fit$mode,
      alpha = fit$alpha,
      curvature = fit$curvature,
      b_cholesky = fit$decomposition$factor,
      jitter = fit$decomposition$jitter,
      cholesky_attempts = fit$decomposition$attempts,
      log_marginal_likelihood = fit$log_marginal_likelihood,
      laplace = fit$diagnostics,
      numerical_policy = c(
        policy,
        list(
          max_iterations = max_iterations,
          convergence_tolerance = convergence_tolerance
        )
      ),
      n_observations = nrow(x),
      n_features = ncol(x),
      schema_version = .model_schema_version,
      package_version = .package_version()
    ),
    class = "gaussianprocesses_latent_model"
  )
}


# Validated data, likelihood, mean, and numerical settings of a latent model.
.latent_model_inputs <- function(
    x,
    y,
    kernel,
    likelihood,
    mean,
    exposure,
    method,
    max_iterations,
    convergence_tolerance,
    initial_jitter,
    fallback_jitter,
    jitter_multiplier,
    max_attempts,
    symmetry_tolerance) {
  x_matrix <- .as_input_matrix(x, name = "x")
  .assert_kernel(kernel)
  .assert_kernel_dimensions(kernel, ncol(x_matrix))
  entry <- .likelihood_entry_of(likelihood)
  .assert_log_concave(entry, likelihood)
  .assert_mean_function(mean)
  .assert_latent_mean(mean)

  if (!identical(method, "laplace")) {
    .gp_abort(
      "'method' must be \"laplace\".",
      "gaussianprocesses_validation_error"
    )
  }

  response_levels <- if (is.factor(y)) levels(y)
  y <- entry$response(y)

  if (length(y) != nrow(x_matrix)) {
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

  max_iterations <- .assert_positive_integer(max_iterations, "max_iterations")
  .assert_numeric_scalar(
    convergence_tolerance,
    "convergence_tolerance",
    lower = 0,
    lower_open = TRUE
  )

  list(
    x = x_matrix,
    y = y,
    response_levels = response_levels,
    mean = .resolve_mean(mean, x_matrix),
    exposure = .likelihood_exposure(likelihood, exposure, length(y)),
    policy = list(
      initial_jitter = initial_jitter,
      fallback_jitter = fallback_jitter,
      jitter_multiplier = jitter_multiplier,
      max_attempts = max_attempts,
      symmetry_tolerance = symmetry_tolerance
    ),
    max_iterations = max_iterations,
    convergence_tolerance = convergence_tolerance
  )
}


.warn_laplace_convergence <- function(model) {
  if (!isTRUE(model$laplace$converged)) {
    warning(
      warningCondition(
        sprintf(
          paste0(
            "The Newton iterations for the posterior mode did not converge ",
            "(%s after %d iterations; Newton decrement %s)."
          ),
          model$laplace$message,
          model$laplace$iterations,
          format(signif(model$laplace$decrement, 3))
        ),
        class = "gaussianprocesses_convergence_warning"
      )
    )
  }

  invisible(model)
}


#' Fit a latent Gaussian-process model
#'
#' Fits a Gaussian process \eqn{f \sim GP(m, k)} observed through a
#' likelihood \eqn{p(y_i \mid f(x_i))}, such as Bernoulli responses for
#' classification or Poisson counts, with the Laplace approximation to the
#' posterior of \eqn{f}. With [gaussian_likelihood()] the approximation is
#' exact and the model equals [fit_gp()].
#'
#' @param x Numeric vector or matrix of training inputs. Matrix rows are
#'   observations and columns are input dimensions.
#' @param y Responses, one per row of `x`, as the likelihood accepts them:
#'   numbers for Gaussian likelihoods; 0/1, logical values, or a two-level
#'   factor whose second level is 1 for Bernoulli likelihoods; and
#'   non-negative integer counts for Poisson likelihoods.
#' @param kernel A Gaussian-process kernel specification.
#' @param likelihood A likelihood specification ([gp_likelihoods]). It must
#'   be log-concave.
#' @param mean Gaussian-process mean specification, with fixed
#'   coefficients, coefficients estimated by maximum likelihood, or
#'   coefficients with a Gaussian [coefficient_prior()].
#' @param exposure Optional positive exposure for Poisson likelihoods: one
#'   value or one per observation (see [poisson_likelihood()]).
#' @param method The approximation: `"laplace"`.
#' @param max_iterations Maximum number of Newton iterations.
#' @param convergence_tolerance The Newton iterations stop when the Newton
#'   decrement is at most this value (see Details).
#' @param initial_jitter,fallback_jitter,jitter_multiplier,max_attempts The
#'   jitter policy for the factorization of \eqn{B}, relative to its mean
#'   diagonal, as in [fit_gp()].
#' @param symmetry_tolerance Relative tolerance for checking that \eqn{B} is
#'   symmetric.
#'
#' @details
#' The posterior of \eqn{f = f(X)} is proportional to
#' \eqn{p(y \mid f) N(f \mid m, K)}. The Laplace approximation replaces it by
#' \eqn{N(\hat f, (K^{-1} + W)^{-1})}, where \eqn{\hat f} maximizes
#' \deqn{\Psi(f) = \log p(y \mid f) - \tfrac12 (f - m)^\top K^{-1} (f - m)}
#' and \eqn{W = -\nabla\nabla \log p(y \mid \hat f)} is diagonal. For a
#' log-concave likelihood \eqn{\Psi} is strictly concave, so the mode is
#' unique; other likelihoods raise an error of class
#' `gaussianprocesses_log_concavity_error`.
#'
#' Newton's method finds the mode (Rasmussen and Williams, 2006, Algorithm
#' 3.1) in the parameterization \eqn{f = m + K a}. Each iteration factorizes
#' \deqn{B = I + W^{1/2} K W^{1/2},}
#' whose eigenvalues are at least 1, so it is well conditioned even when
#' \eqn{K} is nearly singular; the factorization follows the jitter policy
#' of [fit_gp()]. A step-halving line search on \eqn{\Psi} prevents
#' overshooting. The iterations stop when the Newton decrement
#' \deqn{\lambda^2 = g^\top (K^{-1} + W)^{-1} g, \qquad g = \nabla \Psi(f),}
#' which is twice the increase of \eqn{\Psi} that the next Newton step
#' predicts, is at most `convergence_tolerance`. \eqn{\lambda} is about the
#' distance to the mode measured in posterior standard deviations, so the
#' default of \eqn{10^{-14}} places the mode within about \eqn{10^{-7}}
#' standard deviations. `$laplace` records the number of iterations and step
#' halvings, the final decrement, the largest component of \eqn{g}, the last
#' change of \eqn{\Psi}, and whether the iterations converged; if they did
#' not, a warning of class `gaussianprocesses_convergence_warning` is
#' raised.
#'
#' The approximate log marginal likelihood is
#' \deqn{\log q(y \mid X, \theta) = \log p(y \mid \hat f) -
#' \tfrac12 (\hat f - m)^\top K^{-1} (\hat f - m) - \tfrac12 \log |B|,}
#' reported by [log_marginal_likelihood()] and [logLik()], with its gradient
#' in [log_marginal_likelihood_gradient()].
#'
#' A mean with fixed coefficients contributes \eqn{m(x)}. Coefficients with
#' a Gaussian prior \eqn{N(b, B_\beta)} are integrated out: the prior becomes
#' \eqn{m(x) = h(x)^\top b} with covariance
#' \eqn{k(x, x') + h(x)^\top B_\beta h(x')}. Coefficients estimated by
#' maximum likelihood (`estimate_coefficients("ml")`, the default of
#' [linear_mean()]) maximize \eqn{\log q(y \mid X, \theta, \beta)} at the
#' given hyperparameters, by BFGS with the analytical gradient
#' \eqn{H^\top a + [(I + K W)^{-1} H]^\top s}, where \eqn{s} is the
#' derivative of \eqn{-\frac12 \log |B|} in \eqn{\hat f}, starting from the
#' generalized linear model without the latent process. The estimate is
#' plugged in, like the hyperparameters; [coef()] returns it and `$mean_fit`
#' records it. A constant mean estimated this way is the usual mean level of
#' a count or classification model. Restricted likelihoods and vague priors
#' are not supported.
#'
#' @section Accuracy of the Laplace approximation:
#' For a few observations the exact marginal likelihood
#' \eqn{\log p(y \mid X)} can be computed by tensor-product Gauss-Hermite
#' quadrature. With an RBF kernel (variance 2, length scale 1) and the first
#' \eqn{n} of the inputs \eqn{x = 0, 0.5, 2} and of the responses below, the
#' error \eqn{\log q - \log p} of the Laplace approximation is:
#' \tabular{lllll}{
#'   Likelihood \tab Responses \tab \eqn{n = 1} \tab \eqn{n = 2} \tab
#'   \eqn{n = 3} \cr
#'   logit \tab 1, 1, 0 \tab -0.0171 \tab -0.0192 \tab -0.0387 \cr
#'   probit \tab 1, 1, 0 \tab -0.0196 \tab -0.0298 \tab -0.0450 \cr
#'   Poisson \tab 0, 2, 7 \tab 0.0125 \tab -0.0101 \tab -0.0197 \cr
#'   Poisson, mean 5 \tab 120, 150, 200 \tab -0.000681 \tab -0.00118 \tab
#'   -0.00158
#' }
#' For binary responses and small counts the error is a few hundredths and
#' grows with \eqn{n}; the approximation can over- or underestimate. For counts in the hundreds, where the likelihood is nearly
#' Gaussian in \eqn{f}, it is ten times smaller. The test suite reproduces
#' these values.
#'
#' @return An object of class `gaussianprocesses_latent_model` with fields
#'   `x`, `y` (responses in the internal coding), `response_levels` (the
#'   levels of a factor response), `exposure`, `kernel`, `mean`,
#'   `likelihood`, `method`, `mode` (\eqn{\hat f}), `alpha` (\eqn{a}),
#'   `curvature` (the diagonal of \eqn{W}), `b_cholesky` (\eqn{R} with
#'   \eqn{R^\top R = B + \delta I}), `jitter` (\eqn{\delta}),
#'   `log_marginal_likelihood`, `laplace`, `n_observations`, `n_features`,
#'   `schema_version`, and `package_version`. Other fields are internal.
#' @references Rasmussen, C. E. and Williams, C. K. I. (2006). *Gaussian
#'   Processes for Machine Learning*. MIT Press. Chapters 3 and 5.
#' @seealso [predict_latent_gp()], [optimize_latent_gp()],
#'   [sample_gp_posterior()].
#' @examples
#' x <- seq(-3, 3, length.out = 30)
#' y <- as.numeric(sin(x) + 0.3 * cos(5 * x) > 0)
#'
#' model <- fit_latent_gp(x, y, rbf_kernel(variance = 4), bernoulli_likelihood())
#' model
#' predict_latent_gp(model, c(-2, 0, 2))$probability
#' @template experimental
#' @export
fit_latent_gp <- function(
    x,
    y,
    kernel,
    likelihood,
    mean = zero_mean(),
    exposure = NULL,
    method = "laplace",
    max_iterations = 100L,
    convergence_tolerance = 1e-14,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps)) {
  inputs <- .latent_model_inputs(
    x,
    y,
    kernel,
    likelihood,
    mean,
    exposure,
    method,
    max_iterations,
    convergence_tolerance,
    initial_jitter,
    fallback_jitter,
    jitter_multiplier,
    max_attempts,
    symmetry_tolerance
  )

  model <- .fit_latent_model(
    inputs$x,
    inputs$y,
    inputs$response_levels,
    kernel,
    likelihood,
    inputs$mean,
    inputs$exposure,
    inputs$policy,
    inputs$max_iterations,
    inputs$convergence_tolerance
  )

  .warn_laplace_convergence(model)
}


# Gradient of the approximate log marginal likelihood (Rasmussen and
# Williams, 2006, Algorithm 5.1) in the optimizer coordinates of the kernel
# parameters and, with `include_likelihood`, of the likelihood parameters,
# named "likelihood.<name>". With Sigma = (K^{-1} + W)^{-1} and
# Z = W^(1/2) B^{-1} W^(1/2), the explicit part for a kernel parameter with
# dK = d K / d eta is
#   a' dK a / 2 - tr(Z dK) / 2,
# and for a likelihood parameter phi
#   sum(d log p / d phi) + sum(diag(Sigma) d^2 log p / df^2 d phi) / 2.
# The mode moves with the parameters,
#   d f_hat / d eta = (I + K W)^{-1} v = v - K Z v,
# with v = dK d log p / df for kernel parameters and
# v = K d^2 log p / df d phi for likelihood parameters, and log q depends on
# f_hat only through -log |B| / 2, whose derivative in f_hat_i is
# Sigma_ii (d^3 log p / df^3)_i / 2. For estimated mean coefficients beta,
# with m = H beta, the explicit part is H' a and d f_hat / d beta =
# (I + K W)^{-1} H; those derivatives, named "mean.<name>", are returned with
# `include_mean`. They are zero at the estimate, so log q at the estimate has
# the gradient of its other parameters, as for a profile likelihood.
# `implicit = FALSE` leaves out the implicit part; it exists for the tests.
# `covariance` may pass K (with any coefficient-prior term) when the kernel
# derivatives are not needed.
.laplace_gradient <- function(
    model,
    include_kernel = TRUE,
    include_likelihood = TRUE,
    include_mean = FALSE,
    implicit = TRUE,
    covariance = NULL) {
  entry <- .likelihood_entry(model$likelihood)
  parameters <- model$likelihood$parameters
  kernel_terms <- if (include_kernel) {
    .kernel_value_and_gradient(model$kernel, model$x, model$x)
  }

  if (is.null(covariance)) {
    covariance <- if (include_kernel) {
      kernel_terms$value
    } else {
      evaluate_kernel(model$kernel, model$x)
    }
    coefficient <- .coefficient_prior_term(model$mean, model$x)

    if (!is.null(coefficient)) {
      covariance <- covariance + coefficient
    }
  }

  exposure <- if (is.null(model$exposure)) {
    rep(1, model$n_observations)
  } else {
    model$exposure
  }
  root <- sqrt(model$curvature)
  factor <- model$b_cholesky
  inverse_b <- chol2inv(factor)
  z_matrix <- outer(root, root) * inverse_b
  whitened <- backsolve(factor, root * covariance, transpose = TRUE)
  posterior_variance <- diag(covariance) - colSums(whitened^2)
  derivatives <- entry$derivatives(
    model$y,
    model$mode,
    parameters,
    exposure
  )
  sensitivity <- 0.5 * posterior_variance * derivatives$third
  alpha <- model$alpha

  implicit_part <- function(direction) {
    if (!implicit) {
      return(0)
    }

    moved <- direction - drop(covariance %*% drop(z_matrix %*% direction))
    sum(sensitivity * moved)
  }

  gradient <- if (include_kernel) {
    vapply(
      kernel_terms$gradient,
      function(derivative) {
        0.5 * sum(alpha * drop(derivative %*% alpha)) -
          0.5 * sum(z_matrix * derivative) +
          implicit_part(drop(derivative %*% derivatives$first))
      },
      numeric(1)
    )
  } else {
    numeric(0)
  }

  if (include_likelihood && length(parameters) > 0L) {
    terms <- entry$parameter_derivatives(
      model$y,
      model$mode,
      parameters,
      exposure
    )
    likelihood_gradient <- vapply(
      terms,
      function(term) {
        sum(term$log_density) +
          0.5 * sum(posterior_variance * term$second) +
          implicit_part(drop(covariance %*% term$first))
      },
      numeric(1)
    )
    names(likelihood_gradient) <- paste0("likelihood.", names(terms))
    gradient <- c(gradient, likelihood_gradient)
  }

  if (include_mean && !is.null(model$mean_fit)) {
    basis <- .mean_basis(model$mean, model$x)
    mean_gradient <- vapply(
      seq_len(ncol(basis)),
      function(j) sum(basis[, j] * alpha) + implicit_part(basis[, j]),
      numeric(1)
    )
    names(mean_gradient) <- paste0("mean.", colnames(basis))
    gradient <- c(gradient, mean_gradient)
  }

  gradient
}


#' @template stable
#' @rdname gp_coefficients
#' @export
coef.gaussianprocesses_latent_model <- function(object, ...) {
  model <- .prepare_latent_model(object, "object")
  treatment <- .mean_treatment(model$mean)

  if (!identical(treatment$treatment, "marginalized")) {
    return(.model_coefficients(model))
  }

  # Under the Laplace approximation the posterior mean of coefficients with
  # prior N(b, B) is b + B H' a, because f - m = (K + H B H') a.
  basis <- .mean_basis(model$mean, model$x)
  factor <- .coefficient_prior_factor(treatment, ncol(basis))
  posterior_mean <- .coefficient_prior_mean(treatment, ncol(basis)) +
    drop(crossprod(factor, factor %*% crossprod(basis, model$alpha)))

  stats::setNames(posterior_mean, colnames(basis))
}


#' @export
print.gaussianprocesses_latent_model <- function(x, ...) {
  x <- .prepare_latent_model(x, "x")
  cat("Latent Gaussian-process model (Laplace approximation)\n")
  cat("  observations: ", x$n_observations, "\n", sep = "")
  cat("  input dimensions: ", x$n_features, "\n", sep = "")
  cat(
    "  likelihood: ",
    utils::capture.output(print(x$likelihood)),
    "\n",
    sep = ""
  )
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
    "  Newton iterations: ",
    x$laplace$iterations,
    " (",
    if (isTRUE(x$laplace$converged)) "converged" else x$laplace$message,
    ")\n",
    sep = ""
  )
  cat(
    "  approximate log marginal likelihood: ",
    format(signif(x$log_marginal_likelihood, 7), trim = TRUE),
    "\n",
    sep = ""
  )
  cat(
    "  numerical jitter: ",
    format(signif(x$jitter, 3), trim = TRUE),
    "\n",
    sep = ""
  )

  if (!is.null(x$optimization)) {
    cat(
      "  optimization: ",
      if (isTRUE(x$optimization$converged)) "converged" else "not converged",
      " (",
      x$optimization$method,
      ", ",
      x$optimization$gradient,
      " gradient)\n",
      sep = ""
    )
  }

  cat("  kernel:\n")
  cat(paste0("    ", .format_kernel(x$kernel)), sep = "\n")
  cat("\n")
  invisible(x)
}
