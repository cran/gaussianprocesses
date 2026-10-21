# Hyperparameter optimization for latent Gaussian-process models -----------

.likelihood_parameter_prefix <- "likelihood."

.mean_parameter_prefix <- "mean."


# Kernel, likelihood, and estimated mean coefficients at optimizer
# coordinates. Likelihood parameters are named "likelihood.<name>" and mean
# coefficients "mean.<name>" among the optimized parameters.
.decode_latent_parameters <- function(theta, kernel, likelihood, constraints) {
  natural <- .from_coordinates(theta, constraints)
  is_likelihood <- startsWith(names(natural), .likelihood_parameter_prefix)
  is_mean <- startsWith(names(natural), .mean_parameter_prefix)

  for (name in names(natural)[is_likelihood]) {
    parameter <- substring(name, nchar(.likelihood_parameter_prefix) + 1L)
    likelihood$parameters[[parameter]] <- unname(natural[[name]])
  }

  coefficients <- if (any(is_mean)) {
    stats::setNames(
      unname(natural[is_mean]),
      substring(names(natural)[is_mean], nchar(.mean_parameter_prefix) + 1L)
    )
  }

  list(
    kernel = update_kernel_parameters(
      kernel,
      natural[!is_likelihood & !is_mean]
    ),
    likelihood = likelihood,
    coefficients = coefficients,
    parameters = natural
  )
}


# Negative approximate log marginal likelihood and its gradient in optimizer
# coordinates, as .optimization_problem() for exact models. Each fit starts
# its Newton iterations from the mode of the previous successful fit when
# that is better than the prior mean.
.latent_optimization_problem <- function(inputs, kernel, likelihood, setup) {
  last_theta <- NULL
  last_model <- NULL
  warm_alpha <- NULL
  evaluations <- c(objective = 0L, gradient = 0L, fits = 0L)

  fit_at <- function(theta) {
    if (!identical(theta, last_theta)) {
      evaluations[["fits"]] <<- evaluations[["fits"]] + 1L
      decoded <- .decode_latent_parameters(
        theta,
        kernel,
        likelihood,
        setup$constraints
      )

      last_model <<- tryCatch(
        .fit_latent_model(
          inputs$x,
          inputs$y,
          inputs$response_levels,
          decoded$kernel,
          decoded$likelihood,
          inputs$mean,
          inputs$exposure,
          inputs$policy,
          inputs$max_iterations,
          inputs$convergence_tolerance,
          initial_alpha = warm_alpha,
          coefficients = decoded$coefficients
        ),
        error = function(error) NULL
      )
      last_theta <<- theta

      if (!is.null(last_model)) {
        warm_alpha <<- last_model$alpha
      }
    }

    last_model
  }

  objective <- function(theta) {
    evaluations[["objective"]] <<- evaluations[["objective"]] + 1L
    model <- fit_at(theta)

    if (is.null(model) || !is.finite(model$log_marginal_likelihood)) {
      return(1e100)
    }

    -model$log_marginal_likelihood
  }

  gradient <- function(theta) {
    evaluations[["gradient"]] <<- evaluations[["gradient"]] + 1L
    model <- fit_at(theta)
    failed <- rep(0, length(theta))

    if (is.null(model)) {
      return(failed)
    }

    value <- -.laplace_gradient(model, include_mean = TRUE)[names(theta)]

    if (any(!is.finite(value))) {
      return(failed)
    }

    unname(value)
  }

  list(
    objective = objective,
    gradient = gradient,
    evaluations = function() evaluations,
    reset = function() {
      evaluations[] <<- 0L
      invisible(NULL)
    }
  )
}


# Optimized parameters at a bound of the search, as "name (lower)" or
# "name (upper)". L-BFGS-B returns active bounds exactly; the tolerance only
# absorbs rounding in the coordinates.
.parameters_at_bounds <- function(theta, setup, tolerance = 1e-6) {
  lower <- setup$lower[names(theta)]
  upper <- setup$upper[names(theta)]
  at_lower <- is.finite(lower) &
    theta <= lower + tolerance * pmax(1, abs(lower))
  at_upper <- is.finite(upper) &
    theta >= upper - tolerance * pmax(1, abs(upper))

  as.character(c(
    if (any(at_lower)) paste(names(theta)[at_lower], "(lower)"),
    if (any(at_upper)) paste(names(theta)[at_upper], "(upper)")
  ))
}


.warn_parameters_at_bounds <- function(at_bound, likelihood) {
  if (length(at_bound) == 0L) {
    return(invisible(NULL))
  }

  warning(
    warningCondition(
      paste0(
        "Optimized parameter(s) at a bound of the search: ",
        paste(at_bound, collapse = ", "),
        ". The approximate marginal likelihood may increase beyond it.",
        if (identical(likelihood$family, "bernoulli")) {
          paste0(
            " For binary responses a signal variance at its upper bound ",
            "usually means that the classes are separable, so that larger ",
            "variances keep fitting them better; the bound then determines ",
            "the fit."
          )
        }
      ),
      class = "gaussianprocesses_bound_warning"
    )
  )
}


#' Optimize the hyperparameters of a latent Gaussian-process model
#'
#' Estimates kernel parameters and, optionally, likelihood parameters by
#' maximizing the Laplace approximation to the log marginal likelihood, with
#' the optimizer of [optimize_gp()]: L-BFGS-B on unconstrained coordinates,
#' deterministic multiple starts, and analytical gradients.
#'
#' @inheritParams fit_latent_gp
#' @param kernel Initial kernel specification.
#' @param likelihood Likelihood specification; its parameters, if any, are
#'   the starting values of the likelihood parameters.
#' @param optimize_likelihood If `TRUE`, optimize the likelihood parameters,
#'   such as the variance of [gaussian_likelihood()], jointly with the kernel
#'   parameters. Likelihoods without parameters ignore it.
#' @param n_starts,start_spread,lower,upper,control,gradient As in
#'   [optimize_gp()].
#'
#' @details
#' The objective is \eqn{\log q(y \mid X, \theta)} of [fit_latent_gp()]. Its
#' gradient (Rasmussen and Williams, 2006, Algorithm 5.1) has an explicit
#' part, with the mode \eqn{\hat f} held fixed, and an implicit part through
#' the dependence of \eqn{\hat f} on \eqn{\theta}:
#' \deqn{\frac{\partial \hat f}{\partial \theta_j} = (I + K W)^{-1}
#' \frac{\partial K}{\partial \theta_j} \nabla \log p(y \mid \hat f),}
#' which needs the third derivatives of the log likelihood
#' ([evaluate_likelihood()]). Likelihood parameters are named
#' `likelihood.<name>`, for example `likelihood.variance`.
#'
#' Mean coefficients estimated by maximum likelihood
#' (`estimate_coefficients("ml")`) are optimized jointly with the other
#' parameters, on the real line, as `mean.<name>` (for example
#' `mean.intercept`), starting from the generalized linear model fitted
#' without the latent process. Their gradient is
#' \eqn{H^\top \nabla \log p(y \mid \hat f)} plus the implicit part through
#' \eqn{\partial \hat f / \partial \beta = (I + K W)^{-1} H}.
#'
#' Each evaluation starts the Newton iterations from the previous mode when
#' that is better than the prior mean, which saves most of the iterations
#' once the optimizer takes small steps. The returned model is refitted at
#' the optimum from the prior mean, so it equals [fit_latent_gp()] with the
#' optimized parameters.
#'
#' `$optimization` has the fields of [optimize_gp()], with
#' `optimize_likelihood` and `approximation` in place of the noise fields,
#' and `at_bound`, the parameters that ended at `lower` or `upper` (on the
#' natural scale) or at the limits that [optimize_gp()] describes for
#' changepoint locations and spectral frequencies.
#'
#' @section Bounds and separable data:
#' Every positive parameter is searched within `[lower, upper]`. If any
#' parameter ends at a bound, a warning of class
#' `gaussianprocesses_bound_warning` names it, because the likelihood may
#' still increase beyond the bound and the bound, not the data, then
#' determines the estimate. With binary responses this happens when the
#' classes are separable: larger signal variances let the latent function
#' separate them more sharply, and the approximate marginal likelihood
#' keeps increasing up to very large variances. On twenty separable points
#' on a line it peaks only near a variance of 1800, with a nearly flat
#' maximum, and a smaller `upper` ends the search at the bound with this
#' warning. Setting `upper` to the largest plausible variance makes the
#' limit explicit.
#'
#' @return A fitted `gaussianprocesses_latent_model` with optimized
#'   parameters and optimizer diagnostics in `$optimization`.
#' @seealso [fit_latent_gp()].
#' @examples
#' x <- seq(0, 10, length.out = 40)
#' counts <- c(1, 0, 2, 1, 3, 2, 4, 5, 3, 6, 7, 5, 8, 6, 9, 7, 6, 8, 5, 6,
#'   4, 5, 3, 4, 2, 3, 2, 1, 2, 1, 1, 0, 1, 2, 1, 3, 2, 4, 3, 5)
#'
#' model <- optimize_latent_gp(x, counts, rbf_kernel(), poisson_likelihood(),
#'   n_starts = 1)
#' kernel_parameters(model$kernel, flatten = TRUE)
#' @template experimental
#' @export
optimize_latent_gp <- function(
    x,
    y,
    kernel,
    likelihood,
    mean = zero_mean(),
    exposure = NULL,
    optimize_likelihood = TRUE,
    n_starts = 3L,
    start_spread = 1,
    lower = 1e-8,
    upper = 1e8,
    control = list(maxit = 500),
    method = "laplace",
    max_iterations = 100L,
    convergence_tolerance = 1e-14,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps),
    gradient = c("analytical", "numerical")) {
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
  gradient <- match.arg(gradient)

  if (!is.logical(optimize_likelihood) ||
      length(optimize_likelihood) != 1L ||
      is.na(optimize_likelihood)) {
    .gp_abort(
      "'optimize_likelihood' must be TRUE or FALSE.",
      "gaussianprocesses_validation_error"
    )
  }

  n_starts <- .assert_positive_integer(n_starts, "n_starts")
  .assert_numeric_scalar(start_spread, "start_spread", lower = 0)
  bounds <- .validate_optimization_bounds(lower, upper)
  .assert_optimizer_control(control)

  definitions <- .likelihood_entry(likelihood)$parameters
  extra_initial <- numeric(0)
  extra_constraints <- list()

  if (optimize_likelihood) {
    for (parameter in names(definitions)) {
      name <- paste0(.likelihood_parameter_prefix, parameter)
      constraint <- definitions[[parameter]]$constraint
      value <- likelihood$parameters[[parameter]]

      if (identical(constraint$type, "positive")) {
        value <- max(value, bounds[["lower"]])
      }

      extra_initial[[name]] <- value
      extra_constraints[[name]] <- constraint
    }
  }

  if (.estimated_mean(inputs$mean)) {
    start <- .latent_coefficient_start(
      .mean_basis(inputs$mean, inputs$x),
      inputs$y,
      likelihood,
      inputs$exposure
    )

    for (coefficient in names(start)) {
      name <- paste0(.mean_parameter_prefix, coefficient)
      extra_initial[[name]] <- start[[coefficient]]
      extra_constraints[[name]] <- .real_constraint()
    }
  }

  setup <- .optimization_setup(
    kernel,
    inputs$x,
    extra_initial = extra_initial,
    extra_constraints = extra_constraints,
    bounds = bounds,
    n_starts = n_starts,
    start_spread = start_spread
  )
  problem <- .latent_optimization_problem(inputs, kernel, likelihood, setup)
  run <- .run_optimization_starts(problem, setup, control, gradient)
  decoded <- .decode_latent_parameters(
    run$best_theta,
    kernel,
    likelihood,
    setup$constraints
  )

  fitted <- .fit_latent_model(
    inputs$x,
    inputs$y,
    inputs$response_levels,
    decoded$kernel,
    decoded$likelihood,
    inputs$mean,
    inputs$exposure,
    inputs$policy,
    inputs$max_iterations,
    inputs$convergence_tolerance,
    coefficients = decoded$coefficients
  )

  fitted$optimization <- .optimization_record(
    setup,
    run,
    gradient = gradient,
    bounds = bounds,
    n_starts = n_starts,
    start_spread = start_spread,
    specific = list(
      optimize_likelihood = optimize_likelihood && length(definitions) > 0L,
      approximation = "laplace"
    ),
    results = list(
      log_marginal_likelihood = fitted$log_marginal_likelihood,
      at_bound = .parameters_at_bounds(run$best_theta, setup)
    )
  )

  .warn_parameters_at_bounds(fitted$optimization$at_bound, likelihood)
  .warn_laplace_convergence(fitted)
}
