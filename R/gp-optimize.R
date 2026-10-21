# Marginal likelihood and hyperparameter optimization -----------------------

#' Gaussian-process log marginal likelihood
#'
#' Computes the exact Gaussian log marginal likelihood from a fitted model
#' using its stored Cholesky factor and solve vector.
#'
#' @details
#' If numerical jitter was required during fitting, the stored Cholesky factor
#' corresponds to the jitter-stabilized covariance used by the numerical
#' calculation. Jitter remains distinct from statistical observation noise.
#'
#' For a mean with estimated coefficients, the value is the profile
#' likelihood, at the generalized least-squares estimate \eqn{\hat\beta}
#' (`estimate_coefficients("ml")`), or the restricted likelihood
#' (`estimate_coefficients("reml")`). For coefficients with a Gaussian prior
#' it is the marginal likelihood of the Gaussian process with covariance
#' \eqn{k(x, x') + h(x)^\top B h(x')}, and for a vague prior it is the
#' restricted likelihood. [mean_coefficients] gives the formulas.
#'
#' For a latent model from [fit_latent_gp()] or [optimize_latent_gp()] the
#' value is the Laplace approximation \eqn{\log q(y \mid X, \theta)}. For
#' a state-space model from [fit_time_series_gp()] it is the exact log
#' marginal likelihood from the Kalman filter's prediction errors.
#'
#' @param model A fitted `gaussianprocesses_model`,
#'   `gaussianprocesses_latent_model`, or
#'   `gaussianprocesses_state_space_model`.
#'
#' @return A single finite numeric log marginal likelihood.
#' @examples
#' x <- seq(-2, 2, length.out = 15)
#' y <- sin(2 * x)
#'
#' # Compare candidate length scales by their marginal likelihood.
#' sapply(c(0.2, 0.6, 2), function(length_scale) {
#'   log_marginal_likelihood(
#'     fit_gp(x, y, rbf_kernel(length_scale = length_scale), noise_variance = 0.01)
#'   )
#' })
#' @template stable
#' @export
log_marginal_likelihood <- function(model) {
  if (inherits(model, "gaussianprocesses_latent_model")) {
    return(.prepare_latent_model(model)$log_marginal_likelihood)
  }

  if (inherits(model, "gaussianprocesses_state_space_model")) {
    return(.prepare_state_space_model(model)$log_marginal_likelihood)
  }

  model <- .prepare_gp_model(model)

  quadratic <- drop(crossprod(model$centered_y, model$alpha))
  log_determinant_half <- sum(log(diag(model$cholesky)))
  n_observations <- model$n_observations
  fit <- model$mean_fit

  if (identical(.likelihood_type(model), "restricted")) {
    n_coefficients <- length(fit$coefficients)

    return(
      -0.5 * quadratic - log_determinant_half -
        0.5 * .coefficient_log_determinant(model) -
        0.5 * (n_observations - n_coefficients) * log(2 * pi)
    )
  }

  if (identical(fit$method, "gaussian")) {
    # alpha = (C + H B H')^{-1} (y - H b) under the prior N(b, B).
    residual <- model$y - drop(fit$basis %*% fit$prior$mean)

    return(
      -0.5 * sum(residual * model$alpha) - log_determinant_half -
        0.5 * fit$prior$log_determinant -
        0.5 * .coefficient_log_determinant(model) -
        0.5 * n_observations * log(2 * pi)
    )
  }

  -0.5 * quadratic - log_determinant_half - 0.5 * n_observations * log(2 * pi)
}


# d log p / d eta_j = (1/2) tr((alpha alpha' - C^{-1}) dC / d eta_j) with
# eta_j = log(theta_j). The trace needs the full precision matrix C^{-1}; it
# is formed from the stored Cholesky factor with chol2inv() and is never used
# to solve linear systems, which remain triangular solves. Restricted
# likelihoods and Gaussian coefficient priors replace C^{-1} by the
# projection P_A of R/mean-coefficients.R; profile likelihoods keep C^{-1},
# because beta_hat maximizes the likelihood (the envelope theorem).
.log_marginal_likelihood_gradient <- function(
    model,
    include_noise,
    noise_parameter = "noise_variance") {
  derivatives <- .model_covariance_gradient(model)
  precision <- .likelihood_precision(model)
  alpha <- model$alpha

  # tr(A D) = sum(A * D) because both matrices are symmetric.
  gradient <- vapply(
    derivatives,
    function(derivative) {
      0.5 * (
        drop(crossprod(alpha, derivative %*% alpha)) -
          sum(precision * derivative)
      )
    },
    numeric(1)
  )

  if (include_noise && identical(noise_parameter, "noise_variance")) {
    # dC / d log(noise variance) = noise variance * I.
    gradient <- c(
      gradient,
      noise_variance = 0.5 * model$noise_variance *
        (sum(alpha^2) - sum(diag(precision)))
    )
  }

  if (include_noise && identical(noise_parameter, "noise_scale")) {
    # Sigma = s Sigma_0, so dC / d log(s) = s Sigma_0, the training noise.
    gradient <- c(
      gradient,
      noise_scale = 0.5 * sum(
        model$training_noise_variance * (alpha^2 - diag(precision))
      )
    )
  }

  gradient
}


#' Gradient of the Gaussian-process log marginal likelihood
#'
#' Computes the exact gradient of [log_marginal_likelihood()] with respect to
#' the logarithm of every kernel parameter and, for models with a single
#' homoscedastic noise variance, the logarithm of that noise variance.
#'
#' @inheritParams log_marginal_likelihood
#'
#' @details
#' Let \eqn{C} be the observation covariance and
#' \eqn{\alpha = C^{-1}(y - m(X))}. For \eqn{\eta_j = \log \theta_j},
#' \deqn{\frac{\partial \log p(y \mid X, \theta)}{\partial \eta_j} =
#' \frac{1}{2} \mathrm{tr}\left((\alpha \alpha^\top - C^{-1})
#' \frac{\partial C}{\partial \eta_j}\right).}
#' For kernel parameters, \eqn{\partial C / \partial \eta_j} is given by
#' [kernel_gradient()]; for the noise variance \eqn{\sigma_n^2} it is
#' \eqn{\sigma_n^2 I}.
#'
#' The trace requires every entry of the precision matrix \eqn{C^{-1}}. They
#' are computed once from the stored Cholesky factor with [chol2inv()], the
#' Cholesky solve of \eqn{C X = I}; \eqn{C} is never passed to [solve()] or a
#' general-purpose inverse, and the precision matrix is not used to solve
#' linear systems. The alternative, two triangular solves per parameter,
#' never forms \eqn{C^{-1}}; against an exact reference it was equally
#' accurate for condition numbers up to \eqn{10^{12}} and 4 to 18 times
#' slower (`inst/benchmarks/gradient-numerics.R`). As in
#' [log_marginal_likelihood()], any numerical jitter used during fitting is
#' part of the factorized covariance and is treated as a constant.
#'
#' Models fitted with observation-specific noise variances have no noise
#' parameter, so only kernel derivatives are returned for them.
#'
#' For a latent model the gradient is that of the Laplace approximation,
#' including the implicit dependence of the mode on the parameters
#' ([optimize_latent_gp()]), in the optimizer coordinates of the kernel
#' parameters and of the likelihood parameters.
#'
#' @return A named numeric vector. Names follow
#'   `kernel_parameters(model$kernel, flatten = TRUE)`, followed by
#'   `noise_variance` for homoscedastic models, or by
#'   `likelihood.<name>` for the parameters of the likelihood of a latent
#'   model.
#'
#' @examples
#' x <- seq(-2, 2, length.out = 15)
#' model <- fit_gp(
#'   x,
#'   sin(x),
#'   kernel = rbf_kernel(length_scale = 0.8),
#'   noise_variance = 0.01
#' )
#' log_marginal_likelihood_gradient(model)
#' @template stable
#' @export
log_marginal_likelihood_gradient <- function(model) {
  if (inherits(model, "gaussianprocesses_latent_model")) {
    return(.laplace_gradient(.prepare_latent_model(model)))
  }

  if (inherits(model, "gaussianprocesses_state_space_model")) {
    .gp_abort(
      paste0(
        "State-space models have no analytical gradient yet; ",
        "optimize_time_series_gp() uses finite differences."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  model <- .prepare_gp_model(model)

  .log_marginal_likelihood_gradient(
    model,
    include_noise = identical(model$noise_structure, "homoscedastic")
  )
}


.validate_optimization_bounds <- function(lower, upper) {
  .assert_numeric_scalar(lower, "lower", lower = 0, lower_open = TRUE)
  .assert_numeric_scalar(upper, "upper", lower = 0, lower_open = TRUE)

  if (lower >= upper) {
    .gp_abort(
      "'lower' must be strictly smaller than 'upper'.",
      "gaussianprocesses_validation_error"
    )
  }

  invisible(c(lower = lower, upper = upper))
}


# Constraints of the optimized parameters: the flattened kernel parameters,
# followed by the noise parameter when it is optimized: the noise variance,
# or the scale of a known noise pattern.
.optimization_constraints <- function(
    kernel,
    optimize_noise,
    noise_parameter = "noise_variance") {
  constraints <- .flatten_kernel_constraints(kernel)

  if (optimize_noise) {
    noise <- list(.positive_constraint())
    names(noise) <- noise_parameter
    constraints <- c(constraints, noise)
  }

  constraints
}


.noise_parameter <- function(noise_pattern) {
  if (is.null(noise_pattern)) "noise_variance" else "noise_scale"
}


# Natural parameter values to optimizer coordinates and back, matched to
# constraints by name.
.to_coordinates <- function(values, constraints) {
  coordinates <- vapply(
    names(values),
    function(name) constraints[[name]]$to_coordinate(values[[name]]),
    numeric(1)
  )
  names(coordinates) <- names(values)
  coordinates
}


.from_coordinates <- function(coordinates, constraints) {
  values <- vapply(
    names(coordinates),
    function(name) constraints[[name]]$from_coordinate(coordinates[[name]]),
    numeric(1)
  )
  names(values) <- names(coordinates)
  values
}


# Starts in optimizer coordinates. `lower` and `upper` are coordinate bounds,
# either one per parameter or a single value shared by all. Each start after
# the first displaces one parameter by `spread` times its `scale`.
.make_optimization_starts <- function(
    initial_coordinates,
    n_starts,
    spread,
    lower,
    upper,
    scale = 1) {
  n_parameters <- length(initial_coordinates)
  lower <- rep_len(lower, n_parameters)
  upper <- rep_len(upper, n_parameters)
  scale <- rep_len(scale, n_parameters)

  starts <- matrix(
    rep(initial_coordinates, each = n_starts),
    nrow = n_starts,
    byrow = FALSE,
    dimnames = list(
      paste0("start", seq_len(n_starts)),
      names(initial_coordinates)
    )
  )

  if (n_starts == 1L || spread == 0) {
    return(starts)
  }

  for (start_index in 2:n_starts) {
    offset_index <- start_index - 2L
    parameter_index <- (offset_index %% n_parameters) + 1L
    cycle_index <- offset_index %/% n_parameters
    direction <- if ((cycle_index %% 2L) == 0L) 1 else -1
    magnitude <- spread * scale[[parameter_index]] * (1 + cycle_index %/% 2L)

    starts[start_index, parameter_index] <- min(
      upper[[parameter_index]],
      max(
        lower[[parameter_index]],
        initial_coordinates[[parameter_index]] +
          direction * magnitude
      )
    )
  }

  starts
}


# The distinct optima reached by the optimizer starts. Starts share an
# optimum when their log marginal likelihoods differ by at most
# `likelihood_tolerance` and their coordinates by at most
# `parameter_tolerance`; each optimum is represented by its best start.
# Failed starts are left out.
.distinct_optima <- function(
    results,
    likelihood_tolerance = 1e-3,
    parameter_tolerance = 1e-2) {
  objective <- vapply(results, function(result) result$value, numeric(1))
  valid <- which(is.finite(objective) & objective < 1e99)
  representatives <- integer(0)
  members <- list()

  for (index in valid[order(objective[valid])]) {
    match <- 0L

    for (cluster in seq_along(representatives)) {
      representative <- representatives[[cluster]]

      if (abs(objective[[index]] - objective[[representative]]) <=
          likelihood_tolerance &&
          max(abs(results[[index]]$par - results[[representative]]$par)) <=
            parameter_tolerance) {
        match <- cluster
        break
      }
    }

    if (match == 0L) {
      representatives <- c(representatives, index)
      members[[length(members) + 1L]] <- index
    } else {
      members[[match]] <- c(members[[match]], index)
    }
  }

  best <- -objective[[representatives[[1L]]]]

  data.frame(
    optimum = seq_along(representatives),
    starts = vapply(
      members,
      function(starts) paste(sort(starts), collapse = ", "),
      character(1)
    ),
    n_starts = lengths(members),
    log_marginal_likelihood = -objective[representatives],
    gap = best + objective[representatives],
    converged = vapply(
      representatives,
      function(index) identical(results[[index]]$convergence, 0L),
      logical(1)
    ),
    stringsAsFactors = FALSE
  )
}


# `noise_pattern` is NULL for a scalar noise variance; for a known relative
# noise pattern Sigma_0 the optimized noise is s * Sigma_0.
.decode_optimization_parameters <- function(
    theta,
    kernel,
    optimize_noise,
    fixed_noise_variance,
    constraints = .optimization_constraints(
      kernel,
      optimize_noise,
      .noise_parameter(noise_pattern)
    ),
    noise_pattern = NULL) {
  natural <- .from_coordinates(theta, constraints)

  kernel_names <- setdiff(names(natural), c("noise_variance", "noise_scale"))
  updated_kernel <- update_kernel_parameters(
    kernel,
    natural[kernel_names]
  )

  noise_variance <- if (!optimize_noise) {
    fixed_noise_variance
  } else if (is.null(noise_pattern)) {
    unname(natural[["noise_variance"]])
  } else {
    unname(natural[["noise_scale"]]) * noise_pattern
  }

  list(
    kernel = updated_kernel,
    noise_variance = noise_variance,
    parameters = natural
  )
}


# Negative log marginal likelihood and its gradient in optimizer coordinates.
# L-BFGS-B evaluates both at the same points, so the most recent fit is
# reused. Parameters for which the GP cannot be fitted get a large finite
# penalty and a zero gradient.
#
# Every call is counted, including the objective calls optim() makes for
# finite-difference gradients, which its own `counts` omit.
.optimization_problem <- function(
    x,
    y,
    kernel,
    mean,
    optimize_noise,
    fixed_noise_variance,
    initial_jitter,
    fallback_jitter,
    jitter_multiplier,
    max_attempts,
    tolerance,
    noise_pattern = NULL,
    derivative = NULL) {
  last_theta <- NULL
  last_model <- NULL
  evaluations <- c(objective = 0L, gradient = 0L, fits = 0L)
  noise_parameter <- .noise_parameter(noise_pattern)
  constraints <- .optimization_constraints(
    kernel,
    optimize_noise,
    noise_parameter
  )

  fit_at <- function(theta) {
    if (!identical(theta, last_theta)) {
      evaluations[["fits"]] <<- evaluations[["fits"]] + 1L
      decoded <- .decode_optimization_parameters(
        theta,
        kernel = kernel,
        optimize_noise = optimize_noise,
        fixed_noise_variance = fixed_noise_variance,
        constraints = constraints,
        noise_pattern = noise_pattern
      )

      last_model <<- tryCatch(
        fit_gp(
          x,
          y,
          kernel = decoded$kernel,
          noise_variance = decoded$noise_variance,
          mean = mean,
          initial_jitter = initial_jitter,
          fallback_jitter = fallback_jitter,
          jitter_multiplier = jitter_multiplier,
          max_attempts = max_attempts,
          symmetry_tolerance = tolerance,
          derivative = derivative
        ),
        error = function(error) NULL
      )
      last_theta <<- theta
    }

    last_model
  }

  objective <- function(theta) {
    evaluations[["objective"]] <<- evaluations[["objective"]] + 1L
    model <- fit_at(theta)

    if (is.null(model)) {
      return(1e100)
    }

    value <- -log_marginal_likelihood(model)

    if (!is.finite(value)) {
      return(1e100)
    }

    value
  }

  gradient <- function(theta) {
    evaluations[["gradient"]] <<- evaluations[["gradient"]] + 1L
    model <- fit_at(theta)
    failed <- rep(0, length(theta))

    if (is.null(model)) {
      return(failed)
    }

    value <- -.log_marginal_likelihood_gradient(
      model,
      include_noise = optimize_noise,
      noise_parameter = noise_parameter
    )[names(theta)]

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


.assert_optimizer_control <- function(control) {
  if (!is.list(control)) {
    .gp_abort(
      "'control' must be a list.",
      "gaussianprocesses_validation_error"
    )
  }

  invisible(control)
}


# Parameters, constraints, coordinate bounds, and starts of a hyperparameter
# optimization. The parameters are the flattened kernel parameters followed
# by `extra_initial`, such as the noise variance or likelihood parameters,
# whose constraints are `extra_constraints`. `bounds` limits the positive
# parameters on the natural scale.
.optimization_setup <- function(
    kernel,
    x,
    extra_initial,
    extra_constraints,
    bounds,
    n_starts,
    start_spread) {
  kernel_initial <- kernel_parameters(kernel, flatten = TRUE)

  if (length(kernel_initial) == 0L) {
    .gp_abort(
      "The supplied kernel has no optimizable parameters.",
      "gaussianprocesses_validation_error"
    )
  }

  constraints <- c(.flatten_kernel_constraints(kernel), extra_constraints)
  positive <- vapply(
    constraints,
    function(constraint) identical(constraint$type, "positive"),
    logical(1)
  )

  if (any(kernel_initial[positive[names(kernel_initial)]] <= 0)) {
    .gp_abort(
      "All optimized kernel parameters must be strictly positive.",
      "gaussianprocesses_validation_error"
    )
  }

  outside <- !vapply(
    names(kernel_initial),
    function(name) constraints[[name]]$contains(kernel_initial[[name]]),
    logical(1)
  )

  if (any(outside)) {
    .gp_abort(
      sprintf(
        "Kernel parameter(s) outside their constraints: %s.",
        paste(names(kernel_initial)[outside], collapse = ", ")
      ),
      "gaussianprocesses_validation_error"
    )
  }

  initial_parameters <- c(kernel_initial, extra_initial)

  # `lower` and `upper` bound the positive parameters on the natural scale.
  bounded_by_arguments <- initial_parameters[
    positive[names(initial_parameters)]
  ]

  if (any(bounded_by_arguments < bounds[["lower"]]) ||
      any(bounded_by_arguments > bounds[["upper"]])) {
    .gp_abort(
      paste0(
        "Initial optimization parameters must lie within [",
        format(bounds[["lower"]]),
        ", ",
        format(bounds[["upper"]]),
        "]."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  constraints <- constraints[names(initial_parameters)]
  coordinate_bounds <- vapply(
    constraints,
    function(constraint) {
      constraint$coordinate_bounds(bounds[["lower"]], bounds[["upper"]])
    },
    numeric(2)
  )
  lower_coordinates <- coordinate_bounds[1L, ]
  upper_coordinates <- coordinate_bounds[2L, ]
  coordinates <- vapply(
    constraints,
    function(constraint) constraint$coordinate,
    character(1)
  )

  initial_coordinates <- .to_coordinates(initial_parameters, constraints)

  # Changepoint locations are bounded by the training range of their column.
  location_columns <- .changepoint_location_columns(
    kernel,
    seq_len(ncol(x))
  )

  for (path in names(location_columns)) {
    column_range <- range(x[, location_columns[[path]]])
    lower_coordinates[[path]] <- column_range[[1L]]
    upper_coordinates[[path]] <- column_range[[2L]]
    initial_coordinates[[path]] <- min(
      column_range[[2L]],
      max(column_range[[1L]], initial_coordinates[[path]])
    )
  }

  # Spectral-mixture frequencies lie between 0 and the Nyquist frequency of
  # their input column, and later starts move them by one frequency-resolution
  # bin instead of one unit.
  frequency_limits <- .spectral_frequency_limits(
    kernel,
    x,
    seq_len(ncol(x))
  )
  start_scale <- stats::setNames(
    rep(1, length(initial_coordinates)),
    names(initial_coordinates)
  )

  for (path in names(frequency_limits)) {
    limits <- frequency_limits[[path]]
    lower_coordinates[[path]] <- 0
    upper_coordinates[[path]] <- limits$upper
    start_scale[[path]] <- limits$resolution
    initial_coordinates[[path]] <- min(
      limits$upper,
      max(0, initial_coordinates[[path]])
    )
  }

  starts <- .make_optimization_starts(
    initial_coordinates,
    n_starts = n_starts,
    spread = start_spread,
    lower = lower_coordinates,
    upper = upper_coordinates,
    scale = start_scale
  )

  # The likelihood is often multimodal in a changepoint location, so the
  # starts after the first spread it over evenly spaced quantiles.
  if (n_starts > 1L) {
    for (path in names(location_columns)) {
      starts[-1L, path] <- stats::quantile(
        x[, location_columns[[path]]],
        probs = seq_len(n_starts - 1L) / n_starts,
        names = FALSE
      )
    }
  }

  list(
    constraints = constraints,
    initial_parameters = initial_parameters,
    coordinates = coordinates,
    lower = lower_coordinates,
    upper = upper_coordinates,
    starts = starts
  )
}


# Runs L-BFGS-B from every start of `setup` on `problem`, a list of
# `objective`, `gradient`, `evaluations`, and `reset` functions such as
# .optimization_problem() returns, and picks the best converged start, or
# the best start if none converged.
.run_optimization_starts <- function(problem, setup, control, gradient) {
  starts <- setup$starts
  objective_gradient <- if (identical(gradient, "analytical")) {
    problem$gradient
  }

  results <- vector("list", nrow(starts))

  for (index in seq_len(nrow(starts))) {
    start <- starts[index, ]
    names(start) <- colnames(starts)

    problem$reset()
    result <- stats::optim(
      par = start,
      fn = problem$objective,
      gr = objective_gradient,
      method = "L-BFGS-B",
      lower = unname(setup$lower),
      upper = unname(setup$upper),
      control = control
    )

    # With numerical gradients, optim() approximates each gradient itself;
    # the objective calls that costs are included in `objective`.
    used <- problem$evaluations()
    result$evaluations <- c(
      function_evaluations = used[["objective"]],
      gradient_evaluations = if (identical(gradient, "analytical")) {
        used[["gradient"]]
      } else {
        unname(result$counts[["gradient"]])
      },
      model_fits = used[["fits"]]
    )
    results[[index]] <- result
  }

  evaluations <- Reduce(
    `+`,
    lapply(results, function(result) result$evaluations)
  )

  objective_values <- vapply(
    results,
    function(result) result$value,
    numeric(1)
  )
  valid <- is.finite(objective_values) & objective_values < 1e99

  if (!any(valid)) {
    .gp_abort(
      paste0(
        "All hyperparameter optimization starts failed to produce a ",
        "finite marginal likelihood."
      ),
      "gaussianprocesses_numerical_error"
    )
  }

  converged <- vapply(
    results,
    function(result) identical(result$convergence, 0L),
    logical(1)
  )
  preferred <- which(valid & converged)

  best_index <- if (length(preferred) > 0L) {
    preferred[[which.min(objective_values[preferred])]]
  } else {
    valid_indices <- which(valid)
    valid_indices[[which.min(objective_values[valid_indices])]]
  }

  best <- results[[best_index]]
  best_theta <- best$par
  names(best_theta) <- colnames(starts)

  list(
    results = results,
    evaluations = evaluations,
    best_index = best_index,
    best = best,
    best_theta = best_theta
  )
}


# The `$optimization` record of an optimized model. `specific` holds fields
# of one kind of model, listed after the coordinates, and `results` its
# results at the optimum, listed after the distinct optima.
.optimization_record <- function(
    setup,
    run,
    gradient,
    bounds,
    n_starts,
    start_spread,
    specific,
    results) {
  best <- run$best
  coordinates <- setup$coordinates
  constraints <- setup$constraints

  c(
    list(
      method = "L-BFGS-B",
      gradient = gradient,
      parameterization = if (all(coordinates == "log")) "log" else "mixed",
      coordinates = coordinates
    ),
    specific,
    list(
      lower = bounds[["lower"]],
      upper = bounds[["upper"]],
      n_starts = n_starts,
      start_spread = start_spread,
      initial_parameters = setup$initial_parameters,
      optimized_parameters = .from_coordinates(run$best_theta, constraints),
      best_start = run$best_index,
      convergence = best$convergence,
      converged = identical(best$convergence, 0L),
      message = if (is.null(best$message)) "" else best$message,
      counts = best$counts,
      function_evaluations = run$evaluations[["function_evaluations"]],
      gradient_evaluations = run$evaluations[["gradient_evaluations"]],
      model_fits = run$evaluations[["model_fits"]],
      objective = best$value,
      optima = .distinct_optima(run$results)
    ),
    results,
    list(
      starts = setup$starts,
      start_results = lapply(
        seq_along(run$results),
        function(index) {
          result <- run$results[[index]]
          list(
            start = index,
            convergence = result$convergence,
            converged = identical(result$convergence, 0L),
            message = if (is.null(result$message)) "" else result$message,
            objective = result$value,
            counts = result$counts,
            function_evaluations = result$evaluations[["function_evaluations"]],
            gradient_evaluations = result$evaluations[["gradient_evaluations"]],
            model_fits = result$evaluations[["model_fits"]],
            parameters = .from_coordinates(result$par, constraints)
          )
        }
      )
    )
  )
}


#' Optimize Gaussian-process hyperparameters
#'
#' Estimates kernel parameters and, optionally, observation-noise variance by
#' maximizing the exact log marginal likelihood. Optimization is performed in
#' log-parameter space with deterministic multiple starts.
#'
#' @param x Numeric vector or matrix of training inputs.
#' @param y Numeric response vector.
#' @param kernel Initial Gaussian-process kernel specification.
#' @param noise_variance Non-negative observation-noise variance: a scalar,
#'   or one value per observation. With `optimize_noise = TRUE`, a scalar is
#'   the starting value of a homoscedastic noise variance, and a vector is a
#'   known relative noise pattern \eqn{\Sigma_0} whose overall scale
#'   \eqn{s} is estimated, \eqn{\Sigma = s \Sigma_0}, starting from
#'   \eqn{s = 1}. With `optimize_noise = FALSE` the noise is fixed at the
#'   given value or values.
#' @param mean Gaussian-process mean specification. Estimated or
#'   marginalized coefficients are not optimizer parameters: for each value
#'   of the hyperparameters they are computed in closed form, and the
#'   objective is the profile, restricted, or marginal likelihood
#'   ([mean_coefficients]).
#' @param optimize_noise If `TRUE`, optimize the noise variance, or the scale
#'   of a noise pattern, jointly with the kernel parameters.
#' @param n_starts Number of deterministic optimization starts.
#' @param start_spread Log-space displacement used to construct additional
#'   deterministic starts around the supplied parameters.
#' @param lower Strictly positive lower bound applied to all optimized
#'   parameters on the natural scale.
#' @param upper Strictly positive upper bound applied to all optimized
#'   parameters on the natural scale.
#' @param control List passed to [stats::optim()]. The default allows 500
#'   L-BFGS-B iterations per start: simple kernels converge in a few dozen,
#'   but composites with periodic components can need a few hundred. Each
#'   iteration costs about one GP fit with analytical gradients.
#' @param initial_jitter Non-negative jitter tried first, relative to the
#'   covariance scale (the mean of its diagonal).
#' @param fallback_jitter Positive relative jitter tried after the first
#'   failed factorization.
#' @param jitter_multiplier Multiplicative jitter escalation factor.
#' @param max_attempts Maximum Cholesky attempts per objective evaluation.
#' @param symmetry_tolerance Relative tolerance for checking that the
#'   covariance matrix is symmetric.
#' @param derivative Optional observation types for derivative observations,
#'   as in [fit_gp()].
#' @param gradient How L-BFGS-B obtains the gradient of the objective.
#'   `"analytical"` (the default) uses [log_marginal_likelihood_gradient()],
#'   which costs one GP fit per optimizer step. `"numerical"` lets
#'   [stats::optim()] use finite differences, which costs two fits per
#'   parameter per step.
#'
#' @details
#' Kernel parameters are identified through the stable paths returned by
#' `kernel_parameters(kernel, flatten = TRUE)`. Every parameter is optimized
#' on an unconstrained coordinate determined by its constraint. Positive
#' parameters, including the noise variance, are optimized on the log scale
#' and stay positive. Real-valued parameters, such as the location of a
#' [changepoint_kernel()], are optimized on the identity scale.
#' `$optimization$coordinates` records the coordinate of each parameter. The
#' deterministic starts do not depend on `gradient`.
#'
#' A changepoint location is bounded by the training range of its input
#' column, and the starts after the first place it at evenly spaced
#' quantiles of that column, because the likelihood is often multimodal in
#' the location. A spectral-mixture frequency is bounded by 0 and the
#' Nyquist frequency of its column, and later starts displace it by one
#' frequency-resolution bin ([spectral_mixture_kernel()]).
#'
#' The objective is the likelihood that [log_marginal_likelihood()] reports
#' for the mean: `$optimization$likelihood` records which, and
#' `$optimization$mean_coefficients` the coefficients at the optimum.
#'
#' The returned diagnostics retain every start, convergence code, objective
#' value, final parameters, and the gradient method. Evaluation counts are
#' reported per start in `start_results` and summed over all starts:
#' \describe{
#'   \item{`function_evaluations`}{Evaluations of the negative log marginal
#'   likelihood, including those [stats::optim()] makes to approximate
#'   gradients by finite differences.}
#'   \item{`gradient_evaluations`}{Analytical gradient evaluations, or the
#'   number of finite-difference gradient approximations.}
#'   \item{`model_fits`}{GP fits, the dominant cost. With analytical
#'   gradients the value and gradient at the same point share one fit.}
#' }
#' `counts` is the unmodified [stats::optim()] result for the best start;
#' with numerical gradients it excludes the finite-difference evaluations.
#'
#' `$optimization$optima` groups the starts by the optimum they reached:
#' starts whose log marginal likelihoods differ by at most \eqn{10^{-3}} and
#' whose coordinates differ by at most \eqn{10^{-2}} share one. It reports,
#' for each distinct optimum, its starts, log marginal likelihood, and gap to
#' the best. Several optima with large gaps mean that the starts matter;
#' several with nearly equal likelihoods but different parameters suggest a
#' ridge of poorly identified parameters ([gp_hyperparameter_uncertainty()]).
#' [summary()] reports them.
#'
#' @return A fitted `gaussianprocesses_model` with optimized kernel and noise
#'   parameters. Optimization diagnostics are stored in `$optimization`.
#' @examples
#' x <- seq(-2, 2, length.out = 20)
#' y <- sin(2 * x) + 0.1 * cos(9 * x)
#'
#' model <- optimize_gp(x, y, kernel = rbf_kernel(), noise_variance = 0.05, n_starts = 1)
#' model
#'
#' kernel_parameters(model$kernel, flatten = TRUE)
#' model$noise_variance
#' @template stable
#' @export
optimize_gp <- function(
    x,
    y,
    kernel,
    noise_variance = 1e-6,
    mean = zero_mean(),
    optimize_noise = TRUE,
    n_starts = 3L,
    start_spread = 1,
    lower = 1e-8,
    upper = 1e8,
    control = list(maxit = 500),
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps),
    gradient = c("analytical", "numerical"),
    derivative = NULL) {
  data <- .validate_gp_training_data(x, y)
  .assert_kernel(kernel)
  .assert_kernel_dimensions(kernel, ncol(data$x))
  .assert_mean_function(mean)
  training_noise <- .resolve_training_noise_variance(
    noise_variance,
    n_observations = nrow(data$x)
  )
  gradient <- match.arg(gradient)
  derivative <- .validate_derivative_observations(derivative, data$x)

  if (!is.null(derivative)) {
    .assert_differentiable_kernel(kernel)
  }

  if (!is.logical(optimize_noise) ||
      length(optimize_noise) != 1L ||
      is.na(optimize_noise)) {
    .gp_abort(
      "'optimize_noise' must be TRUE or FALSE.",
      "gaussianprocesses_validation_error"
    )
  }

  n_starts <- .assert_positive_integer(n_starts, "n_starts")
  .assert_numeric_scalar(start_spread, "start_spread", lower = 0)
  bounds <- .validate_optimization_bounds(lower, upper)
  .assert_optimizer_control(control)

  # A vector noise with optimize_noise = TRUE is a known relative pattern
  # whose overall scale is estimated.
  noise_pattern <- if (optimize_noise && length(noise_variance) > 1L) {
    training_noise
  }

  if (!is.null(noise_pattern) && !any(noise_pattern > 0)) {
    .gp_abort(
      "A noise pattern whose scale is estimated needs a positive variance.",
      "gaussianprocesses_validation_error"
    )
  }

  noise_parameter <- .noise_parameter(noise_pattern)
  extra_initial <- if (!optimize_noise) {
    numeric(0)
  } else if (is.null(noise_pattern)) {
    c(noise_variance = max(noise_variance, bounds[["lower"]]))
  } else {
    c(noise_scale = 1)
  }

  setup <- .optimization_setup(
    kernel,
    data$x,
    extra_initial = extra_initial,
    extra_constraints = .optimization_constraints(
      kernel,
      optimize_noise,
      noise_parameter
    )[names(extra_initial)],
    bounds = bounds,
    n_starts = n_starts,
    start_spread = start_spread
  )
  constraints <- setup$constraints

  problem <- .optimization_problem(
    x = data$x,
    y = data$y,
    kernel = kernel,
    mean = mean,
    optimize_noise = optimize_noise,
    fixed_noise_variance = noise_variance,
    initial_jitter = initial_jitter,
    fallback_jitter = fallback_jitter,
    jitter_multiplier = jitter_multiplier,
    max_attempts = max_attempts,
    tolerance = symmetry_tolerance,
    noise_pattern = noise_pattern,
    derivative = derivative
  )
  run <- .run_optimization_starts(problem, setup, control, gradient)
  best_theta <- run$best_theta

  decoded <- .decode_optimization_parameters(
    best_theta,
    kernel = kernel,
    optimize_noise = optimize_noise,
    fixed_noise_variance = noise_variance,
    constraints = constraints,
    noise_pattern = noise_pattern
  )

  fitted <- fit_gp(
    data$x,
    data$y,
    kernel = decoded$kernel,
    noise_variance = decoded$noise_variance,
    mean = mean,
    initial_jitter = initial_jitter,
    fallback_jitter = fallback_jitter,
    jitter_multiplier = jitter_multiplier,
    max_attempts = max_attempts,
    symmetry_tolerance = symmetry_tolerance,
    derivative = derivative
  )

  final_log_marginal_likelihood <- log_marginal_likelihood(fitted)

  fitted$optimization <- .optimization_record(
    setup,
    run,
    gradient = gradient,
    bounds = bounds,
    n_starts = n_starts,
    start_spread = start_spread,
    specific = list(
      optimize_noise = optimize_noise,
      noise_parameterization = if (!optimize_noise) {
        "fixed"
      } else if (is.null(noise_pattern)) {
        "variance"
      } else {
        "scale"
      }
    ),
    results = list(
      likelihood = .likelihood_type(fitted),
      log_marginal_likelihood = final_log_marginal_likelihood,
      mean_coefficients = .model_coefficients(fitted)
    )
  )

  fitted
}
