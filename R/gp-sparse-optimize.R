# Hyperparameter optimization and sampling for sparse models ----------------

#' Optimize the hyperparameters of a sparse Gaussian process
#'
#' Estimates kernel parameters, optionally the noise, and optionally the
#' inducing inputs of a sparse model ([fit_sparse_gp()]) by maximizing the
#' VFE bound or the FITC approximate marginal likelihood, with the optimizer
#' of [optimize_gp()]: L-BFGS-B on unconstrained coordinates, deterministic
#' multiple starts, and analytical gradients.
#'
#' @inheritParams fit_sparse_gp
#' @param kernel Initial kernel specification.
#' @param noise_variance Noise variance: a scalar starting value, or one value
#'   per observation forming a known pattern whose scale is estimated, as in
#'   [optimize_gp()]. With `optimize_noise = FALSE` it is fixed.
#' @param method The approximation: `"vfe"` (the default here) or `"fitc"`.
#' @param optimize_inducing If `TRUE`, optimize the inducing inputs jointly
#'   with the hyperparameters, starting from the selected or supplied
#'   points. The kernel must be differentiable in its inputs
#'   ([kernel_input_gradient()]).
#' @param optimize_noise,n_starts,start_spread,lower,upper,control,gradient
#'   As in [optimize_gp()].
#'
#' @details
#' VFE is the default because its objective is a lower bound on the exact
#' log marginal likelihood: maximizing it fits the hyperparameters of the
#' exact model as well as the inducing points allow, and pays for poor
#' coverage through the trace term. FITC's objective is not a bound; its
#' diagonal correction can absorb noise, and maximizing it tends to
#' underestimate the noise variance (Bauer, van der Wilk, and Rasmussen,
#' 2016). The research note on sparse approximations compares the two.
#'
#' The gradient is analytical for both approximations, with the inducing
#' points fixed: with \eqn{C = Q_{ff} + \Lambda} and
#' \eqn{M = C^{-1} r r^\top C^{-1} - C^{-1}},
#' \deqn{\frac{\partial}{\partial \theta_j} \log N(r \mid 0, C) =
#' \frac12 \mathrm{tr}\left(M \left(\frac{\partial Q_{ff}}{\partial \theta_j}
#' + \frac{\partial \Lambda}{\partial \theta_j}\right)\right),}
#' where \eqn{\partial Q_{ff}} follows from the derivatives of
#' \eqn{K_{uf}} and \eqn{K_{uu}}, plus the derivative of the VFE trace term.
#' It is computed in \eqn{O(n m^2)} time without forming an
#' \eqn{n \times n} matrix, processing the observations in blocks.
#'
#' @section Optimizing the inducing inputs:
#' With `optimize_inducing = TRUE` the \eqn{m \times D} inducing inputs
#' \eqn{Z} are optimized as well, on the real line, as parameters named
#' `inducing.<j>.<d>`. Their gradient follows from the same terms:
#' \eqn{Z_{jd}} changes only row \eqn{j} of \eqn{K_{uf}} and row and
#' column \eqn{j} of \eqn{K_{uu}},
#' \deqn{\frac{\partial \mathcal F}{\partial Z_{jd}} = \sum_i
#' G^{uf}_{ji} \frac{\partial k(z_j, x_i)}{\partial z_{jd}} + 2 \sum_l
#' G^{uu}_{jl} \frac{\partial k(z_j, z_l)}{\partial z_{jd}},}
#' with the kernel input derivatives of [kernel_input_gradient()] and the
#' sensitivities \eqn{G^{uf}} and \eqn{G^{uu}} of the objective to
#' \eqn{K_{uf}} and \eqn{K_{uu}}, at the same \eqn{O(n m^2)} cost.
#'
#' Optimizing \eqn{Z} is recommended with VFE: the bound can only move
#' towards the exact marginal likelihood, so better inducing inputs make the
#' approximation better. With FITC it is prone to overfitting: the inducing
#' inputs can move or collapse so that FITC's diagonal correction absorbs the
#' noise, and the noise variance is then underestimated (Bauer, van der Wilk,
#' and Rasmussen, 2016). Inducing points that become nearly identical under
#' the kernel are reported in `$inducing_collapse` with a warning of class
#' `gaussianprocesses_inducing_collapse_warning`; \eqn{K_{uu}} is factorized
#' with the package's jitter policy, so the fit stays defined. Greedy
#' variance selection (`selection = "variance"`) gives a good starting point.
#'
#' @return A fitted `gaussianprocesses_sparse_model` with optimizer
#'   diagnostics in `$optimization`, which has the fields of [optimize_gp()].
#' @seealso [fit_sparse_gp()], [optimize_gp()].
#' @examples
#' set.seed(1)
#' x <- seq(-3, 3, length.out = 300)
#' y <- sin(2 * x) + rnorm(300, sd = 0.2)
#'
#' model <- optimize_sparse_gp(x, y, rbf_kernel(), noise_variance = 0.1,
#'   n_inducing = 20, n_starts = 1)
#' model
#' model$noise_variance
#' @template experimental
#' @export
optimize_sparse_gp <- function(
    x,
    y,
    kernel,
    noise_variance = 1e-6,
    mean = zero_mean(),
    n_inducing = NULL,
    inducing_points = NULL,
    selection = c("farthest", "quantile", "variance"),
    method = c("vfe", "fitc"),
    optimize_inducing = FALSE,
    optimize_noise = TRUE,
    n_starts = 3L,
    start_spread = 1,
    lower = 1e-8,
    upper = 1e8,
    control = list(maxit = 500),
    fitc_diagonal_floor = 1e-10,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps),
    gradient = c("analytical", "numerical")) {
  data <- .validate_gp_training_data(x, y)
  .assert_kernel(kernel)
  .assert_kernel_dimensions(kernel, ncol(data$x))
  .assert_mean_function(mean)
  .assert_fixed_mean(mean, "optimize_sparse_gp()")
  mean <- .resolve_mean(mean, data$x)
  selection <- match.arg(selection)
  method <- match.arg(method)
  gradient <- match.arg(gradient)
  .assert_numeric_scalar(fitc_diagonal_floor, "fitc_diagonal_floor",
                         lower = 0, lower_open = TRUE)
  training_noise <- .resolve_training_noise_variance(
    noise_variance,
    n_observations = nrow(data$x)
  )

  for (flag in c("optimize_noise", "optimize_inducing")) {
    value <- get(flag)

    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      .gp_abort(
        sprintf("'%s' must be TRUE or FALSE.", flag),
        "gaussianprocesses_validation_error"
      )
    }
  }

  if (optimize_inducing) {
    .assert_differentiable_kernel(kernel)
  }

  n_starts <- .assert_positive_integer(n_starts, "n_starts")
  .assert_numeric_scalar(start_spread, "start_spread", lower = 0)
  bounds <- .validate_optimization_bounds(lower, upper)
  .assert_optimizer_control(control)

  noise_pattern <- if (optimize_noise && length(noise_variance) > 1L) {
    training_noise
  }

  if (!is.null(noise_pattern) && !all(noise_pattern > 0)) {
    .gp_abort(
      "A noise pattern whose scale is estimated needs positive variances.",
      "gaussianprocesses_validation_error"
    )
  }

  if (!optimize_noise) {
    .assert_vfe_noise(method, training_noise)
  }

  inducing <- .sparse_inducing_points(data$x, n_inducing, inducing_points,
                                      selection, kernel)
  noise_parameter <- .noise_parameter(noise_pattern)
  extra_initial <- if (!optimize_noise) {
    numeric(0)
  } else if (is.null(noise_pattern)) {
    c(noise_variance = max(noise_variance, bounds[["lower"]]))
  } else {
    c(noise_scale = 1)
  }
  extra_constraints <- .optimization_constraints(
    kernel,
    optimize_noise,
    noise_parameter
  )[names(extra_initial)]
  inducing_names <- .inducing_parameter_names(
    nrow(inducing$points),
    ncol(inducing$points)
  )

  if (optimize_inducing) {
    extra_initial <- c(
      extra_initial,
      stats::setNames(as.vector(inducing$points), inducing_names)
    )
    extra_constraints <- c(
      extra_constraints,
      stats::setNames(
        rep(list(.real_constraint()), length(inducing_names)),
        inducing_names
      )
    )
  }

  setup <- .optimization_setup(
    kernel,
    data$x,
    extra_initial = extra_initial,
    extra_constraints = extra_constraints,
    bounds = bounds,
    n_starts = n_starts,
    start_spread = start_spread
  )
  policy <- list(
    initial_jitter = initial_jitter,
    fallback_jitter = fallback_jitter,
    jitter_multiplier = jitter_multiplier,
    max_attempts = max_attempts,
    symmetry_tolerance = symmetry_tolerance
  )

  fit_with <- function(decoded) {
    .fit_sparse_model(
      data$x,
      data$y,
      decoded$kernel,
      decoded$noise_variance,
      .resolve_training_noise_variance(decoded$noise_variance,
                                       nrow(data$x)),
      mean,
      decoded$inducing,
      if (optimize_inducing) NULL else inducing$indices,
      if (optimize_inducing) "optimized" else inducing$method,
      method,
      fitc_diagonal_floor,
      policy
    )
  }
  # Inducing inputs, when optimized, are the last coordinates.
  decode <- function(theta) {
    is_inducing <- names(theta) %in% inducing_names
    decoded <- .decode_optimization_parameters(
      theta[!is_inducing],
      kernel = kernel,
      optimize_noise = optimize_noise,
      fixed_noise_variance = noise_variance,
      constraints = setup$constraints,
      noise_pattern = noise_pattern
    )
    decoded$inducing <- if (optimize_inducing) {
      matrix(theta[inducing_names], nrow(inducing$points),
             ncol(inducing$points))
    } else {
      inducing$points
    }
    decoded
  }

  last_theta <- NULL
  last_model <- NULL
  evaluations <- c(objective = 0L, gradient = 0L, fits = 0L)
  fit_at <- function(theta) {
    if (!identical(theta, last_theta)) {
      evaluations[["fits"]] <<- evaluations[["fits"]] + 1L
      last_model <<- tryCatch(
        fit_with(decode(theta)),
        error = function(error) NULL
      )
      last_theta <<- theta
    }

    last_model
  }
  problem <- list(
    objective = function(theta) {
      evaluations[["objective"]] <<- evaluations[["objective"]] + 1L
      model <- fit_at(theta)

      if (is.null(model) || !is.finite(model$log_marginal_likelihood)) {
        return(1e100)
      }

      -model$log_marginal_likelihood
    },
    gradient = function(theta) {
      evaluations[["gradient"]] <<- evaluations[["gradient"]] + 1L
      model <- fit_at(theta)

      if (is.null(model)) {
        return(rep(0, length(theta)))
      }

      value <- -.sparse_gradient(
        model,
        include_noise = optimize_noise,
        noise_parameter = noise_parameter,
        include_inducing = optimize_inducing
      )[names(theta)]

      if (any(!is.finite(value))) {
        return(rep(0, length(theta)))
      }

      unname(value)
    },
    evaluations = function() evaluations,
    reset = function() {
      evaluations[] <<- 0L
      invisible(NULL)
    }
  )

  run <- .run_optimization_starts(problem, setup, control, gradient)
  fitted <- fit_with(decode(run$best_theta))

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
      },
      approximation = fitted$approximation,
      optimize_inducing = optimize_inducing
    ),
    results = list(log_marginal_likelihood = fitted$log_marginal_likelihood)
  )

  .warn_inducing_collapse(fitted)
}


#' @rdname sample_gp_posterior
#' @export
sample_gp_posterior.gaussianprocesses_sparse_model <- function(
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
  model <- .prepare_sparse_gp_model(model)
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

  prediction <- predict_sparse_gp(
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
    distribution_model = model$approximation,
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
