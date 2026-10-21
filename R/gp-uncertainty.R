# Hyperparameter uncertainty -------------------------------------------------
#
# The estimated hyperparameters eta, on the optimizer's coordinate scale (log
# for positive parameters, identity for real ones), have the Laplace
# approximation Cov(eta_hat) ~ (-H)^{-1}, with H the Hessian of the log
# marginal likelihood at the optimum: a flat prior on the coordinate scale.
# H comes from central differences of the analytical gradient, one refit per
# displaced coordinate. Under infill asymptotics only some combinations of
# parameters are identified (for Matern kernels sigma^2 / l^(2 nu); Zhang,
# 2004); then -H is nearly singular along the rest, which its condition
# number and smallest eigenvector show.


# Which parameters a model estimated, and their constraints.
.estimated_parameter_setup <- function(model) {
  if (!is.null(model$optimization)) {
    optimized <- names(model$optimization$optimized_parameters)
    include_noise <- any(c("noise_variance", "noise_scale") %in% optimized)
    noise_parameter <- if ("noise_scale" %in% optimized) {
      "noise_scale"
    } else {
      "noise_variance"
    }
  } else {
    include_noise <- identical(model$noise_structure, "homoscedastic")
    noise_parameter <- "noise_variance"
  }

  values <- kernel_parameters(model$kernel, flatten = TRUE)

  if (include_noise) {
    # A noise scale is relative to the fitted noise, so it is 1 there.
    values <- c(
      values,
      stats::setNames(
        if (identical(noise_parameter, "noise_variance")) {
          model$noise_variance
        } else {
          1
        },
        noise_parameter
      )
    )
  }

  constraints <- .optimization_constraints(
    model$kernel,
    include_noise,
    noise_parameter
  )[names(values)]

  list(
    values = values,
    constraints = constraints,
    coordinates = .to_coordinates(values, constraints),
    include_noise = include_noise,
    noise_parameter = noise_parameter
  )
}


# The model refitted at coordinates `eta` with its data, mean, derivative
# observations, and numerical policy.
.refit_at_coordinates <- function(model, eta, setup) {
  natural <- .from_coordinates(eta, setup$constraints)
  noise_names <- c("noise_variance", "noise_scale")
  kernel <- update_kernel_parameters(
    model$kernel,
    natural[setdiff(names(natural), noise_names)]
  )
  noise_variance <- if (!setup$include_noise) {
    model$noise_variance
  } else if (identical(setup$noise_parameter, "noise_variance")) {
    unname(natural[["noise_variance"]])
  } else {
    unname(natural[["noise_scale"]]) * model$training_noise_variance
  }

  policy <- model$numerical_policy

  if (!is.null(policy)) {
    names(policy)[names(policy) == "tolerance"] <- "symmetry_tolerance"
  }

  do.call(
    fit_gp,
    c(
      list(
        x = model$x,
        y = model$y,
        kernel = kernel,
        noise_variance = noise_variance,
        mean = model$mean,
        derivative = model$derivative
      ),
      policy
    )
  )
}


.coordinate_gradient <- function(model, eta, setup) {
  refit <- .refit_at_coordinates(model, eta, setup)

  .log_marginal_likelihood_gradient(
    refit,
    include_noise = setup$include_noise,
    noise_parameter = setup$noise_parameter
  )[names(eta)]
}


# Hessian of the log marginal likelihood in the coordinates, from central
# differences of the analytical gradient, symmetrized.
.coordinate_hessian <- function(model, setup, step) {
  eta <- setup$coordinates
  columns <- lapply(
    seq_along(eta),
    function(j) {
      shift <- stats::setNames(numeric(length(eta)), names(eta))
      shift[[j]] <- step
      (.coordinate_gradient(model, eta + shift, setup) -
        .coordinate_gradient(model, eta - shift, setup)) / (2 * step)
    }
  )
  hessian <- do.call(cbind, columns)
  dimnames(hessian) <- list(names(eta), names(eta))

  (hessian + t(hessian)) / 2
}


#' Uncertainty of estimated hyperparameters
#'
#' Approximates the uncertainty of a model's hyperparameters from the
#' curvature of the log marginal likelihood at its optimum, and flags
#' hyperparameters that the data identify poorly.
#'
#' @param model A fitted `gaussianprocesses_model`, normally from
#'   [optimize_gp()].
#' @param step Step of the central differences of the gradient, on the
#'   optimizer's coordinate scale.
#' @param interval_level Level of the approximate intervals.
#' @param condition_threshold Condition number of the observed information
#'   above which the hyperparameters are flagged as poorly identified.
#' @param log_scale_threshold Standard deviation, on the log scale, of the
#'   least determined combination of log-scale parameters above which they
#'   are flagged as poorly identified. The default, \eqn{\log 2}, flags a
#'   combination known only to within a factor of two at one standard
#'   deviation.
#'
#' @details
#' The hyperparameters are those [optimize_gp()] estimated (all kernel
#' parameters and the noise for other models), on its coordinate scale
#' \eqn{\eta}: logarithms of positive parameters and real parameters
#' themselves. The observed information \eqn{I = -\nabla^2 \ell(\hat\eta)}
#' is computed by central differences of the analytical gradient
#' ([log_marginal_likelihood_gradient()]), one refit per displaced
#' coordinate. Its inverse is the covariance of the Laplace approximation
#' with a flat prior on the coordinate scale:
#' \deqn{\widehat{\mathrm{Cov}}(\hat\eta) \approx I^{-1}.}
#' Standard errors and correlations are reported on the coordinate scale,
#' and approximate intervals \eqn{\hat\eta \pm z\,\mathrm{se}} are mapped back
#' to the natural scale, so that intervals for positive parameters stay
#' positive. Estimated mean coefficients are profiled out (the likelihood
#' is the profile, restricted, or marginal likelihood of the mean), and the
#' approximation conditions on them.
#'
#' The approximation needs a maximum: `gradient` reports the gradient at the
#' estimate, which should be close to zero. Some combinations of
#' hyperparameters can be nearly unidentifiable: under infill asymptotics,
#' only \eqn{\sigma^2 / \ell^{2\nu}} is consistently estimable for a Matérn
#' kernel of smoothness \eqn{\nu} (Zhang, 2004). The likelihood is then
#' nearly flat along a ridge: the information is nearly singular in
#' absolute terms, though not always in relative ones, because the
#' curvature across the ridge can be small too. `poorly_identified` is
#' therefore `TRUE` when the information is not positive definite, when its
#' condition number exceeds `condition_threshold`, or when the standard
#' deviation along the least determined combination of the log-scale
#' parameters, `weakest_log_scale_sd`, exceeds `log_scale_threshold`.
#' `weakest_direction` is the eigenvector of the smallest eigenvalue, the
#' direction along the ridge; for a Matérn kernel under infill it is close
#' to \eqn{(2\nu, 1) / \sqrt{4\nu^2 + 1}} in (log variance, log length
#' scale). Real-valued parameters, such as changepoint locations, have no
#' common scale with the others and enter only the first two criteria.
#' [gp_profile_likelihood()] shows the likelihood along one parameter
#' without the quadratic approximation.
#'
#' @return An object of class `gaussianprocesses_hyperparameter_uncertainty`
#'   with
#' \describe{
#'   \item{`parameters`}{Data frame: `parameter`, `estimate`, `coordinate`
#'   (`"log"` or `"identity"`), `coordinate_estimate`, `standard_error` on
#'   the coordinate scale, and `lower` and `upper` on the natural scale.}
#'   \item{`covariance`, `correlation`}{On the coordinate scale; `NULL` when
#'   the information is not positive definite.}
#'   \item{`information`, `eigenvalues`, `condition_number`,
#'   `weakest_log_scale_sd`, `poorly_identified`, `weakest_direction`}{The
#'   observed information and its diagnostics.}
#'   \item{`gradient`}{The gradient at the estimate.}
#' }
#' @references Zhang, H. (2004). Inconsistent estimation and asymptotically
#'   equal interpolations in model-based geostatistics. *Journal of the
#'   American Statistical Association*, 99(465), 250--261.
#' @seealso [gp_profile_likelihood()], and the article on hyperparameter
#'   uncertainty on the documentation site.
#' @examples
#' set.seed(1)
#' x <- seq(0, 10, length.out = 60)
#' y <- sin(x) + rnorm(60, sd = 0.2)
#' model <- optimize_gp(x, y, rbf_kernel(), noise_variance = 0.1, n_starts = 1)
#'
#' uncertainty <- gp_hyperparameter_uncertainty(model)
#' uncertainty
#' uncertainty$correlation
#' @template experimental
#' @export
gp_hyperparameter_uncertainty <- function(
    model,
    step = 1e-4,
    interval_level = 0.95,
    condition_threshold = 1e6,
    log_scale_threshold = log(2)) {
  model <- .prepare_gp_model(model)
  .assert_numeric_scalar(step, "step", lower = 0, lower_open = TRUE)
  .validate_interval_level(interval_level)
  .assert_numeric_scalar(
    condition_threshold,
    "condition_threshold",
    lower = 1
  )
  .assert_numeric_scalar(
    log_scale_threshold,
    "log_scale_threshold",
    lower = 0,
    lower_open = TRUE
  )

  setup <- .estimated_parameter_setup(model)
  information <- -.coordinate_hessian(model, setup, step)
  eigen_decomposition <- eigen(information, symmetric = TRUE)
  eigenvalues <- eigen_decomposition$values
  positive_definite <- all(eigenvalues > 0)
  condition_number <- if (positive_definite) {
    max(eigenvalues) / min(eigenvalues)
  } else {
    Inf
  }

  covariance <- NULL
  correlation <- NULL
  standard_error <- rep(NA_real_, length(eigenvalues))

  if (positive_definite) {
    vectors <- eigen_decomposition$vectors
    covariance <- vectors %*% (t(vectors) / eigenvalues)
    covariance <- (covariance + t(covariance)) / 2
    dimnames(covariance) <- dimnames(information)
    standard_error <- sqrt(diag(covariance))
    correlation <- stats::cov2cor(covariance)
  }

  # Standard deviation along the least determined combination of the
  # log-scale parameters, from their block of the information.
  log_scale <- vapply(
    setup$constraints,
    function(constraint) identical(constraint$coordinate, "log"),
    logical(1)
  )
  weakest_log_scale_sd <- NA_real_

  if (any(log_scale)) {
    smallest <- min(
      eigen(
        information[log_scale, log_scale, drop = FALSE],
        symmetric = TRUE,
        only.values = TRUE
      )$values
    )
    weakest_log_scale_sd <- if (smallest > 0) 1 / sqrt(smallest) else Inf
  }

  quantile <- stats::qnorm((1 + interval_level) / 2)
  eta <- setup$coordinates
  bound <- function(sign) {
    .from_coordinates(eta + sign * quantile * standard_error, setup$constraints)
  }
  weakest <- eigen_decomposition$vectors[, length(eigenvalues)]
  # Fix the sign so that the largest entry is positive.
  weakest <- weakest * sign(weakest[which.max(abs(weakest))])
  names(weakest) <- names(eta)

  structure(
    list(
      parameters = data.frame(
        parameter = names(eta),
        estimate = unname(setup$values),
        coordinate = vapply(
          setup$constraints,
          function(constraint) constraint$coordinate,
          character(1)
        ),
        coordinate_estimate = unname(eta),
        standard_error = standard_error,
        lower = unname(bound(-1)),
        upper = unname(bound(1)),
        row.names = NULL,
        stringsAsFactors = FALSE
      ),
      covariance = covariance,
      correlation = correlation,
      information = information,
      eigenvalues = eigenvalues,
      condition_number = condition_number,
      condition_threshold = condition_threshold,
      weakest_log_scale_sd = weakest_log_scale_sd,
      log_scale_threshold = log_scale_threshold,
      poorly_identified = !positive_definite ||
        condition_number > condition_threshold ||
        isTRUE(weakest_log_scale_sd > log_scale_threshold),
      weakest_direction = weakest,
      gradient = .coordinate_gradient(model, eta, setup),
      interval_level = interval_level,
      step = step
    ),
    class = "gaussianprocesses_hyperparameter_uncertainty"
  )
}


#' @export
print.gaussianprocesses_hyperparameter_uncertainty <- function(x, ...) {
  cat(
    "Hyperparameter uncertainty (Laplace approximation on the optimizer ",
    "scale)\n",
    sep = ""
  )
  table <- x$parameters
  table$estimate <- signif(table$estimate, 5)
  table$standard_error <- signif(table$standard_error, 3)
  table$lower <- signif(table$lower, 4)
  table$upper <- signif(table$upper, 4)
  names(table)[names(table) %in% c("lower", "upper")] <- paste0(
    c("lower ", "upper "),
    format(100 * x$interval_level),
    "%"
  )
  print(
    table[, setdiff(names(table), "coordinate_estimate")],
    row.names = FALSE
  )
  cat(
    "  condition number of the information: ",
    format(signif(x$condition_number, 3)),
    "\n",
    sep = ""
  )

  if (!is.na(x$weakest_log_scale_sd)) {
    cat(
      "  weakest log-scale combination: standard deviation ",
      format(signif(x$weakest_log_scale_sd, 3)),
      "\n",
      sep = ""
    )
  }

  if (x$poorly_identified) {
    direction <- x$weakest_direction[abs(x$weakest_direction) >= 0.1]
    cat(
      "  poorly identified along: ",
      paste0(
        format(round(direction, 2), trim = TRUE),
        " ",
        names(direction),
        collapse = " + "
      ),
      " (coordinate scale)\n",
      sep = ""
    )
  }

  invisible(x)
}


#' Profile likelihood of one hyperparameter
#'
#' Fixes one estimated hyperparameter at each of a grid of values and
#' maximizes the log marginal likelihood over the others.
#'
#' @inheritParams gp_hyperparameter_uncertainty
#' @param parameter Path of the parameter, as in
#'   `kernel_parameters(model$kernel, flatten = TRUE)`, or `"noise_variance"`.
#' @param values Values of the parameter on its natural scale. If `NULL`, a
#'   grid of `n_values` points spanning three standard errors either side of
#'   the estimate on the coordinate scale, or one unit when the standard
#'   error is not available.
#' @param n_values Number of grid points when `values` is `NULL`.
#'
#' @details
#' At each value the other estimated hyperparameters start from their
#' estimates and are re-optimized with L-BFGS-B and the analytical gradient,
#' within the bounds that [optimize_gp()] used. `deviance` is twice the drop
#' in log marginal likelihood from the estimate, and `converged` whether the
#' other parameters reached a stationary point: the optimizer converged, or
#' the gradient there is below \eqn{10^{-4}}. Under the usual asymptotics
#' values with deviance below `qchisq(0.95, 1)` form an approximate 95%
#' interval. A flat profile shows a poorly identified parameter directly.
#'
#' @return A data frame of class `gaussianprocesses_profile_likelihood` with
#'   one row per value: `value`, `log_marginal_likelihood`, `deviance`, and
#'   `converged`.
#' @examples
#' set.seed(1)
#' x <- seq(0, 10, length.out = 40)
#' model <- optimize_gp(x, sin(x) + rnorm(40, sd = 0.2), rbf_kernel(),
#'   noise_variance = 0.1, n_starts = 1)
#'
#' gp_profile_likelihood(model, "length_scale", n_values = 5)
#' @template experimental
#' @export
gp_profile_likelihood <- function(
    model,
    parameter,
    values = NULL,
    n_values = 11L,
    step = 1e-4) {
  model <- .prepare_gp_model(model)
  setup <- .estimated_parameter_setup(model)

  if (!is.character(parameter) || length(parameter) != 1L ||
      !(parameter %in% names(setup$values))) {
    .gp_abort(
      sprintf(
        "'parameter' must be one of the estimated parameters: %s.",
        paste(names(setup$values), collapse = ", ")
      ),
      "gaussianprocesses_validation_error"
    )
  }

  constraint <- setup$constraints[[parameter]]
  eta <- setup$coordinates

  if (is.null(values)) {
    n_values <- .assert_positive_integer(n_values, "n_values")
    standard_error <- tryCatch(
      {
        uncertainty <- gp_hyperparameter_uncertainty(model, step = step)
        uncertainty$parameters$standard_error[
          uncertainty$parameters$parameter == parameter
        ]
      },
      error = function(error) NA_real_
    )
    half_width <- if (is.finite(standard_error)) 3 * standard_error else 1
    values <- constraint$from_coordinate(
      eta[[parameter]] + seq(-half_width, half_width, length.out = n_values)
    )
  }

  .assert_numeric_vector(values, "values")

  if (!all(constraint$contains(values))) {
    .gp_abort(
      sprintf("'values' must satisfy the constraint of '%s'.", parameter),
      "gaussianprocesses_validation_error"
    )
  }

  others <- setdiff(names(eta), parameter)
  lower <- if (!is.null(model$optimization)) model$optimization$lower else 1e-8
  upper <- if (!is.null(model$optimization)) model$optimization$upper else 1e8
  coordinate_bounds <- vapply(
    setup$constraints[others],
    function(constraint) constraint$coordinate_bounds(lower, upper),
    numeric(2)
  )

  profile <- lapply(
    values,
    function(value) {
      fixed <- eta
      fixed[[parameter]] <- constraint$to_coordinate(value)

      if (length(others) == 0L) {
        refit <- .refit_at_coordinates(model, fixed, setup)
        return(c(log_marginal_likelihood(refit), 1))
      }

      at <- function(theta) {
        point <- fixed
        point[others] <- theta
        point
      }
      objective <- function(theta) {
        refit <- tryCatch(
          .refit_at_coordinates(model, at(theta), setup),
          error = function(error) NULL
        )
        if (is.null(refit)) 1e100 else -log_marginal_likelihood(refit)
      }
      gradient <- function(theta) {
        value <- tryCatch(
          -.coordinate_gradient(model, at(theta), setup)[others],
          error = function(error) rep(0, length(others))
        )
        unname(value)
      }
      result <- stats::optim(
        par = eta[others],
        fn = objective,
        gr = gradient,
        method = "L-BFGS-B",
        lower = coordinate_bounds[1L, ],
        upper = coordinate_bounds[2L, ],
        control = list(maxit = 500)
      )

      # L-BFGS-B can stop abnormally in its line search when it starts at a
      # stationary point, such as the estimate itself; a negligible gradient
      # there still counts as converged.
      stationary <- max(abs(gradient(result$par))) <= 1e-4

      c(-result$value, identical(result$convergence, 0L) || stationary)
    }
  )
  profile <- do.call(rbind, profile)
  maximum <- log_marginal_likelihood(model)

  structure(
    data.frame(
      value = values,
      log_marginal_likelihood = profile[, 1L],
      deviance = 2 * (maximum - profile[, 1L]),
      converged = as.logical(profile[, 2L])
    ),
    parameter = parameter,
    estimate = unname(setup$values[[parameter]]),
    class = c("gaussianprocesses_profile_likelihood", "data.frame")
  )
}
