# Posterior gradient of the latent function -----------------------------------

# Prior variances of the partial derivatives at each input,
# d2k(x, x) / dx_d dy_d, one column per dimension. The diagonals come from
# blocks of at most 2^20 elements in all D^2 cross derivatives, so no large
# matrix between inputs is formed.
.prior_gradient_variance <- function(kernel, x) {
  block_size <- max(1L, min(64L, 1024L %/% ncol(x)))
  variance <- matrix(0, nrow(x), ncol(x))

  for (start in seq.int(1L, nrow(x), by = block_size)) {
    rows <- start:min(start + block_size - 1L, nrow(x))
    block <- x[rows, , drop = FALSE]
    cross <- .kernel_input_derivatives(kernel, block, block, order = 2L)$cross

    for (d in seq_len(ncol(x))) {
      variance[rows, d] <- diag(cross[[d]][[d]])
    }
  }

  variance
}


#' Predict the gradient of the latent function
#'
#' Computes the posterior mean and variance of each partial derivative
#' \eqn{\partial f(x_*) / \partial x_{*,d}} of the latent function of an exact
#' Gaussian-process model. Differentiation is linear, so these derivatives
#' are jointly Gaussian with the observations and the posterior is exact.
#'
#' @param model A fitted `gaussianprocesses_model`.
#' @param x Numeric vector or matrix of prediction inputs.
#' @param interval_level Probability level for the credible intervals.
#' @param variance_tolerance Relative tolerance used when deciding whether a
#'   small negative posterior variance is floating-point error.
#' @param block_size Number of prediction inputs processed together; `NULL`
#'   chooses it from a memory budget, as in [predict_gp()].
#'
#' @details
#' With \eqn{C = R^\top R} the factorized training covariance and
#' \eqn{\alpha = C^{-1}(y - m(X))},
#' \deqn{\mathrm{E}[\partial_d f(x_*) \mid y] = \partial_d m(x_*) +
#' \partial_{x_*,d} k(x_*, X)\,\alpha,}
#' \deqn{\mathrm{Var}[\partial_d f(x_*) \mid y] =
#' \partial_d \partial'_d k(x_*, x_*) - v_d^\top v_d,
#' \qquad v_d = R^{-\top} \partial_{x_*,d} k(X, x_*).}
#' The kernel derivatives come from [kernel_input_gradient()] and its
#' second-order counterpart. Only mean-square differentiable kernels have a
#' gradient process; others raise an error of class
#' `gaussianprocesses_smoothness_error`.
#'
#' The variances are marginal: covariances between partial derivatives, or
#' between points, are not returned.
#'
#' Parametric means add \eqn{\partial h(x_*)^\top \beta / \partial x_{*,d}}
#' to the mean. Marginalized coefficients add their uncertainty to the
#' variance, as in [predict_gp()], with \eqn{H_*} and \eqn{K(x_*, X)}
#' replaced by their derivatives. [basis_mean()] has no input derivatives
#' and raises an error.
#'
#' @return An object of class `gaussianprocesses_gradient_prediction` with
#'   the prediction inputs `x` and, as matrices with one row per input and one
#'   column per input dimension, the posterior `mean`, `variance`, and `sd`
#'   of the partial derivatives and the credible-interval bounds `lower` and
#'   `upper` at `interval_level`.
#' @examples
#' x <- seq(0, 2 * pi, length.out = 25)
#' model <- fit_gp(
#'   x,
#'   sin(x),
#'   rbf_kernel(length_scale = 1),
#'   noise_variance = 1e-4
#' )
#'
#' gradient <- predict_gradient_gp(model, c(0, pi / 2, pi))
#'
#' # The derivative of sin is cos.
#' cbind(estimate = gradient$mean[, 1], truth = cos(c(0, pi / 2, pi)))
#' gradient$sd
#' @template experimental
#' @export
predict_gradient_gp <- function(
    model,
    x,
    interval_level = 0.95,
    variance_tolerance = sqrt(.Machine$double.eps),
    block_size = NULL) {
  model <- .prepare_gp_model(model)
  prediction_x <- .validate_prediction_inputs(model, x)
  .validate_interval_level(interval_level)
  .assert_numeric_scalar(variance_tolerance, "variance_tolerance", lower = 0)
  .assert_differentiable_kernel(model$kernel)

  n_predictions <- nrow(prediction_x)
  n_features <- ncol(prediction_x)
  block_size <- .prediction_block_size(
    block_size,
    n_training = model$n_observations * n_features,
    n_predictions = n_predictions
  )

  mean <- .model_mean_gradient(model, prediction_x)
  explained_variance <- matrix(0, n_predictions, n_features)
  coefficient_variance <- matrix(0, n_predictions, n_features)
  marginalized <- .marginalized_mean(model)

  if (marginalized) {
    whitened_basis <- .whitened_basis(model)
  }

  for (start in seq.int(1L, n_predictions, by = block_size)) {
    rows <- start:min(start + block_size - 1L, n_predictions)
    cross <- .model_gradient_cross_covariance(
      model,
      prediction_x[rows, , drop = FALSE]
    )

    for (d in seq_len(n_features)) {
      mean[rows, d] <- mean[rows, d] + drop(cross[[d]] %*% model$alpha)
      projected <- backsolve(
        model$cholesky,
        t(cross[[d]]),
        upper.tri = TRUE,
        transpose = TRUE
      )
      explained_variance[rows, d] <- colSums(projected^2)

      if (marginalized) {
        uncertainty <- .coefficient_uncertainty(
          model,
          projected,
          .mean_basis_gradient(
            model$mean,
            prediction_x[rows, , drop = FALSE]
          )[[d]],
          whitened_basis
        )
        coefficient_variance[rows, d] <- colSums(uncertainty^2)
      }
    }
  }

  prior_variance <- .prior_gradient_variance(model$kernel, prediction_x)
  variance <- matrix(
    .sanitize_posterior_variance(
      as.vector(prior_variance - explained_variance + coefficient_variance),
      reference_scale = max(abs(prior_variance)),
      tolerance = variance_tolerance
    ),
    n_predictions,
    n_features
  )
  sd <- sqrt(variance)
  quantile <- stats::qnorm((1 + interval_level) / 2)

  dimension_names <- colnames(prediction_x)

  if (is.null(dimension_names)) {
    dimension_names <- paste0("x", seq_len(n_features))
  }

  name_columns <- function(values) {
    colnames(values) <- dimension_names
    values
  }

  structure(
    list(
      x = prediction_x,
      mean = name_columns(mean),
      variance = name_columns(variance),
      sd = name_columns(sd),
      interval_level = interval_level,
      lower = name_columns(mean - quantile * sd),
      upper = name_columns(mean + quantile * sd)
    ),
    class = "gaussianprocesses_gradient_prediction"
  )
}


#' @export
print.gaussianprocesses_gradient_prediction <- function(x, ...) {
  cat("Gaussian-process posterior gradient\n")
  cat("  prediction points: ", nrow(x$mean), "\n", sep = "")
  cat("  input dimensions: ", ncol(x$mean), "\n", sep = "")
  cat(
    "  interval level: ",
    format(100 * x$interval_level, trim = TRUE),
    "%\n",
    sep = ""
  )
  invisible(x)
}
