# Covariance-kernel hyperparameter derivatives -------------------------------
#
# Derivatives are taken with respect to each parameter's optimizer coordinate
# eta_j, defined by the parameter's constraint in the kernel registry. Every
# built-in kernel parameter is positive with eta_j = log(theta_j), so
# dK / d eta_j = theta_j * dK / d theta_j. Leaf derivatives come from the
# registry entries (R/kernel-registry.R); this file combines them through
# sums, products, scale kernels, selections of input columns, which
# differentiate the selected kernel, and changepoints
# (R/kernel-changepoint.R). The derivations are in
# vignette("v07-kernel-derivatives").


.leaf_kernel_value_and_gradient <- function(kernel, x, y, prefix) {
  covariance <- .evaluate_leaf_kernel(kernel, x, y)
  derivatives <- .kernel_entry(kernel$type)$gradient(
    x,
    y,
    kernel$parameters,
    covariance
  )

  gradient <- list()

  for (parameter_name in names(kernel$parameters)) {
    derivative <- derivatives[[parameter_name]]

    if (is.matrix(derivative)) {
      derivative <- list(derivative)
    }

    names(derivative) <- .parameter_paths(
      parameter_name,
      length(kernel$parameters[[parameter_name]]),
      prefix
    )
    gradient <- c(gradient, derivative)
  }

  list(value = covariance, gradient = gradient)
}


# For each factor of a product, the product of all other factors, built from
# prefix and suffix products so that no factor is ever divided out (factors
# such as white noise contain exact zeros).
.products_of_others <- function(values) {
  n_values <- length(values)
  before <- vector("list", n_values)
  after <- vector("list", n_values)
  running <- NULL

  for (index in seq_len(n_values)) {
    before[index] <- list(running)
    running <- if (is.null(running)) {
      values[[index]]
    } else {
      running * values[[index]]
    }
  }

  running <- NULL

  for (index in rev(seq_len(n_values))) {
    after[index] <- list(running)
    running <- if (is.null(running)) {
      values[[index]]
    } else {
      values[[index]] * running
    }
  }

  Map(
    function(left, right) {
      if (is.null(left)) {
        return(right)
      }

      if (is.null(right)) {
        return(left)
      }

      left * right
    },
    before,
    after
  )
}


# Covariance matrix and named derivative matrices of a kernel tree. Values are
# combined exactly as in evaluate_kernel(); derivatives follow the sum,
# product, and scale rules, and paths follow kernel_parameters(flatten = TRUE).
.kernel_value_and_gradient <- function(kernel, x, y, prefix = "") {
  if (length(kernel$children) == 0L) {
    return(.leaf_kernel_value_and_gradient(kernel, x, y, prefix))
  }

  if (identical(kernel$type, "changepoint")) {
    return(.changepoint_value_and_gradient(kernel, x, y, prefix))
  }

  if (identical(kernel$type, "select")) {
    return(
      .kernel_value_and_gradient(
        kernel$children[[1L]],
        .select_input_columns(kernel, x),
        .select_input_columns(kernel, y),
        prefix
      )
    )
  }

  children <- Map(
    function(child, child_prefix) {
      .kernel_value_and_gradient(child, x, y, child_prefix)
    },
    kernel$children,
    .child_prefixes(kernel, prefix)
  )
  values <- lapply(children, `[[`, "value")
  gradients <- lapply(children, `[[`, "gradient")

  switch(
    kernel$type,
    sum = list(
      value = Reduce(`+`, values),
      gradient = do.call(c, gradients)
    ),
    product = {
      others <- .products_of_others(values)
      gradient <- Map(
        function(child_gradient, other) {
          if (is.null(other)) {
            return(child_gradient)
          }

          lapply(child_gradient, function(derivative) derivative * other)
        },
        gradients,
        others
      )

      list(
        value = Reduce(`*`, values),
        gradient = do.call(c, gradient)
      )
    },
    scale = {
      scale <- kernel$parameters$scale
      value <- scale * values[[1L]]
      own <- list(value)
      names(own) <- .parameter_paths("scale", 1L, prefix)

      list(
        value = value,
        gradient = c(
          own,
          lapply(gradients[[1L]], function(derivative) scale * derivative)
        )
      )
    },
    .gp_abort(
      sprintf("Unknown composite kernel type '%s'.", kernel$type),
      "gaussianprocesses_validation_error"
    )
  )
}


#' Derivatives of a covariance matrix with respect to log hyperparameters
#'
#' Computes analytical derivatives of the covariance matrix
#' \eqn{K(x, y; \theta)} with respect to the logarithm of every kernel
#' parameter, \eqn{\partial K / \partial \log \theta_j}. These are the
#' building blocks for gradient-based marginal-likelihood optimization, which
#' works in log-parameter space.
#'
#' @param kernel A Gaussian-process kernel specification.
#' @param x Numeric vector or matrix of observations.
#' @param y Optional numeric vector or matrix. If `NULL`, derivatives of the
#'   covariance of `x` with itself are returned.
#'
#' @details
#' Every kernel parameter is strictly positive, so derivatives are taken with
#' respect to \eqn{\eta_j = \log \theta_j}. The derivative with respect to the
#' parameter on its natural scale is
#' \eqn{\partial K / \partial \theta_j = \theta_j^{-1}
#' \partial K / \partial \eta_j}.
#'
#' Elements are named and ordered exactly like
#' `kernel_parameters(kernel, flatten = TRUE)`, including indexed ARD paths
#' such as `kernel1.length_scale[2]`, so they can be matched to optimizer
#' parameters by name. A scalar length scale shared by several input
#' dimensions has a single derivative.
#'
#' Leaf derivatives are closed-form expressions in the covariance matrix and
#' the scaled input differences. Sums differentiate term by term, products use
#' the product rule (without dividing by any factor, so exact zeros are safe),
#' and scaled kernels use the chain rule. The derivations are in
#' `vignette("v07-kernel-derivatives", package = "gaussianprocesses")`.
#'
#' @return A named list with one numeric matrix per flattened kernel
#'   parameter. Each matrix has `nrow(x)` rows and `nrow(y)` columns (or
#'   `nrow(x)` columns when `y` is `NULL`).
#' @seealso [kernel_input_gradient()] for derivatives with respect to the
#'   inputs.
#'
#' @examples
#' kernel <- sum_kernel(
#'   rbf_kernel(variance = 1.5, length_scale = c(0.5, 2)),
#'   white_noise_kernel(variance = 0.1)
#' )
#' x <- cbind(c(0, 0.5, 1.5), c(1, 0, 2))
#'
#' gradient <- kernel_gradient(kernel, x)
#' names(gradient)
#'
#' # Natural-scale derivative with respect to the first ARD length scale.
#' gradient[["kernel1.length_scale[1]"]] / 0.5
#' @template stable
#' @export
kernel_gradient <- function(kernel, x, y = NULL) {
  .assert_kernel(kernel)
  inputs <- .prepare_kernel_inputs(x, y)

  .kernel_value_and_gradient(kernel, inputs$x, inputs$y)$gradient
}
