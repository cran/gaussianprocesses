# Covariance kernels ----------------------------------------------------------
#
# These functions evaluate one covariance matrix directly from parameter
# values. The formulas are defined once, in the kernel registry
# (R/kernel-registry.R), and shared with kernel specifications.

#' Squared-exponential covariance kernel
#'
#' Computes the squared-exponential (radial basis function) kernel between rows
#' of `x` and `y`. A scalar length scale gives an isotropic kernel; one
#' length scale per input dimension gives automatic relevance determination
#' (ARD).
#'
#' @param x Numeric vector or matrix. Matrix rows are observations and columns
#'   are input dimensions.
#' @param y Optional numeric vector or matrix with the same number of input
#'   dimensions as `x`. If `NULL`, computes the covariance of `x` with
#'   itself.
#' @param variance Positive marginal variance.
#' @param length_scale Positive scalar or numeric vector with one value per
#'   input dimension. A vector enables ARD.
#'
#' @return A numeric covariance matrix with `nrow(x)` rows and `nrow(y)`
#'   columns. Numeric vectors are treated as one-dimensional observations.
#' @seealso [rbf_kernel()] creates the same kernel as a reusable
#'   specification, which can be combined, optimized, and used in models.
#' @examples
#' x <- seq(0, 2, by = 0.5)
#' kernel_rbf(x, variance = 2, length_scale = 0.7)
#'
#' # ARD: one length scale per input dimension.
#' x2 <- cbind(c(0, 1, 2), c(0, 0, 1))
#' kernel_rbf(x2, length_scale = c(0.5, 3))
#'
#' # Cross-covariance between two sets of inputs.
#' kernel_rbf(x, c(0.25, 1.75))
#' @template stable
#' @export
kernel_rbf <- function(x, y = NULL, variance = 1, length_scale = 1) {
  .evaluate_leaf_parameters(
    "rbf",
    list(variance = variance, length_scale = length_scale),
    x,
    y
  )
}


#' Matérn covariance kernels
#'
#' Computes Matérn kernels for smoothness parameters
#' \eqn{\nu = 1/2}, \eqn{3/2}, and \eqn{5/2}. A scalar length scale gives an
#' isotropic distance; one length scale per input dimension gives ARD.
#'
#' @inheritParams kernel_rbf
#'
#' @return A numeric covariance matrix.
#'
#' @seealso [matern_kernel_specs] for the same kernels as reusable
#'   specifications, which can be combined, optimized, and used in models.
#' @template stable
#' @name kernel_matern
#' @examples
#' x <- seq(0, 2, by = 0.5)
#'
#' # Smoothness increases from nu = 1/2 to nu = 5/2.
#' kernel_matern12(x, length_scale = 0.7)[1, ]
#' kernel_matern32(x, length_scale = 0.7)[1, ]
#' kernel_matern52(x, length_scale = 0.7)[1, ]
NULL


#' @rdname kernel_matern
#' @export
kernel_matern12 <- function(x, y = NULL, variance = 1, length_scale = 1) {
  .evaluate_leaf_parameters(
    "matern12",
    list(variance = variance, length_scale = length_scale),
    x,
    y
  )
}


#' @rdname kernel_matern
#' @export
kernel_matern32 <- function(x, y = NULL, variance = 1, length_scale = 1) {
  .evaluate_leaf_parameters(
    "matern32",
    list(variance = variance, length_scale = length_scale),
    x,
    y
  )
}


#' @rdname kernel_matern
#' @export
kernel_matern52 <- function(x, y = NULL, variance = 1, length_scale = 1) {
  .evaluate_leaf_parameters(
    "matern52",
    list(variance = variance, length_scale = length_scale),
    x,
    y
  )
}


#' Rational-quadratic covariance kernel
#'
#' The rational-quadratic kernel can be interpreted as a scale mixture of
#' squared-exponential kernels with different length scales. Scalar and ARD
#' length-scale parameterizations are both supported.
#'
#' @inheritParams kernel_rbf
#' @param alpha Positive scale-mixture parameter.
#'
#' @return A numeric covariance matrix.
#' @seealso [rational_quadratic_kernel()] creates the same kernel as a reusable
#'   specification, which can be combined, optimized, and used in models.
#' @examples
#' x <- seq(0, 2, by = 0.5)
#' kernel_rational_quadratic(x, length_scale = 0.7, alpha = 0.5)
#'
#' # As alpha grows, the kernel approaches the RBF kernel.
#' all.equal(
#'   kernel_rational_quadratic(x, alpha = 1e6),
#'   kernel_rbf(x),
#'   tolerance = 1e-5
#' )
#' @template stable
#' @export
kernel_rational_quadratic <- function(
    x,
    y = NULL,
    variance = 1,
    length_scale = 1,
    alpha = 1) {
  .evaluate_leaf_parameters(
    "rational_quadratic",
    list(variance = variance, length_scale = length_scale, alpha = alpha),
    x,
    y
  )
}


#' Periodic covariance kernel
#'
#' Computes a periodic kernel with a shared period and either an isotropic
#' length scale or one ARD length scale per input dimension.
#'
#' @inheritParams kernel_rbf
#' @param period Positive period shared across input dimensions.
#'
#' @return A numeric covariance matrix.
#' @seealso [periodic_kernel()] creates the same kernel as a reusable
#'   specification, which can be combined, optimized, and used in models.
#' @examples
#' x <- c(0, 0.25, 1, 2)
#' covariance <- kernel_periodic(x, length_scale = 0.8, period = 1)
#'
#' # Inputs a whole number of periods apart are perfectly correlated.
#' covariance[1, ]
#' @template stable
#' @export
kernel_periodic <- function(
    x,
    y = NULL,
    variance = 1,
    length_scale = 1,
    period = 1) {
  .evaluate_leaf_parameters(
    "periodic",
    list(variance = variance, length_scale = length_scale, period = period),
    x,
    y
  )
}


#' Linear covariance kernel
#'
#' Computes a linear kernel from inner products between observations.
#'
#' @inheritParams kernel_rbf
#'
#' @return A numeric covariance matrix.
#' @seealso [linear_kernel()] creates the same kernel as a reusable
#'   specification, which can be combined, optimized, and used in models.
#' @examples
#' x <- cbind(c(1, 2, 3), c(0, 1, 0))
#' kernel_linear(x, variance = 0.5)
#' @template stable
#' @export
kernel_linear <- function(x, y = NULL, variance = 1) {
  .evaluate_leaf_parameters("linear", list(variance = variance), x, y)
}


#' White-noise covariance kernel
#'
#' Returns `variance` when two observations are exactly identical and zero
#' otherwise. For `kernel_white_noise(x)`, the diagonal therefore equals the
#' requested white-noise variance.
#'
#' @inheritParams kernel_rbf
#'
#' @return A numeric covariance matrix.
#' @seealso [white_noise_kernel()] creates the same kernel as a reusable
#'   specification, which can be combined, optimized, and used in models.
#' @examples
#' # Only coincident inputs covary.
#' kernel_white_noise(c(0, 0, 1), variance = 0.1)
#' @template stable
#' @export
kernel_white_noise <- function(x, y = NULL, variance = 1) {
  .evaluate_leaf_parameters("white_noise", list(variance = variance), x, y)
}
