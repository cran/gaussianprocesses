# Mean functions --------------------------------------------------------------
#
# A mean specification has a `type`. Parametric means, m(x) = h(x)' beta,
# have a basis h and a coefficient treatment in `coefficients`:
#   fixed         known values of beta;
#   estimated     generalized least squares, with the profile or restricted
#                 likelihood for the hyperparameters;
#   marginalized  a Gaussian prior beta ~ N(b, B), or its vague limit.
# zero_mean() has no basis. constant_mean(value) with a number keeps the
# representation of earlier versions, list(parameters = list(value)), which
# saved models use. The fitting computations are in R/mean-coefficients.R.

.assert_mean_function <- function(mean_function, name = "mean") {
  if (!inherits(mean_function, "gaussianprocesses_mean")) {
    .gp_abort(
      sprintf("'%s' must be a Gaussian-process mean specification.", name),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(mean_function)
}


# Coefficient treatments -------------------------------------------------------

#' Treatments of mean-function coefficients
#'
#' A parametric mean \eqn{m(x) = h(x)^\top \beta} has coefficients
#' \eqn{\beta} that are fixed, estimated, or marginalized. Pass the result of
#' these functions as the `coefficients` argument of [linear_mean()],
#' [polynomial_mean()], [basis_mean()], or as the `value` of
#' [constant_mean()]; a numeric vector fixes the coefficients instead.
#'
#' @param method `"ml"` or `"reml"`: the likelihood that
#'   [log_marginal_likelihood()] reports and [optimize_gp()] maximizes.
#' @param mean Prior mean \eqn{b} of the coefficients: one value, recycled, or
#'   one per coefficient.
#' @param covariance Prior covariance \eqn{B}: `NULL` for the vague limit
#'   \eqn{B^{-1} \to 0}; a positive number or vector, for a diagonal
#'   covariance; or a symmetric positive-definite matrix.
#'
#' @details
#' Write \eqn{H} for the matrix of basis rows \eqn{h(x_i)^\top},
#' \eqn{C = K + \Sigma} for the observation covariance, and \eqn{p} for the
#' number of coefficients.
#'
#' **Estimated.** The coefficients are the generalized least-squares estimate
#' \eqn{\hat\beta = (H^\top C^{-1} H)^{-1} H^\top C^{-1} y}, which maximizes
#' the likelihood for given hyperparameters. Predictions plug in
#' \eqn{\hat\beta} and do not include its uncertainty. With `method = "ml"`
#' the likelihood is the profile likelihood, the Gaussian likelihood at
#' \eqn{\hat\beta}. With `method = "reml"` it is the restricted likelihood
#' \deqn{\ell_R = -\tfrac12 y^\top P y - \tfrac12 \log|C| -
#' \tfrac12 \log|H^\top C^{-1} H| - \tfrac{n - p}{2} \log 2\pi,}
#' with \eqn{P = C^{-1} - C^{-1} H (H^\top C^{-1} H)^{-1} H^\top C^{-1}},
#' which accounts for the \eqn{p} degrees of freedom used by
#' \eqn{\hat\beta} and so estimates variances with less bias.
#'
#' **Marginalized.** With the prior \eqn{\beta \sim N(b, B)}, the model is
#' exactly a Gaussian process with mean \eqn{h(x)^\top b} and covariance
#' \eqn{k(x, x') + h(x)^\top B h(x')} (Rasmussen and Williams, 2006,
#' section 2.7). Predictions use the posterior mean
#' \eqn{\bar\beta = (B^{-1} + H^\top C^{-1} H)^{-1}(H^\top C^{-1} y + B^{-1} b)}
#' and add the coefficient uncertainty
#' \eqn{R_* (B^{-1} + H^\top C^{-1} H)^{-1} R_*^\top},
#' \eqn{R_* = H_* - K_{*X} C^{-1} H}, to the predictive covariance; posterior
#' draws include it too. The likelihood is the marginal likelihood of that
#' Gaussian process. In the vague limit \eqn{\bar\beta = \hat\beta} and the
#' likelihood, up to a constant, is the restricted likelihood above.
#'
#' Rank-deficient bases make estimated coefficients, and those with a vague
#' prior, unidentifiable, and raise an error of class
#' `gaussianprocesses_rank_deficiency_error`. A proper prior identifies them.
#'
#' @return An object of class `gaussianprocesses_coefficients`.
#' @references Rasmussen, C. E. and Williams, C. K. I. (2006). *Gaussian
#'   Processes for Machine Learning*. MIT Press. Section 2.7.
#' @seealso [linear_mean()], and the universal-kriging section of
#'   `vignette("v03-marginal-likelihood-and-prediction")`.
#' @template stable
#' @name mean_coefficients
#' @examples
#' x <- seq(0, 4, length.out = 25)
#' y <- 1 + 0.5 * x + sin(2 * x)
#'
#' estimated <- fit_gp(x, y, rbf_kernel(length_scale = 0.5),
#'   noise_variance = 0.01, mean = linear_mean())
#' coef(estimated)
#'
#' # A vague prior gives the same coefficients and adds their uncertainty to
#' # predictions away from the data.
#' marginalized <- fit_gp(x, y, rbf_kernel(length_scale = 0.5),
#'   noise_variance = 0.01, mean = linear_mean(coefficient_prior()))
#' predict_gp(estimated, 8)$latent_sd
#' predict_gp(marginalized, 8)$latent_sd
NULL


#' @rdname mean_coefficients
#' @export
estimate_coefficients <- function(method = c("ml", "reml")) {
  method <- match.arg(method)

  structure(
    list(treatment = "estimated", method = method),
    class = "gaussianprocesses_coefficients"
  )
}


#' @rdname mean_coefficients
#' @export
coefficient_prior <- function(mean = 0, covariance = NULL) {
  if (!is.numeric(mean) || !is.null(dim(mean)) || length(mean) == 0L ||
      any(!is.finite(mean))) {
    .gp_abort(
      "'mean' must be a non-empty finite numeric vector.",
      "gaussianprocesses_validation_error"
    )
  }

  if (!is.null(covariance)) {
    if (is.matrix(covariance)) {
      .assert_numeric_matrix(covariance, "covariance")
      .assert_square_matrix(covariance, "covariance")

      if (!isTRUE(all.equal(covariance, t(covariance), tolerance = 1e-12)) ||
          inherits(try(chol(covariance), silent = TRUE), "try-error")) {
        .gp_abort(
          "'covariance' must be a symmetric positive-definite matrix.",
          "gaussianprocesses_validation_error"
        )
      }
    } else if (!is.numeric(covariance) || length(covariance) == 0L ||
               any(!is.finite(covariance)) || any(covariance <= 0)) {
      .gp_abort(
        paste0(
          "'covariance' must be NULL, positive variances, or a symmetric ",
          "positive-definite matrix."
        ),
        "gaussianprocesses_validation_error"
      )
    }
  }

  structure(
    list(treatment = "marginalized", mean = mean, covariance = covariance),
    class = "gaussianprocesses_coefficients"
  )
}


# A treatment from a `coefficients` argument: a treatment object, or numeric
# values that fix the coefficients.
.as_coefficient_treatment <- function(coefficients, name = "coefficients") {
  if (inherits(coefficients, "gaussianprocesses_coefficients")) {
    return(coefficients)
  }

  if (!is.numeric(coefficients) || !is.null(dim(coefficients)) ||
      length(coefficients) == 0L || any(!is.finite(coefficients))) {
    .gp_abort(
      sprintf(
        paste0(
          "'%s' must be finite numeric values, estimate_coefficients(), or ",
          "coefficient_prior()."
        ),
        name
      ),
      "gaussianprocesses_validation_error"
    )
  }

  list(treatment = "fixed", values = as.numeric(coefficients))
}


# Models other than exact ones need known mean values.
.assert_fixed_mean <- function(mean, caller) {
  if (!(.mean_treatment(mean)$treatment %in% c("none", "fixed"))) {
    .gp_abort(
      sprintf(
        paste0(
          "%s supports means with fixed coefficients only. Estimated and ",
          "marginalized coefficients need an exact model (fit_gp())."
        ),
        caller
      ),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(mean)
}


# Treatment of a mean: "none" for the zero mean, otherwise fixed, estimated,
# or marginalized.
.mean_treatment <- function(mean) {
  if (!is.null(mean$coefficients)) {
    return(mean$coefficients)
  }

  switch(
    mean$type,
    zero = list(treatment = "none"),
    constant = list(treatment = "fixed", values = mean$parameters$value),
    .gp_abort(
      sprintf("Unknown mean-function type '%s'.", mean$type),
      "gaussianprocesses_validation_error"
    )
  )
}


.new_mean <- function(type, coefficients, ...) {
  structure(
    c(
      list(type = type, parameters = list()),
      list(...),
      list(coefficients = .as_coefficient_treatment(coefficients))
    ),
    class = "gaussianprocesses_mean"
  )
}


# Constructors -----------------------------------------------------------------

#' Gaussian-process mean functions
#'
#' Mean specifications separate the deterministic part of the model from the
#' covariance kernel. `zero_mean()` is zero everywhere and `constant_mean()`
#' is a constant, fixed or with an estimated or marginalized value. For
#' regression trends see [linear_mean()], [polynomial_mean()], and
#' [basis_mean()].
#'
#' @param value A finite number that fixes the constant, or
#'   [estimate_coefficients()] or [coefficient_prior()] for a constant that
#'   is estimated from the data or has a prior (ordinary kriging).
#'
#' @return An object of class `gaussianprocesses_mean`.
#' @template stable
#' @name gp_mean_functions
#' @examples
#' zero_mean()
#' constant_mean(2.5)
#'
#' evaluate_mean(constant_mean(2.5), 1:3)
#'
#' # Estimate the constant instead.
#' x <- seq(0, 3, length.out = 15)
#' model <- fit_gp(x, 4 + sin(2 * x), rbf_kernel(), noise_variance = 0.01,
#'   mean = constant_mean(estimate_coefficients()))
#' coef(model)
NULL


#' @rdname gp_mean_functions
#' @export
zero_mean <- function() {
  structure(
    list(type = "zero", parameters = list()),
    class = "gaussianprocesses_mean"
  )
}


#' @rdname gp_mean_functions
#' @export
constant_mean <- function(value = 0) {
  if (inherits(value, "gaussianprocesses_coefficients")) {
    return(.new_mean("constant", value))
  }

  .assert_numeric_scalar(value, "value")

  structure(
    list(
      type = "constant",
      parameters = list(value = value)
    ),
    class = "gaussianprocesses_mean"
  )
}


#' Parametric mean functions
#'
#' Regression trends \eqn{m(x) = h(x)^\top \beta}, the mean of universal
#' kriging, with fixed, estimated, or marginalized coefficients.
#'
#' @param coefficients [estimate_coefficients()] (the default),
#'   [coefficient_prior()], or a numeric vector that fixes the coefficients,
#'   in the order of the basis columns.
#' @param degree Non-negative integer polynomial degree.
#' @param center,scale Optional centring and scaling of the inputs for
#'   `polynomial_mean()`: one value, recycled, or one per input column. If
#'   `NULL`, they are set from the training inputs when a model is fitted.
#' @param basis A function of the input matrix that returns a numeric matrix
#'   with one row per input and one column per coefficient.
#'
#' @details
#' The bases are:
#' \describe{
#'   \item{`linear_mean()`}{\eqn{h(x) = (1, x_1, \ldots, x_D)}, on the raw
#'   inputs, so the coefficients are an intercept and slopes. Coefficients
#'   are named `intercept`, `x1`, ..., `xD`.}
#'   \item{`polynomial_mean()`}{An intercept and the powers
#'   \eqn{z_d, z_d^2, \ldots, z_d^{k}} of each standardized input
#'   \eqn{z_d = (x_d - c_d) / s_d}, without interactions. Powers of raw
#'   inputs make badly conditioned bases, so by default \eqn{c_d} is the
#'   midpoint and \eqn{s_d} half the range of the training inputs, which
#'   maps them to \eqn{[-1, 1]}. A fitted model stores these values in its
#'   mean specification and uses them for prediction, and the coefficients,
#'   named `intercept`, `x1`, `x1^2`, ..., refer to the standardized
#'   inputs. `polynomial_mean(0)` is an unknown constant.}
#'   \item{`basis_mean()`}{The columns returned by `basis`, named by its
#'   column names or `h1`, `h2`, ....}
#' }
#'
#' [mean_coefficients] describes the three treatments of \eqn{\beta}.
#' [coef()] and [vcov()] return the coefficients of a fitted model and
#' their covariance.
#'
#' Parametric means are supported by exact models: [fit_gp()],
#' [optimize_gp()], and the functions that use them. Sparse and
#' heteroscedastic models accept them only with fixed coefficients.
#'
#' @return An object of class `gaussianprocesses_mean`.
#' @examples
#' x <- cbind(seq(0, 4, length.out = 30), rep(c(-1, 1), 15))
#' y <- 2 + 0.5 * x[, 1] - x[, 2] + sin(2 * x[, 1])
#'
#' model <- fit_gp(x, y, rbf_kernel(length_scale = c(0.5, 5)),
#'   noise_variance = 0.01, mean = linear_mean())
#' coef(model)
#' sqrt(diag(vcov(model)))
#'
#' # A quadratic trend in one input; coefficients refer to inputs scaled to
#' # [-1, 1].
#' x1 <- seq(-3, 3, length.out = 30)
#' quadratic <- fit_gp(x1, x1^2 + sin(3 * x1), rbf_kernel(length_scale = 0.4),
#'   noise_variance = 0.01, mean = polynomial_mean(2))
#' coef(quadratic)
#' quadratic$mean$scale
#'
#' # Any basis, here a periodic one with fixed coefficients.
#' seasonal <- basis_mean(
#'   function(x) cbind(sin = sin(x[, 1]), cos = cos(x[, 1])),
#'   coefficients = c(1, 0.5)
#' )
#' evaluate_mean(seasonal, c(0, pi / 2))
#' @template stable
#' @name parametric_means
NULL


#' @rdname parametric_means
#' @export
linear_mean <- function(coefficients = estimate_coefficients()) {
  .new_mean("linear", coefficients)
}


#' @rdname parametric_means
#' @export
polynomial_mean <- function(
    degree = 2L,
    coefficients = estimate_coefficients(),
    center = NULL,
    scale = NULL) {
  degree <- .assert_non_negative_integer(degree, "degree")

  if (!is.null(center) &&
      (!is.numeric(center) || length(center) == 0L ||
       any(!is.finite(center)))) {
    .gp_abort(
      "'center' must be NULL or finite numeric values.",
      "gaussianprocesses_validation_error"
    )
  }

  if (!is.null(scale) &&
      (!is.numeric(scale) || length(scale) == 0L ||
       any(!is.finite(scale)) || any(scale <= 0))) {
    .gp_abort(
      "'scale' must be NULL or positive numeric values.",
      "gaussianprocesses_validation_error"
    )
  }

  .new_mean(
    "polynomial",
    coefficients,
    degree = degree,
    center = center,
    scale = scale
  )
}


#' @rdname parametric_means
#' @export
basis_mean <- function(basis, coefficients = estimate_coefficients()) {
  if (!is.function(basis)) {
    .gp_abort(
      "'basis' must be a function of the input matrix.",
      "gaussianprocesses_validation_error"
    )
  }

  .new_mean("basis", coefficients, basis = basis)
}


# Bases ------------------------------------------------------------------------

.assert_non_negative_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0 || value != round(value)) {
    .gp_abort(
      sprintf("'%s' must be a non-negative integer.", name),
      "gaussianprocesses_validation_error"
    )
  }

  as.integer(value)
}


# One value per input column from a scalar or a vector of that length.
.per_column <- function(values, n_features, name) {
  if (length(values) == 1L) {
    return(rep(values, n_features))
  }

  if (length(values) != n_features) {
    .gp_abort(
      sprintf(
        "'%s' must have one value or one per input column (%d); it has %d.",
        name,
        n_features,
        length(values)
      ),
      "gaussianprocesses_validation_error"
    )
  }

  as.numeric(values)
}


# Fills in the standardization of a polynomial mean from training inputs:
# midpoints and half-ranges, so that the inputs map to [-1, 1]. Inputs with
# no range keep a scale of 1.
.resolve_mean <- function(mean, x) {
  if (!identical(mean$type, "polynomial")) {
    return(mean)
  }

  ranges <- apply(x, 2L, range)

  if (is.null(mean$center)) {
    mean$center <- (ranges[1L, ] + ranges[2L, ]) / 2
  }

  if (is.null(mean$scale)) {
    half_range <- (ranges[2L, ] - ranges[1L, ]) / 2
    half_range[half_range == 0] <- 1
    mean$scale <- half_range
  }

  mean
}


.standardized_inputs <- function(mean, x) {
  if (is.null(mean$center) || is.null(mean$scale)) {
    .gp_abort(
      paste0(
        "This polynomial_mean() has no 'center' and 'scale'. Fitting a model ",
        "sets them from the training inputs; supply them to evaluate the ",
        "mean directly."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  center <- .per_column(mean$center, ncol(x), "center")
  scale <- .per_column(mean$scale, ncol(x), "scale")

  sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
}


.polynomial_basis_names <- function(degree, n_features) {
  powers <- unlist(
    lapply(
      seq_len(n_features),
      function(dimension) {
        if (degree == 0L) {
          return(character(0))
        }

        paste0("x", dimension, c("", if (degree > 1L) paste0("^", 2:degree)))
      }
    )
  )

  c("intercept", powers)
}


.user_basis <- function(mean, x) {
  basis <- mean$basis(x)

  if (is.numeric(basis) && is.null(dim(basis))) {
    basis <- matrix(basis, ncol = 1L)
  }

  if (!is.matrix(basis) || !is.numeric(basis) || nrow(basis) != nrow(x) ||
      ncol(basis) == 0L || any(!is.finite(basis))) {
    .gp_abort(
      paste0(
        "The basis function must return a finite numeric matrix with one ",
        "row per input and at least one column."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (is.null(colnames(basis))) {
    colnames(basis) <- paste0("h", seq_len(ncol(basis)))
  }

  basis
}


# Basis matrix H with one row h(x_i)' per input and named columns. The zero
# mean has no columns.
.mean_basis <- function(mean, x) {
  n <- nrow(x)

  switch(
    mean$type,
    zero = matrix(0, n, 0L),
    constant = matrix(1, n, 1L, dimnames = list(NULL, "intercept")),
    linear = {
      basis <- cbind(1, x)
      dimnames(basis) <- list(
        NULL,
        c("intercept", paste0("x", seq_len(ncol(x))))
      )
      basis
    },
    polynomial = {
      z <- .standardized_inputs(mean, x)
      columns <- lapply(
        seq_len(ncol(x)),
        function(dimension) {
          outer(z[, dimension], seq_len(mean$degree), "^")
        }
      )
      basis <- do.call(cbind, c(list(rep(1, n)), columns))
      dimnames(basis) <- list(
        NULL,
        .polynomial_basis_names(mean$degree, ncol(x))
      )
      basis
    },
    basis = .user_basis(mean, x),
    .gp_abort(
      sprintf("Unknown mean-function type '%s'.", mean$type),
      "gaussianprocesses_validation_error"
    )
  )
}


# Derivatives of the basis with respect to the inputs: a list over input
# dimensions d of matrices dh(x_i)' / dx_d.
.mean_basis_gradient <- function(mean, x) {
  n <- nrow(x)
  dimensions <- seq_len(ncol(x))

  switch(
    mean$type,
    zero = lapply(dimensions, function(d) matrix(0, n, 0L)),
    constant = lapply(dimensions, function(d) matrix(0, n, 1L)),
    linear = lapply(
      dimensions,
      function(d) {
        gradient <- matrix(0, n, ncol(x) + 1L)
        gradient[, d + 1L] <- 1
        gradient
      }
    ),
    polynomial = {
      z <- .standardized_inputs(mean, x)
      scale <- .per_column(mean$scale, ncol(x), "scale")
      degree <- mean$degree
      lapply(
        dimensions,
        function(d) {
          gradient <- matrix(0, n, 1L + ncol(x) * degree)

          if (degree > 0L) {
            columns <- 1L + (d - 1L) * degree + seq_len(degree)
            gradient[, columns] <- outer(
              z[, d],
              seq_len(degree) - 1L,
              "^"
            ) %*% diag(seq_len(degree) / scale[[d]], degree)
          }

          gradient
        }
      )
    },
    .gp_abort(
      paste0(
        "basis_mean() has no derivatives with respect to the inputs, so ",
        "it cannot be used for gradient predictions or derivative ",
        "observations."
      ),
      "gaussianprocesses_validation_error"
    )
  )
}


# Fixed coefficients must match the basis.
.checked_coefficient_values <- function(values, basis) {
  if (length(values) != ncol(basis)) {
    .gp_abort(
      sprintf(
        "The mean has %d basis column(s) (%s) but %d fixed coefficient(s).",
        ncol(basis),
        paste(colnames(basis), collapse = ", "),
        length(values)
      ),
      "gaussianprocesses_validation_error"
    )
  }

  values
}


#' Evaluate a Gaussian-process mean function
#'
#' @param mean A Gaussian-process mean specification.
#' @param x Numeric vector or matrix of observations.
#'
#' @details
#' Means with fixed coefficients evaluate to \eqn{h(x)^\top \beta}. A mean
#' whose coefficients have a Gaussian prior evaluates to its prior mean
#' \eqn{h(x)^\top b}. Estimated coefficients, and those with a vague prior,
#' have no value before a model is fitted; [predict_gp()] and [coef()] use
#' the fitted values.
#'
#' @return A numeric vector with one mean value per observation.
#' @examples
#' evaluate_mean(constant_mean(1), matrix(0, nrow = 4, ncol = 2))
#' evaluate_mean(linear_mean(c(1, 2)), c(0, 0.5, 1))
#' @template stable
#' @export
evaluate_mean <- function(mean, x) {
  .assert_mean_function(mean)
  x_matrix <- .as_input_matrix(x, name = "x")
  treatment <- .mean_treatment(mean)

  if (identical(treatment$treatment, "none")) {
    return(rep(0, nrow(x_matrix)))
  }

  if (identical(mean$type, "constant") && is.null(mean$coefficients)) {
    return(rep(mean$parameters$value, nrow(x_matrix)))
  }

  basis <- .mean_basis(mean, x_matrix)

  values <- switch(
    treatment$treatment,
    fixed = .checked_coefficient_values(treatment$values, basis),
    marginalized = if (!is.null(treatment$covariance)) {
      .coefficient_prior_mean(treatment, ncol(basis))
    },
    NULL
  )

  if (is.null(values)) {
    .gp_abort(
      paste0(
        "The mean's coefficients are estimated, or have a vague prior, so ",
        "it has no value before a model is fitted. Use predict_gp() or ",
        "coef() on a fitted model."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  drop(basis %*% values)
}


# Formatting -------------------------------------------------------------------

.format_coefficient_values <- function(values) {
  formatted <- vapply(signif(values, 6), format, character(1), trim = TRUE)

  if (length(formatted) == 1L) {
    return(formatted)
  }

  paste0("c(", paste(formatted, collapse = ", "), ")")
}


.format_treatment <- function(treatment) {
  switch(
    treatment$treatment,
    fixed = paste0(
      "coefficients=",
      .format_coefficient_values(treatment$values)
    ),
    estimated = paste0("estimated by ", toupper(treatment$method)),
    marginalized = if (is.null(treatment$covariance)) {
      "marginalized, vague prior"
    } else {
      "marginalized, Gaussian prior"
    }
  )
}


.format_mean <- function(mean) {
  if (identical(mean$type, "zero")) {
    return("ZeroMean()")
  }

  if (identical(mean$type, "constant") && is.null(mean$coefficients)) {
    return(
      paste0(
        "ConstantMean(value=",
        format(signif(mean$parameters$value, 6), trim = TRUE),
        ")"
      )
    )
  }

  name <- switch(
    mean$type,
    constant = "ConstantMean",
    linear = "LinearMean",
    polynomial = "PolynomialMean",
    basis = "BasisMean",
    return(mean$type)
  )
  degree <- if (identical(mean$type, "polynomial")) {
    paste0("degree=", mean$degree, ", ")
  }

  paste0(name, "(", degree, .format_treatment(mean$coefficients), ")")
}


#' @export
print.gaussianprocesses_mean <- function(x, ...) {
  cat(.format_mean(x), "\n", sep = "")
  invisible(x)
}
