# Coefficients of parametric means --------------------------------------------
#
# For m(x) = h(x)' beta with basis matrix H (n x p) and the Cholesky factor
# C = R'R of the observation covariance, every computation uses the whitened
# basis H~ = R^{-T} H and the whitened response y~ = R^{-T} y.
#
#   Estimated, or vague prior: beta_hat solves min ||y~ - H~ beta||^2 through
#   the QR decomposition H~ = Q G.
#   Gaussian prior N(b, B) with B = U'U: beta_bar solves
#   min ||y~ - H~ beta||^2 + ||U^{-T} (beta - b)||^2 through the QR
#   decomposition of [H~; U^{-T}].
#
# Either way A = G'G equals H'C^{-1}H, plus B^{-1} under a Gaussian prior: the
# precision of beta given the hyperparameters. (H'C^{-1}H)^{-1} is never
# formed. Models store beta_bar, G, and H in `mean_fit`; `alpha` is
# C^{-1}(y - H beta_bar), which also equals P_A (y - H b) for the projection
# P_A = C^{-1} - C^{-1} H A^{-1} H' C^{-1} (with b = 0 without a prior).


.per_coefficient <- function(values, n_coefficients, name) {
  if (length(values) == 1L) {
    return(rep(as.numeric(values), n_coefficients))
  }

  if (length(values) != n_coefficients) {
    .gp_abort(
      sprintf(
        paste0(
          "The coefficient prior '%s' must have one value or one per ",
          "coefficient (%d); it has %d."
        ),
        name,
        n_coefficients,
        length(values)
      ),
      "gaussianprocesses_validation_error"
    )
  }

  as.numeric(values)
}


.coefficient_prior_mean <- function(treatment, n_coefficients) {
  .per_coefficient(treatment$mean, n_coefficients, "mean")
}


# Upper-triangular U with B = U'U.
.coefficient_prior_factor <- function(treatment, n_coefficients) {
  covariance <- treatment$covariance

  if (is.matrix(covariance)) {
    if (nrow(covariance) != n_coefficients) {
      .gp_abort(
        sprintf(
          paste0(
            "The coefficient prior covariance must be %d x %d, one row per ",
            "coefficient; it is %d x %d."
          ),
          n_coefficients,
          n_coefficients,
          nrow(covariance),
          ncol(covariance)
        ),
        "gaussianprocesses_validation_error"
      )
    }

    return(chol(covariance))
  }

  diag(
    sqrt(.per_coefficient(covariance, n_coefficients, "covariance")),
    n_coefficients
  )
}


# Prior covariance h(x)' B h(x') of the coefficient term, for prior sampling.
# NULL unless the coefficients have a Gaussian prior.
.mean_prior_covariance <- function(mean, x) {
  treatment <- .mean_treatment(mean)

  if (!identical(treatment$treatment, "marginalized") ||
      is.null(treatment$covariance)) {
    return(NULL)
  }

  basis <- .mean_basis(mean, x)
  factor <- .coefficient_prior_factor(treatment, ncol(basis))

  tcrossprod(basis %*% t(factor))
}


.rank_deficiency_error <- function(rank, basis) {
  .gp_abort(
    sprintf(
      paste0(
        "The mean basis has rank %d but %d columns (%s), so its ",
        "coefficients are not identified. Remove redundant basis columns, ",
        "lower the polynomial degree, or give the coefficients a proper ",
        "coefficient_prior()."
      ),
      rank,
      ncol(basis),
      paste(colnames(basis), collapse = ", ")
    ),
    c(
      "gaussianprocesses_rank_deficiency_error",
      "gaussianprocesses_validation_error"
    )
  )
}


# The mean at the training inputs and, for estimated or marginalized
# coefficients, their fit. `factor` is the Cholesky factor R of C, and
# `types` the observation types of derivative observations, or NULL.
.fit_mean <- function(mean, x, y, factor, types = NULL) {
  treatment <- .mean_treatment(mean)

  if (treatment$treatment %in% c("none", "fixed")) {
    return(
      list(fitted_mean = .typed_mean_values(mean, x, types), mean_fit = NULL)
    )
  }

  basis <- .typed_mean_basis(mean, x, types)
  n_coefficients <- ncol(basis)
  design <- backsolve(factor, basis, upper.tri = TRUE, transpose = TRUE)
  response <- backsolve(factor, y, upper.tri = TRUE, transpose = TRUE)
  marginalized <- identical(treatment$treatment, "marginalized")
  proper <- marginalized && !is.null(treatment$covariance)
  prior <- NULL

  if (marginalized) {
    # A vague prior does not use its mean, but it must still match the basis.
    prior_mean <- .coefficient_prior_mean(treatment, n_coefficients)
  }

  if (proper) {
    prior_factor <- .coefficient_prior_factor(treatment, n_coefficients)
    root_precision <- backsolve(
      prior_factor,
      diag(n_coefficients),
      upper.tri = TRUE,
      transpose = TRUE
    )
    design <- rbind(design, root_precision)
    response <- c(response, drop(root_precision %*% prior_mean))
    prior <- list(
      mean = prior_mean,
      log_determinant = 2 * sum(log(diag(prior_factor)))
    )
  }

  decomposition <- qr(design)

  # With full rank, qr() keeps the columns in order.
  if (decomposition$rank < n_coefficients ||
      !identical(decomposition$pivot, seq_len(n_coefficients))) {
    .rank_deficiency_error(decomposition$rank, basis)
  }

  coefficient_factor <- qr.R(decomposition)
  coefficients <- backsolve(
    coefficient_factor,
    qr.qty(decomposition, response)[seq_len(n_coefficients)]
  )
  names(coefficients) <- colnames(basis)
  dimnames(coefficient_factor) <- list(colnames(basis), colnames(basis))

  list(
    fitted_mean = drop(basis %*% coefficients),
    mean_fit = list(
      treatment = treatment$treatment,
      method = if (identical(treatment$treatment, "estimated")) {
        treatment$method
      } else if (proper) {
        "gaussian"
      } else {
        "vague"
      },
      coefficients = coefficients,
      coefficient_factor = coefficient_factor,
      basis = basis,
      prior = prior
    )
  )
}


# Which likelihood a model reports: "marginal" for fixed means and Gaussian
# coefficient priors, "profile" for coefficients estimated by ML, and
# "restricted" for REML estimates and vague priors.
.likelihood_type <- function(model) {
  fit <- model$mean_fit

  if (is.null(fit)) {
    return("marginal")
  }

  switch(
    fit$method,
    ml = "profile",
    reml = "restricted",
    vague = "restricted",
    gaussian = "marginal"
  )
}


.coefficient_log_determinant <- function(model) {
  2 * sum(log(abs(diag(model$mean_fit$coefficient_factor))))
}


# H~ = R^{-T} H.
.whitened_basis <- function(model) {
  backsolve(
    model$cholesky,
    model$mean_fit$basis,
    upper.tri = TRUE,
    transpose = TRUE
  )
}


# G = C^{-1} H A^{-1/2} = R^{-1} H~ G^{-1}, so that P_A = C^{-1} - G G'.
.coefficient_projection <- function(model) {
  orthonormal <- t(
    backsolve(
      model$mean_fit$coefficient_factor,
      t(.whitened_basis(model)),
      upper.tri = TRUE,
      transpose = TRUE
    )
  )

  backsolve(model$cholesky, orthonormal, upper.tri = TRUE)
}


# The matrix in the likelihood gradient (1/2) tr((alpha alpha' - P) dC):
# P = C^{-1} for fixed means and profile likelihoods, P_A otherwise.
.likelihood_precision <- function(model) {
  precision <- chol2inv(model$cholesky)

  if (identical(.likelihood_type(model), "profile") ||
      is.null(model$mean_fit)) {
    return(precision)
  }

  precision - tcrossprod(.coefficient_projection(model))
}


# Whether predictions add the uncertainty of the coefficients.
.marginalized_mean <- function(model) {
  identical(model$mean_fit$treatment, "marginalized")
}


# Coefficient values used by predictions: fixed values, or the fitted
# estimate or posterior mean. The zero mean has none.
.model_coefficients <- function(model) {
  if (!is.null(model$mean_fit)) {
    return(model$mean_fit$coefficients)
  }

  treatment <- .mean_treatment(model$mean)

  if (identical(treatment$treatment, "none")) {
    return(stats::setNames(numeric(0), character(0)))
  }

  names <- colnames(.mean_basis(model$mean, model$x[1L, , drop = FALSE]))
  stats::setNames(treatment$values, names)
}


# The fitted mean function at new inputs.
.model_mean <- function(model, x) {
  if (is.null(model$mean_fit)) {
    return(evaluate_mean(model$mean, x))
  }

  drop(.mean_basis(model$mean, x) %*% model$mean_fit$coefficients)
}


# Its derivatives, one column per input dimension.
.model_mean_gradient <- function(model, x) {
  coefficients <- .model_coefficients(model)
  gradients <- vapply(
    .mean_basis_gradient(model$mean, x),
    function(gradient) drop(gradient %*% coefficients),
    numeric(nrow(x))
  )

  matrix(gradients, nrow(x), ncol(x))
}


# Z = G^{-T} R_*' for the coefficient uncertainty Z'Z of marginalized means,
# with R_* = H_* - V' H~ and V = R^{-T} K(X, x*) (n x m). `basis` is H_* or,
# for gradients, dH_* / dx_d.
.coefficient_uncertainty <- function(model, projected, basis, whitened_basis) {
  backsolve(
    model$mean_fit$coefficient_factor,
    t(basis - crossprod(projected, whitened_basis)),
    upper.tri = TRUE,
    transpose = TRUE
  )
}


#' Coefficients of a Gaussian-process mean
#'
#' `coef()` returns the coefficients \eqn{\beta} of a model's mean function,
#' and `vcov()` their covariance given the hyperparameters.
#'
#' @param object A fitted `gaussianprocesses_model` or, for `coef()`,
#'   `gaussianprocesses_latent_model`.
#' @param ... Unused.
#'
#' @details
#' For estimated coefficients, `coef()` is the generalized least-squares
#' estimate \eqn{\hat\beta} and `vcov()` its covariance
#' \eqn{(H^\top C^{-1} H)^{-1}}. For marginalized coefficients they are the
#' posterior mean \eqn{\bar\beta} and covariance
#' \eqn{(B^{-1} + H^\top C^{-1} H)^{-1}}. Fixed coefficients have zero
#' covariance. Both treat the kernel and noise parameters as known. The zero
#' mean has no coefficients.
#'
#' For latent models ([fit_latent_gp()]), `coef()` returns the fixed
#' coefficients, the estimate that maximizes the approximate marginal
#' likelihood, or, under a Gaussian prior \eqn{N(b, B)}, the approximate
#' posterior mean \eqn{b + B H^\top a} with \eqn{a = \nabla \log p(y \mid
#' \hat f)}. They have no `vcov()` method.
#'
#' @return `coef()`: a named numeric vector. `vcov()`: a matrix with matching
#'   row and column names.
#' @seealso [mean_coefficients], [linear_mean()].
#' @examples
#' x <- seq(0, 4, length.out = 25)
#' model <- fit_gp(x, 1 + 0.5 * x + sin(2 * x), rbf_kernel(length_scale = 0.5),
#'   noise_variance = 0.01, mean = linear_mean())
#'
#' coef(model)
#' vcov(model)
#' @name gp_coefficients
#' @importFrom stats coef vcov
NULL


#' @rdname gp_coefficients
#' @export
coef.gaussianprocesses_model <- function(object, ...) {
  .model_coefficients(.prepare_gp_model(object))
}


#' @rdname gp_coefficients
#' @export
vcov.gaussianprocesses_model <- function(object, ...) {
  model <- .prepare_gp_model(object)
  coefficients <- .model_coefficients(model)
  names <- list(names(coefficients), names(coefficients))

  if (is.null(model$mean_fit)) {
    return(
      matrix(0, length(coefficients), length(coefficients), dimnames = names)
    )
  }

  covariance <- chol2inv(model$mean_fit$coefficient_factor)
  dimnames(covariance) <- names
  covariance
}
