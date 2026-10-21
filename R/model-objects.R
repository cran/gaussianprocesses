# Fitted model objects: schema, upgrades, and accessors -----------------------
#
# Every fitted model records the schema version of its structure and the
# package version that created it. Exported functions that accept a model
# call one of the .prepare_*() functions below first. They check the class
# and pass the model through .upgrade_model(), which converts every known
# older schema to the current one and refuses unknown or newer schemas.
#
# Schema history:
#   0  Models from before versioning (up to 0.1.0.9000 at commit 4e638aa).
#      Exact models stored the kernel and observation covariance matrices.
#   1  Exact models store only the Cholesky factor of the covariance; the
#      matrices are recomputed on demand. Every model records
#      `schema_version` and `package_version`, and exact models record the
#      numerical policy used to fit them.
#   2  Exact models record `mean_fit`, the fit of estimated or marginalized
#      mean coefficients (R/mean-coefficients.R). Schema 1 had only fixed
#      means, for which `mean_fit` is NULL, so upgrading changes nothing
#      else.
#   3  Exact models record `derivative`, the observation types of derivative
#      observations (R/derivative-observations.R); NULL when every
#      observation is a function value, as in every older model. Latent
#      models (R/gp-laplace.R) were added at this schema.
#   4  Sparse models record `method` ("fitc" or "vfe"), `trace_term`, and
#      `numerical_policy` (R/gp-sparse.R). Older sparse models are FITC; the
#      upgrade computes their trace term and leaves the policy unknown.
#      State-space models (R/gp-state-space.R) were added at this schema.
#
# Schema 4 is the format of version 1.0 and is frozen for 1.x (#39): every
# 1.x release reads schemas 0 to 4, and a later schema needs an upgrade path
# from schema 4. tests/testthat/test-model-objects.R reads models saved by
# the 0.1.0 release and schema-4 models written before 1.0
# (tests/testthat/fixtures) and checks their results.

.model_schema_version <- 4L

.supported_model_schemas <- c(0L, 1L, 2L, 3L, 4L)


.package_version <- function() {
  unname(getNamespaceVersion("gaussianprocesses"))
}


.model_descriptions <- c(
  gaussianprocesses_model = "a fitted Gaussian-process model",
  gaussianprocesses_sparse_model = "a fitted sparse Gaussian-process model",
  gaussianprocesses_heteroscedastic_model = "a fitted heteroscedastic GP model",
  gaussianprocesses_latent_model = "a fitted latent Gaussian-process model",
  gaussianprocesses_state_space_model =
    "a fitted state-space Gaussian-process model"
)


# Fields every schema of a class has, checked before any upgrade.
.model_required_fields <- list(
  gaussianprocesses_model = c(
    "x", "y", "centered_y", "kernel", "mean", "noise_variance",
    "training_noise_variance", "noise_structure", "cholesky", "alpha",
    "jitter", "n_observations", "n_features"
  ),
  gaussianprocesses_sparse_model = c(
    "x", "y", "kernel", "mean", "inducing_points", "kuu_cholesky",
    "b_cholesky", "beta", "n_observations", "n_features"
  ),
  gaussianprocesses_heteroscedastic_model = c(
    "x", "y", "mean_model", "noise_model", "kernel", "noise_kernel",
    "training_noise_variance"
  ),
  gaussianprocesses_latent_model = c(
    "x", "y", "kernel", "mean", "likelihood", "prior_mean", "mode", "alpha",
    "curvature", "b_cholesky", "jitter", "log_marginal_likelihood",
    "n_observations", "n_features"
  ),
  gaussianprocesses_state_space_model = c(
    "x", "y", "kernel", "mean", "prior_mean", "noise_variance",
    "training_noise_variance", "log_marginal_likelihood", "n_observations",
    "n_features"
  )
)


.model_format_error <- function(message, found_schema, package_version) {
  .gp_abort(
    message,
    "gaussianprocesses_model_format_error",
    found_schema = found_schema,
    supported_schemas = .supported_model_schemas,
    package_version = package_version
  )
}


# The schema of a model: 0 for models from before versioning.
.model_schema_of <- function(model) {
  version <- model[["schema_version"]]

  if (is.null(version)) {
    return(0L)
  }

  if (!is.numeric(version) || length(version) != 1L || !is.finite(version) ||
      version != round(version) || version < 1) {
    .model_format_error(
      sprintf(
        "The model has an invalid schema version (%s).",
        paste(format(version), collapse = ", ")
      ),
      found_schema = version,
      package_version = model[["package_version"]]
    )
  }

  as.integer(version)
}


.upgrade_model <- function(model, class) {
  schema <- .model_schema_of(model)
  creator <- model[["package_version"]]

  if (!(schema %in% .supported_model_schemas)) {
    .model_format_error(
      sprintf(
        paste0(
          "The model uses schema %d%s, but gaussianprocesses %s reads ",
          "schemas %s. Update the package or refit the model."
        ),
        schema,
        if (is.character(creator) && length(creator) == 1L &&
            !is.na(creator)) {
          sprintf(" (created by gaussianprocesses %s)", creator)
        } else {
          ""
        },
        .package_version(),
        paste(.supported_model_schemas, collapse = " and ")
      ),
      found_schema = schema,
      package_version = creator
    )
  }

  missing_fields <- setdiff(.model_required_fields[[class]], names(model))

  if (length(missing_fields) > 0L) {
    .model_format_error(
      sprintf(
        "The model is not a valid %s object; it has no field(s): %s.",
        class,
        paste(missing_fields, collapse = ", ")
      ),
      found_schema = schema,
      package_version = creator
    )
  }

  if (schema == .model_schema_version) {
    return(model)
  }

  .upgrade_model_from_old_schema(model, class, schema)
}


# Schema 0 to 1: exact models drop the stored covariance matrices. The
# creating package version and the numerical policy are unknown. Schemas 1
# to 3: nothing changes, because models without `mean_fit` have fixed means
# and models without `derivative` observe only function values. Schema 3 to
# 4: sparse models are FITC, and their trace term follows from the stored
# residual diagonal.
.upgrade_model_from_old_schema <- function(model, class, schema) {
  if (schema == 0L) {
    model$package_version <- NA_character_

    if (identical(class, "gaussianprocesses_model")) {
      model$kernel_covariance <- NULL
      model$observation_covariance <- NULL
    }
  }

  if (schema < 4L && identical(class, "gaussianprocesses_sparse_model")) {
    noise <- model$training_noise_variance
    model$method <- "fitc"
    model$trace_term <- if (!is.null(noise) && all(noise > 0)) {
      sum(model$residual_diagonal / noise)
    } else {
      NA_real_
    }
  }

  if (identical(class, "gaussianprocesses_heteroscedastic_model")) {
    model$mean_model <- .upgrade_model(
      model$mean_model,
      "gaussianprocesses_model"
    )
    model$noise_model <- .upgrade_model(
      model$noise_model,
      "gaussianprocesses_model"
    )
  }

  model$schema_version <- .model_schema_version
  model
}


.prepare_model <- function(model, class, name) {
  if (!inherits(model, class)) {
    .gp_abort(
      sprintf("'%s' must be %s.", name, .model_descriptions[[class]]),
      "gaussianprocesses_validation_error"
    )
  }

  .upgrade_model(model, class)
}


.prepare_gp_model <- function(model, name = "model") {
  .prepare_model(model, "gaussianprocesses_model", name)
}


.prepare_sparse_gp_model <- function(model, name = "model") {
  .prepare_model(model, "gaussianprocesses_sparse_model", name)
}


.prepare_heteroscedastic_model <- function(model, name = "model") {
  .prepare_model(model, "gaussianprocesses_heteroscedastic_model", name)
}


# Covariance matrices of exact models -----------------------------------------
#
# fit_gp() factorizes C = K + Sigma and stores only the factor. These
# accessors recompute K and C with the same code, so they reproduce the
# matrices that were factorized exactly.

.observation_covariance <- function(
    kernel_covariance,
    training_noise_variance) {
  diag(kernel_covariance) <- diag(kernel_covariance) + training_noise_variance
  kernel_covariance
}


.model_kernel_covariance <- function(model) {
  .model_typed_kernel_covariance(model)
}


.model_observation_covariance <- function(model) {
  .observation_covariance(
    .model_kernel_covariance(model),
    model$training_noise_variance
  )
}


# Covariance scale of C from its diagonal alone, without forming C.
.model_covariance_scale <- function(model) {
  diagonal <- if (is.null(model$derivative)) {
    kernel_diagonal(model$kernel, model$x)
  } else {
    .typed_diagonal(model$kernel, model$x, model$derivative)
  }

  .diagonal_scale(diagonal + model$training_noise_variance)
}


# Model specifications ---------------------------------------------------------
#
# A plain list that describes an exact model without its data or matrices:
# kernel, mean, hyperparameters, noise, numerical policy, and versions.
# Together with the training data it refits the model exactly. It is internal
# until the API review decides on its name and export.

.gp_model_specification <- function(model) {
  model <- .prepare_gp_model(model)

  list(
    schema_version = model$schema_version,
    package_version = model$package_version,
    model = "exact",
    kernel = model$kernel,
    mean = model$mean,
    derivative = model$derivative,
    hyperparameters = kernel_parameters(model$kernel, flatten = TRUE),
    noise_variance = model$noise_variance,
    noise_structure = model$noise_structure,
    numerical_policy = model$numerical_policy,
    jitter = model$jitter,
    time = if (!is.null(model$time)) {
      model$time[c("origin", "unit", "input_type")]
    },
    optimization = if (!is.null(model$optimization)) {
      model$optimization[
        c("method", "gradient", "converged", "log_marginal_likelihood")
      ]
    }
  )
}


# Refits an exact model from its specification and the training data. A
# specification of an upgraded schema-0 model has no numerical policy, so the
# defaults of fit_gp() are used.
.fit_gp_from_specification <- function(specification, x, y) {
  policy <- specification$numerical_policy

  # Models fitted before the argument was renamed record `tolerance`.
  if (!is.null(policy)) {
    names(policy)[names(policy) == "tolerance"] <- "symmetry_tolerance"
  }

  do.call(
    fit_gp,
    c(
      list(
        x = x,
        y = y,
        kernel = specification$kernel,
        noise_variance = specification$noise_variance,
        mean = specification$mean,
        derivative = specification$derivative
      ),
      policy
    )
  )
}


#' Fitted model objects
#'
#' The structure, versioning, and persistence of the models returned by
#' [fit_gp()], [optimize_gp()], [fit_time_series_gp()], [fit_sparse_gp()],
#' [fit_heteroscedastic_gp()], [fit_latent_gp()], [optimize_latent_gp()], and
#' [optimize_sparse_gp()].
#'
#' @section Exact models:
#' Objects of class `gaussianprocesses_model`, from [fit_gp()],
#' [optimize_gp()], and [fit_time_series_gp()], document these fields:
#' \describe{
#'   \item{`x`, `y`}{Training inputs as a matrix and responses.}
#'   \item{`derivative`}{For derivative observations, the type of each
#'   observation: 0 for a function value and \eqn{d} for
#'   \eqn{\partial f / \partial x_d}. `NULL` when every observation is a
#'   function value.}
#'   \item{`kernel`, `mean`}{The kernel and mean specifications. A
#'   polynomial mean records the centring and scaling of its inputs.}
#'   \item{`mean_fit`}{For estimated or marginalized mean coefficients, a
#'   list with `treatment`, `method` (`"ml"`, `"reml"`, `"gaussian"`, or
#'   `"vague"`), the `coefficients` \eqn{\bar\beta}, the basis matrix
#'   `basis` (\eqn{H}), and `coefficient_factor`, an upper-triangular
#'   \eqn{G} with \eqn{G^\top G = H^\top C^{-1} H} plus \eqn{B^{-1}} under a
#'   Gaussian prior. `NULL` for other means. Use [coef()] and [vcov()] to
#'   read the coefficients.}
#'   \item{`noise_variance`, `training_noise_variance`, `noise_structure`}{The
#'   noise variance as supplied, one value per training observation, and
#'   `"homoscedastic"` or `"known_heteroscedastic"`.}
#'   \item{`cholesky`}{Upper-triangular \eqn{R} with
#'   \eqn{R^\top R = C + \delta I}, where \eqn{C = K + \Sigma_\varepsilon} and
#'   \eqn{\delta} is `jitter`.}
#'   \item{`alpha`}{\eqn{\alpha = (C + \delta I)^{-1}(y - m(X))}, with
#'   \eqn{m(X) = H \bar\beta} for estimated or marginalized coefficients.}
#'   \item{`jitter`}{Diagonal jitter \eqn{\delta} that the factorization
#'   needed; usually 0. It is a multiple of the mean diagonal of \eqn{C}
#'   (see [fit_gp()]).}
#'   \item{`n_observations`, `n_features`}{Training size and input
#'   dimension.}
#'   \item{`optimization`}{Optimizer diagnostics, for [optimize_gp()]
#'   models.}
#'   \item{`time`}{Time metadata, for [fit_time_series_gp()] models.}
#'   \item{`schema_version`, `package_version`}{See below.}
#' }
#' The kernel matrix \eqn{K} and the observation covariance \eqn{C} are not
#' stored, because they can be recomputed exactly:
#' `evaluate_kernel(model$kernel, model$x)` is \eqn{K}, and adding
#' `diag(model$training_noise_variance)` gives \eqn{C}. Storing only
#' \eqn{R} keeps a model with \eqn{n} observations at about \eqn{8 n^2}
#' bytes. Other fields, such as `centered_y`, `prior_mean`,
#' `cholesky_attempts`, and `numerical_policy`, are internal and may change
#' without notice.
#'
#' @section Sparse and heteroscedastic models:
#' Sparse models (`gaussianprocesses_sparse_model`) document `x`, `y`,
#' `kernel`, `mean`, `noise_variance`, `inducing_points`, `n_inducing`,
#' `method` (`"fitc"` or `"vfe"`), `approximation` (`"FITC"` or `"VFE"`),
#' `log_marginal_likelihood` (the FITC objective or the VFE bound),
#' `trace_term`, and, for [optimize_sparse_gp()] models, `optimization`;
#' their other fields are internal. Heteroscedastic models
#' (`gaussianprocesses_heteroscedastic_model`) document `x`, `y`, `kernel`,
#' `noise_kernel`, `training_noise_variance`, `converged`, `iterations`, and
#' `history`, and contain two exact models, `mean_model` and `noise_model`.
#' Both record `schema_version` and `package_version`.
#'
#' @section Latent models:
#' Latent models (`gaussianprocesses_latent_model`) document the fields
#' listed in [fit_latent_gp()]. Like exact models, they store the Cholesky
#' factor of \eqn{B = I + W^{1/2} K W^{1/2}} but not \eqn{K}, and record
#' `schema_version` and `package_version`; they were added at schema 3.
#'
#' @section State-space models:
#' State-space models (`gaussianprocesses_state_space_model`), from
#' [fit_time_series_gp()] with `method = "state_space"`, store the training
#' data, kernel, mean, noise, log marginal likelihood, and the smoothed
#' fitted means and variances, but no covariance matrix or factor: every
#' prediction filters and smooths again. They were added at schema 4.
#'
#' @section Schema versions:
#' `schema_version` identifies the structure of a model, and
#' `package_version` is the version of gaussianprocesses that created it.
#' The current schema is 4: schema 2 added `mean_fit`, schema 3 added
#' `derivative`, and schema 4 added the `method` and `trace_term` of sparse
#' models. Schema 4 is the format of version 1.0. Models saved by
#' versions that predate schema versions, which stored \eqn{K} and \eqn{C},
#' are schema 0, and models saved before parametric means were added are
#' schema 1. Every function
#' that accepts a model first converts older schemas to schema 4; the
#' converted model gives the same results, and a converted sparse model is
#' a FITC model. A converted schema-0 model has
#' `package_version` `NA`, because the creating version is unknown.
#'
#' A model with a newer or unknown schema, or without the fields its class
#' requires, raises an error of class
#' `gaussianprocesses_model_format_error`. The condition object has fields
#' `found_schema`, `supported_schemas`, and `package_version`.
#'
#' @section Persistence:
#' Models can be saved with [saveRDS()] and restored with [readRDS()]; a
#' restored model gives the same predictions, likelihoods, gradients,
#' diagnostics, and seeded draws as the original.
#'
#' Every 1.x release reads models saved by any earlier version, back to
#' 0.1.0. A 1.x release may introduce a new schema only together with an
#' upgrade from schema 4, so models saved by version 1.0 stay readable and
#' give the same results. A model saved by a newer version than the one
#' reading it raises `gaussianprocesses_model_format_error` rather than
#' giving wrong results. The tests read models saved by the 0.1.0 release
#' and schema-4 models written before 1.0, and check that they give the
#' results they were saved with. Keeping the training data with the kernel,
#' mean, and noise also reproduces a model exactly: refitting them with the
#' same numerical settings gives the same model.
#'
#' @template stable
#' @name gaussianprocesses_model
#' @aliases gaussianprocesses_sparse_model
#'   gaussianprocesses_heteroscedastic_model gaussianprocesses_latent_model
#'   gaussianprocesses_state_space_model
#' @examples
#' x <- seq(-2, 2, length.out = 10)
#' model <- fit_gp(x, sin(x), rbf_kernel(), noise_variance = 0.01)
#'
#' model$schema_version
#'
#' # C is recomputed from the kernel, and R'R reproduces it.
#' C <- evaluate_kernel(model$kernel, model$x) +
#'   diag(model$training_noise_variance)
#' all.equal(crossprod(model$cholesky), C + diag(model$jitter, nrow(C)))
#'
#' # A saved model gives identical predictions.
#' path <- tempfile(fileext = ".rds")
#' saveRDS(model, path)
#' identical(predict_gp(readRDS(path), 0.5), predict_gp(model, 0.5))
NULL
