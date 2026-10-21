# Fitted model objects: schema versions, slim exact models, persistence, and
# upgrades of models saved before versioning (fixtures/legacy-models.rds,
# from a development version, and fixtures/release-0.1.0-models.rds, from
# the 0.1.0 release), and the frozen 1.0 format (fixtures/schema-4-models.rds).

legacy <- readRDS(test_path("fixtures", "legacy-models.rds"))
current_schema <- gaussianprocesses:::.model_schema_version

# The legacy fixture was computed on one platform; BLAS and math-library
# differences elsewhere change the last bits of results computed from it.
# Comparisons with it therefore allow a small relative difference, while
# results computed on the running platform are compared exactly.
cross_platform_tolerance <- 1e-10

model_x <- seq(-2, 2, length.out = 12)
model_y <- sin(2 * model_x) + 0.1 * cos(5 * model_x)
model_x_new <- c(-2.5, -0.3, 0.4, 2.2)


# One freshly fitted model of every class and construction route.
current_models <- function() {
  simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 20)

  list(
    exact = fit_gp(model_x, model_y, rbf_kernel(1.1, 0.6), 0.02),
    known_noise = fit_gp(
      model_x,
      model_y,
      matern52_kernel(1, 0.8),
      noise_variance = seq(0.01, 0.06, length.out = 12)
    ),
    optimized = optimize_gp(
      model_x,
      model_y,
      kernel = rbf_kernel(),
      noise_variance = 0.05,
      n_starts = 1
    ),
    time = fit_time_series_gp(1:12, model_y, rbf_kernel(1, 3), noise_variance = 0.02),
    sparse = fit_sparse_gp(
      model_x,
      model_y,
      rbf_kernel(1.1, 0.6),
      noise_variance = 0.02,
      n_inducing = 5
    ),
    heteroscedastic = fit_heteroscedastic_gp(
      simulation$x,
      simulation$observed,
      kernel = matern52_kernel(length_scale = 0.2),
      noise_kernel = rbf_kernel(length_scale = 0.4),
      max_iterations = 30
    )
  )
}


# What each class's exported functions return, for comparisons.
model_outputs <- function(model) {
  if (inherits(model, "gaussianprocesses_sparse_model")) {
    return(
      list(
        prediction = predict_sparse_gp(
          model,
          model_x_new,
          include_covariance = TRUE
        )
      )
    )
  }

  if (inherits(model, "gaussianprocesses_heteroscedastic_model")) {
    x_new <- c(0.1, 0.5, 0.9)
    return(
      list(
        prediction = predict_heteroscedastic_gp(
          model,
          x_new,
          include_covariance = TRUE
        ),
        samples = sample_gp_posterior(
          model,
          x_new,
          n_draws = 3,
          type = "observation",
          seed = 5
        )$observation
      )
    )
  }

  noise <- if (identical(model$noise_structure, "known_heteroscedastic")) {
    rep(0.03, length(model_x_new))
  }
  outputs <- list(
    prediction = predict_gp(
      model,
      model_x_new,
      include_covariance = TRUE,
      observation_noise_variance = noise
    ),
    log_marginal_likelihood = log_marginal_likelihood(model),
    gradient = log_marginal_likelihood_gradient(model),
    diagnostics = gp_diagnostics(model),
    samples = sample_gp_posterior(
      model,
      model_x_new,
      n_draws = 3,
      seed = 5
    )$latent
  )

  if (!is.null(model$time)) {
    outputs$forecast <- forecast_gp(model, 13:15)
  }

  outputs
}


test_that("every model records its schema and package version", {
  version <- as.character(utils::packageVersion("gaussianprocesses"))
  models <- current_models()

  for (name in names(models)) {
    expect_identical(models[[name]]$schema_version, current_schema, info = name)
    expect_identical(models[[name]]$package_version, version, info = name)
  }

  heteroscedastic <- models$heteroscedastic
  expect_identical(heteroscedastic$mean_model$schema_version, current_schema)
  expect_identical(heteroscedastic$noise_model$schema_version, current_schema)
})


test_that("exact models store one dense matrix, the Cholesky factor", {
  n <- 300
  x <- seq(0, 10, length.out = n)
  model <- fit_gp(x, sin(x), rbf_kernel(1, 0.8), noise_variance = 0.01)

  square <- vapply(
    model,
    function(field) is.matrix(field) && all(dim(field) == n),
    logical(1)
  )
  expect_identical(names(model)[square], "cholesky")
  expect_false(
    any(c("kernel_covariance", "observation_covariance") %in% names(model))
  )

  # The model is about the size of its factor, 8 n^2 bytes.
  expect_lt(
    as.numeric(utils::object.size(model)),
    1.1 * as.numeric(utils::object.size(model$cholesky))
  )
})


test_that("covariance accessors reproduce the factorized matrices exactly", {
  for (name in c("exact", "known_noise", "optimized", "time")) {
    old <- legacy$models[[name]]
    upgraded <- gaussianprocesses:::.prepare_gp_model(old)

    # The schema-0 models stored K and C; recomputing reproduces them.
    expect_equal(
      gaussianprocesses:::.model_kernel_covariance(upgraded),
      old$kernel_covariance,
      tolerance = 1e-14,
      info = name
    )
    expect_equal(
      gaussianprocesses:::.model_observation_covariance(upgraded),
      old$observation_covariance,
      tolerance = 1e-14,
      info = name
    )
  }

  model <- fit_gp(c(0, 0, 1), c(1, 1, 0), rbf_kernel(), noise_variance = 0)
  covariance <- gaussianprocesses:::.model_observation_covariance(model)
  expect_gt(model$jitter, 0)
  expect_equal(
    crossprod(model$cholesky),
    covariance + diag(model$jitter, 3),
    tolerance = 1e-14
  )
})


test_that("saved models give identical results after reading", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  models <- current_models()

  for (name in names(models)) {
    model <- models[[name]]
    saveRDS(model, path)
    restored <- readRDS(path)

    expect_identical(restored, model, info = name)
    expect_identical(model_outputs(restored), model_outputs(model), info = name)
  }
})


test_that("schema-0 models are upgraded and give the same results", {
  for (name in names(legacy$models)) {
    old <- legacy$models[[name]]

    expect_null(old$schema_version)
    expected <- legacy$outputs[[name]]
    actual <- model_outputs(old)[names(expected)]

    # gp_numerical_diagnostics() later replaced
    # jitter_relative_to_covariance_scale (jitter over the larger of 1 and the
    # largest diagonal entry) with relative_jitter and covariance_scale
    # (jitter over the mean diagonal), and later gained observation_types;
    # its other fields are unchanged.
    if (!is.null(expected$diagnostics)) {
      expected$diagnostics$numerical$jitter_relative_to_covariance_scale <- NULL
      actual$diagnostics$numerical[
        c("relative_jitter", "covariance_scale", "observation_types")
      ] <- NULL
    }

    # Heteroscedastic predictions later gained the shared prediction fields
    # observation_noise_variance and observation_covariance; the others are
    # unchanged.
    # They also lost `noise_variance` in 1.0.0 (#40), the copy of
    # observation_noise_variance kept until then.
    if (identical(name, "heteroscedastic")) {
      expect_equal(actual$prediction$observation_noise_variance,
                   expected$prediction$noise_variance,
                   tolerance = cross_platform_tolerance)
      expected$prediction$noise_variance <- NULL
      actual$prediction[
        c("observation_noise_variance", "observation_covariance")
      ] <- NULL
    }

    # Sparse predictions later gained `approximation`, the label of the
    # sparse method.
    if (identical(name, "sparse")) {
      expect_identical(actual$prediction$approximation, "FITC")
      actual$prediction$approximation <- NULL
    }

    # At an optimum the gradient is a near-cancelling difference of terms of
    # order one, so last-bit differences in those terms appear as absolute
    # differences of about 1e-11 in a result of about 1e-7. It is compared on
    # an absolute scale.
    if (!is.null(expected$gradient)) {
      expect_lt(
        max(abs(actual$gradient - expected$gradient)),
        1e-8,
        label = paste(name, "gradient")
      )
      expect_identical(names(actual$gradient), names(expected$gradient))
      actual$gradient <- NULL
      expected$gradient <- NULL
    }

    expect_equal(
      actual,
      expected,
      tolerance = cross_platform_tolerance,
      info = name
    )
  }

  upgraded <- gaussianprocesses:::.prepare_gp_model(legacy$models$exact)
  expect_identical(upgraded$schema_version, current_schema)
  expect_identical(upgraded$package_version, NA_character_)
  expect_false(
    any(c("kernel_covariance", "observation_covariance") %in% names(upgraded))
  )

  heteroscedastic <- gaussianprocesses:::.prepare_heteroscedastic_model(
    legacy$models$heteroscedastic
  )
  expect_identical(heteroscedastic$schema_version, current_schema)
  expect_null(heteroscedastic$mean_model$observation_covariance)
  expect_identical(heteroscedastic$noise_model$schema_version, current_schema)
})


# Every exported function that accepts a model, by class.
model_entry_points <- list(
  gaussianprocesses_model = list(
    predict_gp = function(model) predict_gp(model, 0.5),
    predict_gradient_gp = function(model) predict_gradient_gp(model, 0.5),
    log_marginal_likelihood = log_marginal_likelihood,
    log_marginal_likelihood_gradient = log_marginal_likelihood_gradient,
    loo_gp = loo_gp,
    gp_calibration_diagnostics = gp_calibration_diagnostics,
    gp_numerical_diagnostics = gp_numerical_diagnostics,
    gp_diagnostics = gp_diagnostics,
    sample_gp_posterior = function(model) {
      sample_gp_posterior(model, 0.5, seed = 1)
    },
    print = function(model) utils::capture.output(print(model))
  ),
  time = list(
    forecast_gp = function(model) forecast_gp(model, 13)
  ),
  gaussianprocesses_sparse_model = list(
    predict_sparse_gp = function(model) predict_sparse_gp(model, 0.5),
    print = function(model) utils::capture.output(print(model))
  ),
  gaussianprocesses_heteroscedastic_model = list(
    predict_heteroscedastic_gp = function(model) {
      predict_heteroscedastic_gp(model, 0.5)
    },
    sample_gp_posterior = function(model) {
      sample_gp_posterior(model, 0.5, seed = 1)
    },
    print = function(model) utils::capture.output(print(model))
  )
)

entry_point_cases <- function() {
  list(
    list(model = legacy$models$exact, functions = model_entry_points[[1]]),
    list(model = legacy$models$time, functions = model_entry_points$time),
    list(model = legacy$models$sparse, functions = model_entry_points[[3]]),
    list(
      model = legacy$models$heteroscedastic,
      functions = model_entry_points[[4]]
    )
  )
}


test_that("the model format is frozen at schema 4 for 1.x", {
  # Schema 4 is the format of version 1.0 (#39). A later schema needs an
  # upgrade path from schema 4, and every 1.x release must keep reading the
  # fixtures of 0.1.0 and schema-4 models below; see the versioning policy
  # in CONTRIBUTING.md.
  frozen <- readRDS(test_path("fixtures", "schema-4-models.rds"))

  expect_identical(current_schema, 4L)
  expect_identical(gaussianprocesses:::.supported_model_schemas, 0:4)
  expect_identical(frozen$schema, current_schema)
})


# Fields of `expected` that `actual` must still have, compared one level
# deep: results keep their class and may gain fields but never lose them.
expect_same_fields <- function(actual, expected, label) {
  expect_identical(class(actual), class(expected),
                   label = paste(label, "class"))
  # Compare the fields only: data frames also carry row names.
  fields <- function(result) {
    result <- unclass(result)
    attributes(result) <- list(names = names(result))
    result
  }
  actual <- fields(actual)
  expected <- fields(expected)

  # Fields removed in 1.0.0 (#40) are compared with their replacements.
  removed <- gaussianprocesses:::.removed_in_1.0.0$fields
  for (field in intersect(names(removed), names(expected))) {
    expect_equal(actual[[removed[[field]]]], expected[[field]],
                 tolerance = cross_platform_tolerance,
                 label = paste(label, field))
    expected[[field]] <- NULL
  }

  expect_identical(setdiff(names(expected), names(actual)), character(0),
                   label = label)
  expect_equal(actual[names(expected)], expected,
               tolerance = cross_platform_tolerance, label = label)
}


test_that("models saved by release 0.1.0 are read and give its results", {
  release <- readRDS(test_path("fixtures", "release-0.1.0-models.rds"))
  x_new <- release$x_new

  for (name in names(release$models)) {
    model <- release$models[[name]]
    expected <- release$outputs[[name]]

    expect_null(model$schema_version)

    if (inherits(model, "gaussianprocesses_sparse_model")) {
      expect_same_fields(
        predict_sparse_gp(model, x_new, include_covariance = TRUE),
        expected$prediction, paste(name, "prediction")
      )
      next
    }

    if (inherits(model, "gaussianprocesses_heteroscedastic_model")) {
      expect_same_fields(
        predict_heteroscedastic_gp(model, c(0.1, 0.5, 0.9),
                                   include_covariance = TRUE),
        expected$prediction, paste(name, "prediction")
      )
      next
    }

    noise <- if (identical(name, "known_noise")) rep(0.03, length(x_new))
    expect_same_fields(
      predict_gp(model, x_new, include_covariance = TRUE,
                 observation_noise_variance = noise),
      expected$prediction, paste(name, "prediction")
    )
    expect_equal(log_marginal_likelihood(model),
                 expected$log_marginal_likelihood,
                 tolerance = cross_platform_tolerance, label = name)
    expect_same_fields(loo_gp(model), expected$loo, paste(name, "loo"))

    if (!is.null(expected$forecast)) {
      expect_same_fields(
        forecast_gp(model, as.Date("2020-01-13") + 0:2),
        expected$forecast, paste(name, "forecast")
      )
    }
  }
})


test_that("schema-4 models, the 1.0 format, give the results they were saved with", {
  frozen <- readRDS(test_path("fixtures", "schema-4-models.rds"))
  x_new <- frozen$x_new

  for (name in names(frozen$models)) {
    model <- frozen$models[[name]]
    expected <- frozen$outputs[[name]]

    expect_identical(model$schema_version, 4L)

    actual <- if (inherits(model, "gaussianprocesses_sparse_model")) {
      list(
        prediction = predict_sparse_gp(model, x_new, include_covariance = TRUE),
        log_marginal_likelihood = model$log_marginal_likelihood
      )
    } else if (inherits(model, "gaussianprocesses_heteroscedastic_model")) {
      list(prediction = predict_heteroscedastic_gp(
        model, c(0.1, 0.5, 0.9), include_covariance = TRUE
      ))
    } else if (inherits(model, "gaussianprocesses_latent_model")) {
      list(prediction = predict_latent_gp(model, x_new),
           log_marginal_likelihood = log_marginal_likelihood(model))
    } else if (!is.null(model$time)) {
      list(forecast = forecast_gp(model, as.Date("2020-01-24") + 0:3,
                                  include_covariance = TRUE),
           log_marginal_likelihood = log_marginal_likelihood(model))
    } else {
      list(prediction = predict_gp(model, x_new, include_covariance = TRUE),
           log_marginal_likelihood = log_marginal_likelihood(model))
    }

    for (part in names(expected)) {
      if (is.list(expected[[part]])) {
        expect_same_fields(actual[[part]], expected[[part]],
                           paste(name, part))
      } else {
        expect_equal(actual[[part]], expected[[part]],
                     tolerance = cross_platform_tolerance,
                     label = paste(name, part))
      }
    }
  }
})


test_that("every function that accepts a model reads schema-0 models", {
  for (case in entry_point_cases()) {
    for (name in names(case$functions)) {
      error <- tryCatch(
        {
          case$functions[[name]](case$model)
          NULL
        },
        error = conditionMessage
      )
      expect_null(error, label = name)
    }
  }
})


test_that("every function that accepts a model refuses unknown schemas", {
  for (case in entry_point_cases()) {
    newer <- case$model
    newer$schema_version <- 7L
    newer$package_version <- "9.0.0"

    for (name in names(case$functions)) {
      condition <- tryCatch(
        case$functions[[name]](newer),
        error = function(condition) condition
      )

      expect_s3_class(condition, "gaussianprocesses_model_format_error")
      expect_identical(condition$found_schema, 7L, info = name)
      expect_identical(
        condition$supported_schemas,
        c(0L, 1L, 2L, 3L, 4L),
        info = name
      )
      expect_identical(condition$package_version, "9.0.0", info = name)
      expect_match(conditionMessage(condition), "schema 7", info = name)
    }
  }
})


test_that("malformed model objects raise format errors", {
  model <- fit_gp(model_x, model_y, rbf_kernel(), noise_variance = 0.02)

  for (version in list("1", 1.5, NA_integer_, 0L, c(1L, 1L))) {
    invalid <- model
    invalid$schema_version <- version

    expect_error(
      predict_gp(invalid, 0.5),
      class = "gaussianprocesses_model_format_error"
    )
  }

  without_factor <- model
  without_factor$cholesky <- NULL
  expect_error(
    log_marginal_likelihood(without_factor),
    "no field\\(s\\): cholesky",
    class = "gaussianprocesses_model_format_error"
  )

  # Class checks are unchanged.
  expect_error(
    predict_gp(list(), 0.5),
    "'model' must be a fitted Gaussian-process model",
    class = "gaussianprocesses_validation_error"
  )
})


# TRUE if `value` contains a matrix, function, or environment at any depth.
contains_heavy_values <- function(value) {
  if (is.matrix(value) || is.function(value) || is.environment(value)) {
    return(TRUE)
  }

  is.list(value) && any(vapply(value, contains_heavy_values, logical(1)))
}


test_that("a model specification is a plain list that refits the model", {
  models <- current_models()

  for (name in c("exact", "known_noise", "optimized", "time")) {
    model <- models[[name]]
    specification <- gaussianprocesses:::.gp_model_specification(model)

    expect_false(contains_heavy_values(specification), info = name)
    expect_identical(specification$schema_version, current_schema)
    expect_identical(
      specification$hyperparameters,
      kernel_parameters(model$kernel, flatten = TRUE)
    )

    refit <- gaussianprocesses:::.fit_gp_from_specification(
      specification,
      model$x,
      model$y
    )
    expect_identical(refit$cholesky, model$cholesky, info = name)
    expect_identical(refit$alpha, model$alpha, info = name)
  }

  specification <- gaussianprocesses:::.gp_model_specification(
    models$optimized
  )
  expect_true(specification$optimization$converged)
  expect_identical(
    gaussianprocesses:::.gp_model_specification(models$time)$time$unit,
    models$time$time$unit
  )

  # Schema-0 models did not record a numerical policy; the defaults refit
  # them exactly because they were fitted with the defaults.
  old <- legacy$models$exact
  specification <- gaussianprocesses:::.gp_model_specification(old)
  expect_null(specification$numerical_policy)
  expect_identical(specification$package_version, NA_character_)
  refit <- gaussianprocesses:::.fit_gp_from_specification(
    specification,
    old$x,
    old$y
  )
  expect_equal(
    refit$cholesky,
    old$cholesky,
    tolerance = cross_platform_tolerance
  )
})
