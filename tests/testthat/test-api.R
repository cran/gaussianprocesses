# The public API: methods for base R generics and the prediction contract
# shared by every model class. inst/notes/api-audit.md
# records the decisions.

api_x <- seq(-2, 2, length.out = 20)
api_y <- sin(2 * api_x) + 0.1 * cos(7 * api_x)
api_x_new <- c(-2.5, -0.3, 0.4, 2.2)

api_models <- function() {
  simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 25)

  list(
    exact = fit_gp(api_x, api_y, rbf_kernel(1.1, 0.6), noise_variance = 0.02),
    optimized = optimize_gp(
      api_x,
      api_y,
      kernel = rbf_kernel(),
      noise_variance = 0.05,
      n_starts = 1
    ),
    sparse = fit_sparse_gp(
      api_x,
      api_y,
      rbf_kernel(1.1, 0.6),
      noise_variance = 0.02,
      n_inducing = 8
    ),
    heteroscedastic = fit_heteroscedastic_gp(
      simulation$x,
      simulation$observed,
      kernel = matern52_kernel(length_scale = 0.2),
      noise_kernel = rbf_kernel(length_scale = 0.4),
      max_iterations = 30
    ),
    time = fit_time_series_gp(
      as.Date("2026-01-01") + 0:23,
      sin((0:23) / 3),
      kernel = rbf_kernel(1, 4),
      noise_variance = 0.01
    )
  )
}


test_that("predict() returns the model-specific prediction", {
  models <- api_models()

  expect_identical(
    predict(models$exact, api_x_new, interval_level = 0.9),
    predict_gp(models$exact, api_x_new, interval_level = 0.9)
  )
  expect_identical(
    predict(models$sparse, api_x_new, include_covariance = TRUE),
    predict_sparse_gp(models$sparse, api_x_new, include_covariance = TRUE)
  )
  expect_identical(
    predict(models$heteroscedastic, c(0.2, 0.8)),
    predict_heteroscedastic_gp(models$heteroscedastic, c(0.2, 0.8))
  )

  future <- as.Date("2026-01-25") + 0:2
  expect_identical(
    predict(models$time, future),
    forecast_gp(models$time, future)
  )

  # Without newdata, predictions are at the training inputs.
  expect_identical(predict(models$exact), predict_gp(models$exact, api_x))
  expect_identical(
    predict(models$sparse)$mean,
    predict_sparse_gp(models$sparse, api_x)$mean
  )
})


test_that("logLik() is the type-II likelihood with estimated-parameter df", {
  models <- api_models()

  exact <- logLik(models$exact)
  expect_s3_class(exact, "logLik")
  expect_equal(
    as.numeric(exact),
    log_marginal_likelihood(models$exact)
  )
  expect_identical(attr(exact, "df"), 0L)
  expect_identical(attr(exact, "nobs"), 20L)

  optimized <- logLik(models$optimized)
  expect_identical(
    attr(optimized, "df"),
    length(models$optimized$optimization$optimized_parameters)
  )
  expect_equal(
    AIC(models$optimized),
    -2 * as.numeric(optimized) + 2 * attr(optimized, "df")
  )
  expect_equal(
    BIC(models$optimized),
    -2 * as.numeric(optimized) + log(20) * attr(optimized, "df")
  )

  sparse <- logLik(models$sparse)
  expect_equal(as.numeric(sparse), models$sparse$log_marginal_likelihood)
  expect_true(is.finite(AIC(models$sparse)))

  expect_error(
    logLik(models$heteroscedastic),
    "plug-in estimate",
    class = "gaussianprocesses_validation_error"
  )
})


test_that("nobs(), fitted(), and residuals() are consistent", {
  models <- api_models()

  for (name in c("exact", "optimized", "sparse", "heteroscedastic")) {
    model <- models[[name]]
    fitted_values <- fitted(model)
    reference <- switch(
      name,
      sparse = predict_sparse_gp(model, model$x)$mean,
      heteroscedastic = predict_heteroscedastic_gp(model, model$x)$mean,
      predict_gp(model, model$x)$mean
    )

    expect_identical(nobs(model), length(model$y), info = name)
    expect_equal(fitted_values, reference, tolerance = 1e-10, info = name)
    expect_identical(residuals(model), model$y - fitted_values, info = name)
  }

  # The fitted-value identity holds when the factorization needed jitter.
  jittered <- fit_gp(c(0, 0, 1, 2), c(1, 1, 0, 0.5), rbf_kernel(), 0)
  expect_gt(jittered$jitter, 0)
  expect_equal(
    fitted(jittered),
    predict_gp(jittered, jittered$x)$mean,
    tolerance = 1e-8
  )
})


test_that("summary() reports the model compactly", {
  models <- api_models()

  for (name in names(models)) {
    summary_object <- summary(models[[name]])
    printed <- utils::capture.output(print(summary_object))

    expect_s3_class(summary_object, "summary.gaussianprocesses_model")
    expect_true(any(grepl("observations:", printed, fixed = TRUE)), info = name)
    expect_true(any(grepl("kernel:", printed, fixed = TRUE)), info = name)
    expect_true(any(grepl("numerical jitter:", printed, fixed = TRUE)), info = name)
  }

  expect_true(
    any(grepl(
      "optimization: converged",
      utils::capture.output(print(summary(models$optimized))),
      fixed = TRUE
    ))
  )
  expect_null(summary(models$heteroscedastic)$log_marginal_likelihood)
})


test_that("every prediction shares the documented fields", {
  models <- api_models()
  contract <- c(
    "x", "mean", "latent_variance", "latent_sd",
    "observation_noise_variance", "observation_variance", "observation_sd",
    "interval_level", "latent_interval", "prediction_interval"
  )
  predictions <- list(
    exact = predict_gp(models$exact, api_x_new, include_covariance = TRUE),
    sparse = predict_sparse_gp(models$sparse, api_x_new, include_covariance = TRUE),
    heteroscedastic = predict_heteroscedastic_gp(
      models$heteroscedastic,
      c(0.1, 0.5, 0.9),
      include_covariance = TRUE
    ),
    forecast = forecast_gp(
      models$time,
      as.Date("2026-01-25") + 0:2,
      include_covariance = TRUE
    )
  )

  for (name in names(predictions)) {
    prediction <- predictions[[name]]
    quantile <- stats::qnorm(0.975)

    expect_true(all(contract %in% names(prediction)), info = name)
    expect_true(
      all(c("latent_covariance", "observation_covariance") %in% names(prediction)),
      info = name
    )

    # The common fields mean the same thing for every model class.
    expect_equal(prediction$latent_sd, sqrt(prediction$latent_variance))
    expect_equal(
      prediction$observation_variance,
      prediction$latent_variance + prediction$observation_noise_variance
    )
    expect_equal(prediction$observation_sd, sqrt(prediction$observation_variance))
    expect_equal(
      unname(prediction$prediction_interval[, "upper"]),
      prediction$mean + quantile * prediction$observation_sd
    )
    expect_equal(diag(prediction$latent_covariance), prediction$latent_variance)
    expect_equal(
      prediction$observation_covariance,
      prediction$latent_covariance +
        diag(prediction$observation_noise_variance, length(prediction$mean))
    )
  }

  # The deprecated copy `noise_variance` was removed in 1.0.0.
  expect_false("noise_variance" %in% names(predictions$heteroscedastic))

  # Without covariances, no prediction carries them.
  marginal <- predict_heteroscedastic_gp(models$heteroscedastic, 0.5)
  expect_null(marginal$latent_covariance)
  expect_null(marginal$observation_covariance)
})
