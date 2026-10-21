# Generates schema-4-models.rds, used by test-model-objects.R.
#
# Schema 4 is the model format of version 1.0 (#39). The committed fixture
# holds one model of every class and construction route in that format,
# written before 1.0, together with outputs computed when it was written.
# The tests show that every later version reads these models and gives the
# same results. Do not regenerate it once schema 4 has been released: the
# point is that the files stay as they were written. From the repository
# root:
#
#   Rscript -e 'pkgload::load_all(); source("tests/testthat/fixtures/make-schema-4-models.R")'

stopifnot(identical(gaussianprocesses:::.model_schema_version, 4L))

set.seed(4)
x <- seq(-2, 2, length.out = 15)
y <- sin(2 * x) + 0.1 * cos(5 * x) + stats::rnorm(15, sd = 0.05)
x_new <- c(-2.5, -0.3, 0.4, 2.2)
time <- as.Date("2020-01-01") + c(0, 1, 3, 4, 7, 8, 9, 12, 13, 15, 18, 19,
                                   20, 22, 25)
binary <- as.integer(y > 0)
counts <- stats::rpois(15, exp(1 + y))

models <- list(
  estimated_mean = fit_gp(x, y, rbf_kernel(1.1, 0.6), noise_variance = 0.02,
                          mean = linear_mean(estimate_coefficients("reml"))),
  vague_mean = fit_gp(x, y, matern52_kernel(1, 0.8), noise_variance = 0.02,
                      mean = linear_mean(coefficient_prior())),
  derivative = fit_gp(c(x, -1, 0, 1), c(y, 2 * cos(c(-2, 0, 2))),
                      rbf_kernel(1, 0.7), noise_variance = 0.01,
                      derivative = rep(c(0L, 1L), c(15, 3))),
  optimized = optimize_gp(x, y, sum_kernel(rbf_kernel(), white_noise_kernel(0.01)),
                          noise_variance = 0.05, n_starts = 1),
  time_exact = fit_time_series_gp(time, y, matern32_kernel(1, 3),
                                  noise_variance = 0.02),
  time_state_space = fit_time_series_gp(time, y, matern32_kernel(1, 3),
                                        noise_variance = 0.02,
                                        method = "state_space"),
  fitc = fit_sparse_gp(x, y, rbf_kernel(1.1, 0.6), noise_variance = 0.02,
                       n_inducing = 6, method = "fitc"),
  vfe = optimize_sparse_gp(x, y, rbf_kernel(1, 0.6), noise_variance = 0.05,
                           n_inducing = 6, n_starts = 1),
  heteroscedastic = local({
    simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 20)
    fit_heteroscedastic_gp(
      simulation$x, simulation$observed,
      kernel = matern52_kernel(length_scale = 0.2),
      noise_kernel = rbf_kernel(length_scale = 0.4),
      max_iterations = 30
    )
  }),
  bernoulli = fit_latent_gp(x, binary, rbf_kernel(2, 0.8),
                            bernoulli_likelihood("probit"),
                            mean = constant_mean(estimate_coefficients())),
  poisson = fit_latent_gp(x, counts, matern52_kernel(0.5, 1),
                          poisson_likelihood(), exposure = rep(c(1, 2), length.out = 15))
)

outputs <- lapply(models, function(model) {
  if (inherits(model, "gaussianprocesses_sparse_model")) {
    return(list(
      prediction = predict_sparse_gp(model, x_new, include_covariance = TRUE),
      log_marginal_likelihood = model$log_marginal_likelihood
    ))
  }

  if (inherits(model, "gaussianprocesses_heteroscedastic_model")) {
    return(list(
      prediction = predict_heteroscedastic_gp(model, c(0.1, 0.5, 0.9),
                                              include_covariance = TRUE)
    ))
  }

  if (inherits(model, "gaussianprocesses_latent_model")) {
    return(list(
      prediction = predict_latent_gp(model, x_new),
      log_marginal_likelihood = log_marginal_likelihood(model)
    ))
  }

  if (!is.null(model$time)) {
    return(list(
      forecast = forecast_gp(model, as.Date("2020-01-24") + 0:3,
                             include_covariance = TRUE),
      log_marginal_likelihood = log_marginal_likelihood(model)
    ))
  }

  list(
    prediction = predict_gp(model, x_new, include_covariance = TRUE),
    log_marginal_likelihood = log_marginal_likelihood(model)
  )
})

saveRDS(
  list(schema = 4L, x_new = x_new, models = models, outputs = outputs),
  file.path("tests", "testthat", "fixtures", "schema-4-models.rds"),
  version = 2
)
