# Generates legacy-models.rds, used by test-model-objects.R.
#
# The committed fixture was generated from commit 4e638aa, before model
# objects had a schema version. It holds one model of every class in that
# unversioned format (schema 0), together with outputs computed by the code
# of that commit, so the tests can show that old models are upgraded and give
# the same results. Do not regenerate it: the point is that it predates the
# schema. From the repository root:
#
#   Rscript -e 'pkgload::load_all(); source("tests/testthat/fixtures/make-legacy-models.R")'

x <- seq(-2, 2, length.out = 12)
y <- sin(2 * x) + 0.1 * cos(5 * x)
x_new <- c(-2.5, -0.3, 0.4, 2.2)

exact <- fit_gp(x, y, rbf_kernel(1.1, 0.6), noise_variance = 0.02)
known_noise <- fit_gp(
  x,
  y,
  matern52_kernel(1, 0.8),
  noise_variance = seq(0.01, 0.06, length.out = 12)
)
optimized <- optimize_gp(
  x,
  y,
  kernel = rbf_kernel(),
  noise_variance = 0.05,
  n_starts = 1
)
time_model <- fit_gp_time(
  1:12,
  y,
  kernel = rbf_kernel(1, 3),
  noise_variance = 0.02
)
sparse <- fit_sparse_gp(
  x,
  y,
  rbf_kernel(1.1, 0.6),
  noise_variance = 0.02,
  n_inducing = 5
)
simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 20)
heteroscedastic <- fit_heteroscedastic_gp(
  simulation$x,
  simulation$observed,
  kernel = matern52_kernel(length_scale = 0.2),
  noise_kernel = rbf_kernel(length_scale = 0.4),
  max_iterations = 30
)

exact_outputs <- function(model, observation_noise_variance = NULL) {
  list(
    prediction = predict_gp(
      model,
      x_new,
      include_covariance = TRUE,
      observation_noise_variance = observation_noise_variance
    ),
    log_marginal_likelihood = log_marginal_likelihood(model),
    gradient = log_marginal_likelihood_gradient(model),
    diagnostics = gp_diagnostics(model),
    samples = sample_gp_posterior(model, x_new, n_draws = 3, seed = 5)$latent
  )
}

legacy <- list(
  x_new = x_new,
  models = list(
    exact = exact,
    known_noise = known_noise,
    optimized = optimized,
    time = time_model,
    sparse = sparse,
    heteroscedastic = heteroscedastic
  ),
  outputs = list(
    exact = exact_outputs(exact),
    known_noise = exact_outputs(known_noise, rep(0.03, length(x_new))),
    optimized = exact_outputs(optimized),
    time = c(
      exact_outputs(time_model),
      list(forecast = forecast_gp(time_model, 13:15))
    ),
    sparse = list(
      prediction = predict_sparse_gp(sparse, x_new, include_covariance = TRUE)
    ),
    heteroscedastic = list(
      prediction = predict_heteroscedastic_gp(
        heteroscedastic,
        c(0.1, 0.5, 0.9),
        include_covariance = TRUE
      ),
      samples = sample_gp_posterior(
        heteroscedastic,
        c(0.1, 0.5, 0.9),
        n_draws = 3,
        type = "observation",
        seed = 5
      )$observation
    )
  )
)

saveRDS(
  legacy,
  file.path("tests", "testthat", "fixtures", "legacy-models.rds"),
  version = 2
)
