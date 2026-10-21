# Generates release-0.1.0-models.rds, used by test-model-objects.R.
#
# The committed fixture holds one model of every class saved by the released
# gaussianprocesses 0.1.0 (tag v0.1.0), the only release before the model
# format was versioned, together with outputs computed by that release. The
# tests show that models saved by 0.1.0 are read by every later version and
# give the same results (#39). Do not regenerate it with a later version: the
# point is that it was written by 0.1.0. To reproduce it, install 0.1.0 from
# the tag into a separate library and run, from the repository root,
#
#   git worktree add ../gp-0.1.0 v0.1.0
#   R CMD INSTALL -l <library> ../gp-0.1.0
#   R_LIBS=<library> Rscript tests/testthat/fixtures/make-release-models.R

library(gaussianprocesses)
stopifnot(identical(as.character(utils::packageVersion("gaussianprocesses")),
                    "0.1.0"))

x <- seq(-2, 2, length.out = 12)
y <- sin(2 * x) + 0.1 * cos(5 * x)
x_new <- c(-2.5, -0.3, 0.4, 2.2)

exact <- fit_gp(x, y, rbf_kernel(1.1, 0.6), noise_variance = 0.02)
constant <- fit_gp(x, y, matern32_kernel(0.8, 0.7), noise_variance = 0.03,
                   mean = constant_mean(0.4))
known_noise <- fit_gp(
  x,
  y,
  matern52_kernel(1, 0.8),
  noise_variance = seq(0.01, 0.06, length.out = 12)
)
optimized <- optimize_gp(
  x,
  y,
  kernel = sum_kernel(rbf_kernel(), periodic_kernel(0.3, 1, 2)),
  noise_variance = 0.05,
  n_starts = 1
)
time_model <- fit_gp_time(
  as.Date("2020-01-01") + 0:11,
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
    loo = loo_gp(model)
  )
}

release <- list(
  version = "0.1.0",
  x_new = x_new,
  models = list(
    exact = exact,
    constant = constant,
    known_noise = known_noise,
    optimized = optimized,
    time = time_model,
    sparse = sparse,
    heteroscedastic = heteroscedastic
  ),
  outputs = list(
    exact = exact_outputs(exact),
    constant = exact_outputs(constant),
    known_noise = exact_outputs(known_noise, rep(0.03, length(x_new))),
    optimized = exact_outputs(optimized),
    time = c(
      exact_outputs(time_model),
      list(forecast = forecast_gp(time_model, as.Date("2020-01-13") + 0:2))
    ),
    sparse = list(
      prediction = predict_sparse_gp(sparse, x_new, include_covariance = TRUE)
    ),
    heteroscedastic = list(
      prediction = predict_heteroscedastic_gp(
        heteroscedastic,
        c(0.1, 0.5, 0.9),
        include_covariance = TRUE
      )
    )
  )
)

saveRDS(
  release,
  file.path("tests", "testthat", "fixtures", "release-0.1.0-models.rds"),
  version = 2
)
