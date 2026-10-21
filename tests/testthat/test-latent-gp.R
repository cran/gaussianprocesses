# Prediction, sampling, optimization, and methods of latent models.

latent_data <- function() {
  set.seed(33)
  x <- seq(0, 6, length.out = 30)
  f <- 1.5 * sin(x)

  list(
    x = x,
    binary = as.numeric(f + stats::rnorm(30, sd = 0.5) > 0),
    counts = stats::rpois(30, 2 * exp(f)),
    exposure = rep(c(1, 2), length.out = 30)
  )
}


test_that("latent predictions follow Algorithm 3.2 and the likelihood", {
  data <- latent_data()
  kernel <- rbf_kernel(variance = 2, length_scale = 1)
  model <- fit_latent_gp(data$x, data$binary, kernel, bernoulli_likelihood())
  x_new <- c(-1, 0.5, 2.25, 5.9, 12)

  # Direct computation with solve(): mean k_*' d log p / df and covariance
  # K_** - K_*' (K + W^{-1})^{-1} K_*. The model predicts with
  # a = K^{-1} (f_hat - m), which equals d log p / df at the mode up to the
  # Newton convergence tolerance, hence the tolerance of the mean.
  covariance <- evaluate_kernel(kernel, data$x)
  cross <- evaluate_kernel(kernel, data$x, x_new)
  gradient <- evaluate_likelihood(bernoulli_likelihood(), data$binary,
                                  model$mode)$first
  expected_mean <- drop(crossprod(cross, gradient))
  expected_covariance <- evaluate_kernel(kernel, x_new) -
    crossprod(cross, solve(covariance + diag(1 / model$curvature), cross))

  prediction <- predict_latent_gp(model, x_new, include_covariance = TRUE)
  expect_s3_class(prediction, "gaussianprocesses_latent_prediction")
  expect_equal(prediction$mean, expected_mean, tolerance = 1e-6)
  expect_equal(prediction$latent_covariance, expected_covariance,
               tolerance = 1e-8)
  expect_equal(prediction$latent_variance, diag(expected_covariance),
               tolerance = 1e-8)
  expect_equal(
    prediction$probability,
    likelihood_predictive(bernoulli_likelihood(), prediction$mean,
                          prediction$latent_variance)$probability
  )
  expect_equal(prediction$response_mean, prediction$probability)
  expect_equal(prediction$latent_interval[, "upper"],
               prediction$mean + stats::qnorm(0.975) * prediction$latent_sd)
  expect_null(prediction$exposure)

  # Far from the data the probability returns to 1/2, while plogis() of the
  # latent mean does not account for the latent variance.
  expect_lt(abs(prediction$probability[5] - 0.5), 1e-3)

  # Marginal predictions, in any block size, agree with the full covariance.
  for (block_size in list(NULL, 1L, 2L)) {
    marginal <- predict_latent_gp(model, x_new, block_size = block_size)
    expect_equal(marginal$mean, prediction$mean, tolerance = 1e-12)
    expect_equal(marginal$latent_variance, prediction$latent_variance,
                 tolerance = 1e-12)
    expect_null(marginal$latent_covariance)
  }

  # Poisson: the predictive mean is that of a lognormal rate with exposure.
  counts <- fit_latent_gp(data$x, data$counts, kernel, poisson_likelihood(),
                          exposure = data$exposure)
  count_prediction <- predict_latent_gp(counts, x_new, exposure = 3)
  expect_identical(count_prediction$exposure, rep(3, 5))
  expect_equal(
    count_prediction$response_mean,
    3 * exp(count_prediction$mean + count_prediction$latent_variance / 2)
  )
  expect_equal(
    predict_latent_gp(counts, x_new)$response_mean,
    count_prediction$response_mean / 3
  )
  expect_null(count_prediction$probability)

  for (call in list(
    function() predict_latent_gp(model, cbind(x_new, x_new)),
    function() predict_latent_gp(model, x_new, interval_level = 1),
    function() predict_latent_gp(model, x_new, include_covariance = NA),
    function() predict_latent_gp(model, x_new, exposure = 2),
    function() predict_latent_gp(counts, x_new, exposure = c(1, 2)),
    function() predict_latent_gp(fit_gp(data$x, data$binary, kernel), x_new)
  )) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }

  output <- utils::capture.output(print(prediction))
  expect_match(output[1], "Latent Gaussian-process prediction")
  expect_true(any(grepl("probability range", output)))
  expect_true(any(grepl(
    "response mean range",
    utils::capture.output(print(count_prediction))
  )))
})


test_that("posterior draws pass moment tests against the Gaussian approximation", {
  data <- latent_data()
  kernel <- rbf_kernel(variance = 2, length_scale = 1)
  model <- fit_latent_gp(data$x, data$binary, kernel, bernoulli_likelihood())
  x_new <- c(-0.5, 1.2, 1.4, 4)
  prediction <- predict_latent_gp(model, x_new, include_covariance = TRUE)

  draws <- sample_gp_posterior(model, x_new, n_draws = n_statistical_draws,
                               type = "observation", seed = 11)
  expect_s3_class(draws, "gaussianprocesses_samples")
  expect_identical(draws$model, "Laplace approximation")
  expect_equal(draws$mean, prediction$mean)
  expect_equal(draws$latent_covariance, prediction$latent_covariance)
  expect_gaussian_moments(draws$latent, prediction$mean,
                          prediction$latent_covariance, "Laplace latent")

  # Bernoulli draws have mean Pr(y = 1), with standard error
  # sqrt(p (1 - p) / N).
  probability <- prediction$probability
  expect_true(all(draws$observation %in% c(0, 1)))
  expect_lt(
    max(abs(rowMeans(draws$observation) - probability) /
          sqrt(probability * (1 - probability) / n_statistical_draws)),
    5
  )
  # Observations are drawn from the latent draws they are paired with.
  expect_lt(
    max(abs(rowMeans(draws$observation - stats::plogis(draws$latent))) /
          sqrt(0.25 / n_statistical_draws)),
    5
  )

  counts <- fit_latent_gp(data$x, data$counts, kernel, poisson_likelihood())
  count_draws <- sample_gp_posterior(counts, x_new,
                                     n_draws = n_statistical_draws,
                                     type = "observation", exposure = 2,
                                     seed = 12)
  count_prediction <- predict_latent_gp(counts, x_new, exposure = 2)
  expect_true(all(count_draws$observation == round(count_draws$observation)))
  expect_lt(
    max(abs(rowMeans(count_draws$observation) -
              count_prediction$response_mean) /
          sqrt(count_prediction$response_variance / n_statistical_draws)),
    5
  )

  # Gaussian observations add the likelihood variance; probit observations
  # are binary.
  gaussian <- fit_latent_gp(data$x, sin(data$x), kernel,
                            gaussian_likelihood(0.04))
  gaussian_draws <- sample_gp_posterior(gaussian, x_new,
                                        n_draws = n_statistical_draws,
                                        type = "observation", seed = 13)
  expect_gaussian_moments(
    gaussian_draws$observation - gaussian_draws$latent,
    rep(0, 4),
    diag(0.04, 4),
    "Gaussian likelihood noise"
  )
  probit <- fit_latent_gp(data$x, data$binary, kernel,
                          bernoulli_likelihood("probit"))
  probit_draws <- sample_gp_posterior(probit, x_new, n_draws = 10,
                                      type = "observation", seed = 14)
  expect_true(all(probit_draws$observation %in% c(0, 1)))

  # Latent draws only, reproducible, and the caller's random state is kept.
  set.seed(99)
  before <- .Random.seed
  latent_only <- sample_gp_posterior(model, x_new, n_draws = 3, seed = 4)
  expect_identical(.Random.seed, before)
  expect_null(latent_only$observation)
  expect_identical(
    sample_gp_posterior(model, x_new, n_draws = 3, seed = 4)$latent,
    latent_only$latent
  )
  expect_match(
    utils::capture.output(print(draws))[4],
    "latent function draws and observation draws"
  )
})


test_that("latent hyperparameters are optimized with the shared optimizer", {
  data <- latent_data()
  x <- data$x
  y <- sin(x) + stats::rnorm(30, sd = 0.2)

  # With a Gaussian likelihood the optimum is that of optimize_gp().
  latent <- optimize_latent_gp(x, y, rbf_kernel(), gaussian_likelihood(0.1),
                               n_starts = 1)
  exact <- optimize_gp(x, y, rbf_kernel(), noise_variance = 0.1, n_starts = 1)
  expect_equal(
    unname(latent$optimization$optimized_parameters),
    unname(exact$optimization$optimized_parameters),
    tolerance = 1e-5
  )
  expect_equal(log_marginal_likelihood(latent), log_marginal_likelihood(exact),
               tolerance = 1e-8)
  expect_identical(
    names(latent$optimization$optimized_parameters),
    c("variance", "length_scale", "likelihood.variance")
  )
  expect_equal(latent$likelihood$parameters$variance,
               exact$noise_variance, tolerance = 1e-5)

  # Binary responses: the analytical gradient reaches the optimum that
  # finite-difference gradients reach.
  analytical <- optimize_latent_gp(x, data$binary, rbf_kernel(),
                                   bernoulli_likelihood(), n_starts = 2)
  numerical <- optimize_latent_gp(x, data$binary, rbf_kernel(),
                                  bernoulli_likelihood(), n_starts = 1,
                                  gradient = "numerical")
  expect_equal(log_marginal_likelihood(numerical),
               log_marginal_likelihood(analytical), tolerance = 1e-6)

  record <- analytical$optimization
  expect_true(record$converged)
  expect_identical(record$approximation, "laplace")
  expect_false(record$optimize_likelihood)
  expect_identical(record$n_starts, 2L)
  expect_identical(nrow(record$starts), 2L)
  expect_length(record$start_results, 2L)
  expect_s3_class(record$optima, "data.frame")
  expect_gt(record$model_fits, 0L)
  expect_identical(record$log_marginal_likelihood,
                   log_marginal_likelihood(analytical))
  expect_named(record$optimized_parameters, c("variance", "length_scale"))

  # The returned model is refitted from the prior mean at the optimum.
  refit <- fit_latent_gp(x, data$binary, analytical$kernel,
                         bernoulli_likelihood())
  expect_identical(refit$mode, analytical$mode)
  expect_identical(refit$log_marginal_likelihood,
                   analytical$log_marginal_likelihood)

  # A likelihood parameter can be held fixed.
  fixed <- optimize_latent_gp(x, y, rbf_kernel(), gaussian_likelihood(0.1),
                              optimize_likelihood = FALSE, n_starts = 1)
  expect_identical(fixed$likelihood$parameters$variance, 0.1)
  expect_named(fixed$optimization$optimized_parameters,
               c("variance", "length_scale"))

  for (call in list(
    function() optimize_latent_gp(x, y, rbf_kernel(), gaussian_likelihood(),
                                  optimize_likelihood = NA),
    function() optimize_latent_gp(x, y, rbf_kernel(), gaussian_likelihood(),
                                  n_starts = 0),
    function() optimize_latent_gp(x, y, rbf_kernel(), gaussian_likelihood(),
                                  control = 1),
    function() optimize_latent_gp(x, y, rbf_kernel(), gaussian_likelihood(),
                                  lower = 2, upper = 1)
  )) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }
})


test_that("latent models support the base R methods and persistence", {
  data <- latent_data()
  model <- optimize_latent_gp(data$x, data$counts, rbf_kernel(),
                              poisson_likelihood(), exposure = data$exposure,
                              n_starts = 1)

  expect_identical(nobs(model), 30L)
  expect_identical(fitted(model), model$mode)
  expect_identical(
    predict(model, c(1, 2)),
    predict_latent_gp(model, c(1, 2))
  )
  expect_identical(predict(model)$mean, predict_latent_gp(model, data$x)$mean)

  log_lik <- logLik(model)
  expect_s3_class(log_lik, "logLik")
  expect_identical(as.numeric(log_lik), model$log_marginal_likelihood)
  expect_identical(attr(log_lik, "df"), 2L)
  expect_identical(
    attr(logLik(fit_latent_gp(data$x, data$counts, rbf_kernel(),
                              poisson_likelihood())), "df"),
    0L
  )

  output <- utils::capture.output(print(model))
  expect_match(output[1], "Latent Gaussian-process model")
  expect_true(any(grepl("likelihood: PoissonLikelihood\\(link=log\\)",
                        output)))
  expect_true(any(grepl("optimization: converged", output)))

  summary_output <- utils::capture.output(print(summary(model)))
  expect_match(summary_output[1], "Latent Gaussian-process \\(Laplace\\) model")
  expect_false(any(grepl("noise:", summary_output)))
  expect_true(any(grepl("newton converged: TRUE", summary_output)))

  # A saved model gives identical results.
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(model, path)
  restored <- readRDS(path)
  expect_identical(predict_latent_gp(restored, c(0.5, 7)),
                   predict_latent_gp(model, c(0.5, 7)))
  expect_identical(log_marginal_likelihood_gradient(restored),
                   log_marginal_likelihood_gradient(model))
  expect_identical(
    sample_gp_posterior(restored, 1, n_draws = 2, seed = 1),
    sample_gp_posterior(model, 1, n_draws = 2, seed = 1)
  )
  expect_identical(model$schema_version, .model_schema_version)

  broken <- model
  broken$b_cholesky <- NULL
  expect_error(predict_latent_gp(broken, 1),
               class = "gaussianprocesses_model_format_error")
  expect_error(log_marginal_likelihood(structure(list(), class =
    "gaussianprocesses_latent_model")),
    class = "gaussianprocesses_model_format_error")
})
