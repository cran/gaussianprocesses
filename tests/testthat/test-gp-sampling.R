# Prior and posterior draws. The moment checks are in helper-sampling.R.

sampling_training_data <- function() {
  x <- c(-2, -1, 0, 1.5)
  list(x = x, y = sin(2 * x))
}


test_that("prior draws reproduce the prior mean and covariance", {
  x <- c(0, 0.5, 2)
  kernel <- rbf_kernel(variance = 1.3, length_scale = 0.8)
  draws <- sample_gp_prior(
    x,
    kernel,
    mean = constant_mean(1.5),
    n_draws = n_statistical_draws,
    seed = 3
  )

  expect_equal(draws$mean, rep(1.5, 3))
  expect_equal(draws$latent_covariance, evaluate_kernel(kernel, x))
  expect_gaussian_moments(
    draws$latent,
    draws$mean,
    draws$latent_covariance,
    "1D prior"
  )
})


test_that("prior draws support multidimensional ARD and composite kernels", {
  x <- cbind(c(0, 0.4, 1, 1.5), c(1, 0, -0.5, 2))
  kernel <- sum_kernel(
    rbf_kernel(variance = 1, length_scale = c(0.6, 2)),
    periodic_kernel(variance = 0.5, length_scale = 1, period = 1.3)
  )
  draws <- sample_gp_prior(
    x,
    kernel,
    n_draws = n_statistical_draws,
    seed = 5
  )

  expect_identical(dim(draws$latent), c(4L, n_statistical_draws))
  expect_gaussian_moments(
    draws$latent,
    rep(0, 4),
    evaluate_kernel(kernel, x),
    "2D composite prior"
  )
})


test_that("one prior draw equals the latent values of simulate_gp_data", {
  x <- seq(0, 3, length.out = 7)
  kernel <- rbf_kernel(variance = 1.3, length_scale = 0.8)

  expect_identical(
    drop(
      sample_gp_prior(x, kernel, mean = constant_mean(0.5), seed = 42)$latent
    ),
    simulate_gp_data(x, kernel, mean = constant_mean(0.5), seed = 42)$latent
  )
})


test_that("latent posterior draws reproduce predict_gp()", {
  data <- sampling_training_data()
  model <- fit_gp(
    data$x,
    data$y,
    kernel = rbf_kernel(variance = 1, length_scale = 0.7),
    noise_variance = 0.01
  )
  x_new <- c(-0.5, 0.25, 3)

  draws <- sample_gp_posterior(
    model,
    x_new,
    n_draws = n_statistical_draws,
    seed = 4
  )
  prediction <- predict_gp(model, x_new, include_covariance = TRUE)

  expect_null(draws$observation)
  expect_equal(draws$mean, prediction$mean)
  expect_equal(draws$latent_covariance, prediction$latent_covariance)
  expect_gaussian_moments(
    draws$latent,
    prediction$mean,
    prediction$latent_covariance,
    "1D latent posterior"
  )
})


test_that("posterior draws support multidimensional ARD composite models", {
  x <- cbind(seq(-1.5, 1.5, length.out = 12), cos(1:12))
  model <- fit_gp(
    x,
    sin(2 * x[, 1]) + 0.3 * x[, 2],
    kernel = product_kernel(
      rbf_kernel(variance = 1.2, length_scale = c(0.8, 1.5)),
      linear_kernel(variance = 0.5)
    ),
    noise_variance = 0.02
  )
  x_new <- cbind(c(-1, 0.2, 2), c(0.5, -0.3, 1))

  draws <- sample_gp_posterior(
    model,
    x_new,
    n_draws = n_statistical_draws,
    seed = 6
  )
  prediction <- predict_gp(model, x_new, include_covariance = TRUE)

  expect_gaussian_moments(
    draws$latent,
    prediction$mean,
    prediction$latent_covariance,
    "2D latent posterior"
  )
})


test_that("observation draws add independent noise to the latent draws", {
  data <- sampling_training_data()
  model <- fit_gp(
    data$x,
    data$y,
    kernel = rbf_kernel(variance = 1, length_scale = 0.7),
    noise_variance = 0.01
  )
  x_new <- c(-0.5, 0.25, 3)

  draws <- sample_gp_posterior(
    model,
    x_new,
    n_draws = n_statistical_draws,
    type = "observation",
    seed = 7
  )
  prediction <- predict_gp(model, x_new, include_covariance = TRUE)

  expect_identical(draws$observation_noise_variance, rep(0.01, 3))
  expect_gaussian_moments(
    draws$observation,
    prediction$mean,
    prediction$observation_covariance,
    "observation posterior"
  )

  # The noise is the difference between paired observation and latent draws.
  noise <- draws$observation - draws$latent
  expect_gaussian_moments(noise, rep(0, 3), diag(0.01, 3), "noise")
  expect_lt(
    max(abs(diag(stats::cor(t(noise), t(draws$latent))))),
    5 / sqrt(n_statistical_draws)
  )
})


test_that("known heteroscedastic models take prediction noise explicitly", {
  data <- sampling_training_data()
  model <- fit_gp(
    data$x,
    data$y,
    kernel = rbf_kernel(variance = 1, length_scale = 0.7),
    noise_variance = c(0.01, 0.02, 0.03, 0.04)
  )
  x_new <- c(0, 1)

  # Latent draws need no noise variance.
  expect_identical(
    dim(sample_gp_posterior(model, x_new, n_draws = 2, seed = 1)$latent),
    c(2L, 2L)
  )
  expect_error(
    sample_gp_posterior(model, x_new, type = "observation"),
    class = "gaussianprocesses_validation_error"
  )

  draws <- sample_gp_posterior(
    model,
    x_new,
    n_draws = n_statistical_draws,
    type = "observation",
    observation_noise_variance = c(0.1, 0.3),
    seed = 8
  )

  expect_identical(draws$observation_noise_variance, c(0.1, 0.3))
  expect_gaussian_moments(
    draws$observation - draws$latent,
    c(0, 0),
    diag(c(0.1, 0.3)),
    "supplied noise"
  )
})


test_that("heteroscedastic draws use the plug-in noise estimate", {
  simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 40)
  model <- fit_heteroscedastic_gp(
    simulation$x,
    simulation$observed,
    kernel = matern52_kernel(length_scale = 0.2),
    noise_kernel = rbf_kernel(length_scale = 0.4),
    max_iterations = 30
  )
  x_new <- c(0.1, 0.5, 0.9)

  draws <- sample_gp_posterior(
    model,
    x_new,
    n_draws = n_statistical_draws,
    type = "observation",
    seed = 9
  )
  prediction <- predict_heteroscedastic_gp(
    model,
    x_new,
    include_covariance = TRUE
  )

  expect_identical(draws$model, "heteroscedastic plug-in")
  expect_equal(draws$observation_noise_variance, prediction$observation_noise_variance)
  expect_gaussian_moments(
    draws$latent,
    prediction$mean,
    prediction$latent_covariance,
    "heteroscedastic latent"
  )
  expect_gaussian_moments(
    draws$observation - draws$latent,
    rep(0, 3),
    diag(prediction$observation_noise_variance),
    "heteroscedastic noise"
  )
})


test_that("draws are reproducible and leave the caller's random state alone", {
  data <- sampling_training_data()
  model <- fit_gp(
    data$x,
    data$y,
    kernel = rbf_kernel(),
    noise_variance = 0.01
  )
  x_new <- c(0, 0.5, 1)

  first <- sample_gp_posterior(model, x_new, n_draws = 3, seed = 11)
  second <- sample_gp_posterior(model, x_new, n_draws = 3, seed = 11)
  other <- sample_gp_posterior(model, x_new, n_draws = 3, seed = 12)

  expect_identical(first$latent, second$latent)
  expect_false(identical(first$latent, other$latent))

  set.seed(99)
  before <- .Random.seed
  sample_gp_posterior(model, x_new, n_draws = 3, seed = 11)
  sample_gp_prior(x_new, rbf_kernel(), n_draws = 3, seed = 11)
  expect_identical(.Random.seed, before)
})


test_that("draw matrices have one row per input and one column per draw", {
  single <- sample_gp_prior(1.5, rbf_kernel(), seed = 1)
  matrix_input <- sample_gp_prior(
    cbind(1:5, 5:1),
    rbf_kernel(length_scale = c(1, 2)),
    n_draws = 4,
    seed = 1
  )

  expect_identical(dim(single$latent), c(1L, 1L))
  expect_identical(dim(matrix_input$latent), c(5L, 4L))
  expect_identical(matrix_input$n_draws, 4L)
})


test_that("nearly singular covariances use the jitter policy", {
  # Closely spaced inputs make the prior covariance numerically singular.
  dense <- sample_gp_prior(
    seq(0, 1, length.out = 200),
    rbf_kernel(length_scale = 2),
    n_draws = 2,
    seed = 1
  )

  expect_gt(dense$jitter, 0)
  expect_true(all(is.finite(dense$latent)))

  # Repeated prediction inputs make the posterior covariance singular.
  data <- sampling_training_data()
  model <- fit_gp(data$x, data$y, rbf_kernel(), noise_variance = 0.01)
  repeated <- sample_gp_posterior(
    model,
    c(0.3, 0.3, 1),
    n_draws = 3,
    seed = 1
  )

  expect_gt(repeated$jitter, 0)
  expect_lt(max(abs(repeated$latent[1, ] - repeated$latent[2, ])), 1e-3)
})


test_that("sampling validates its inputs", {
  data <- sampling_training_data()
  model <- fit_gp(data$x, data$y, rbf_kernel(), noise_variance = 0.01)

  expect_error(
    sample_gp_prior(1:3, rbf_kernel(), n_draws = 0),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    sample_gp_prior(1:3, rbf_kernel(), seed = -1),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    sample_gp_posterior(list(), 0),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    sample_gp_posterior(model, 0, observation_noise_variance = 0.1),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    sample_gp_posterior(model, cbind(0, 1)),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("samples print a concise summary", {
  data <- sampling_training_data()
  model <- fit_gp(data$x, data$y, rbf_kernel(), noise_variance = 0.01)

  output <- capture.output(
    print(
      sample_gp_posterior(model, c(0, 1), n_draws = 2, type = "observation")
    )
  )

  expect_true(any(grepl("posterior (exact)", output, fixed = TRUE)))
  expect_true(any(grepl("noisy observation draws", output, fixed = TRUE)))
})
