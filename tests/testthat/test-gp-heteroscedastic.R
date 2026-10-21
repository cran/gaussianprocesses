# Some tests stop the iteration early on purpose, to keep them fast; the
# convergence warning is tested separately.
suppress_convergence_warning <- function(code) {
  withCallingHandlers(
    code,
    gaussianprocesses_convergence_warning = function(warning) {
      invokeRestart("muffleWarning")
    }
  )
}


test_that("fit_gp supports known observation-specific noise", {
  x <- c(0, 1, 2)
  y <- c(1, -0.5, 0.25)
  kernel <- rbf_kernel(variance = 1.4, length_scale = 0.8)
  noise <- c(0.1, 0.2, 0.4)

  model <- fit_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = noise
  )

  covariance <- evaluate_kernel(kernel, x) + diag(noise, 3)

  expect_equal(
    model$alpha,
    solve(covariance, y),
    tolerance = 1e-10
  )
  expect_equal(model$training_noise_variance, noise)
  expect_equal(model$noise_structure, "known_heteroscedastic")
})


test_that("vector-noise models require prediction noise explicitly", {
  model <- fit_gp(
    c(0, 1, 2),
    c(0, 1, 0),
    kernel = rbf_kernel(),
    noise_variance = c(0.1, 0.2, 0.3)
  )

  expect_error(
    predict_gp(model, c(0.5, 1.5)),
    class = "gaussianprocesses_validation_error"
  )

  prediction <- predict_gp(
    model,
    c(0.5, 1.5),
    observation_noise_variance = c(0.15, 0.25)
  )

  expect_equal(
    prediction$observation_variance -
      prediction$latent_variance,
    c(0.15, 0.25),
    tolerance = 1e-12
  )
})


test_that("log squared residual correction has the Gaussian expectation", {
  set.seed(1729)
  sigma2 <- 0.7
  residual <- stats::rnorm(200000, sd = sqrt(sigma2))

  target <- gaussianprocesses:::.log_squared_residual_target(
    residual,
    noise_floor = 1e-12
  )

  expect_equal(
    mean(target),
    log(sigma2),
    tolerance = 0.02
  )
})


test_that("heteroscedastic fit returns separate mean and noise processes", {
  set.seed(42)

  x <- seq(-2, 2, length.out = 40)
  true_noise <- 0.03 + 0.18 * ((x + 2) / 4)^2
  y <- sin(1.5 * x) +
    stats::rnorm(length(x), sd = sqrt(true_noise))

  # Four iterations are deliberately too few to converge, which warns.
  expect_warning(
    model <- fit_heteroscedastic_gp(
      x,
      y,
      kernel = rbf_kernel(
        variance = 1,
        length_scale = 0.7
      ),
      noise_kernel = rbf_kernel(
        variance = 1,
        length_scale = 1
      ),
      initial_noise_variance = 0.08,
      noise_floor = 1e-4,
      damping = 0.6,
      max_iterations = 4
    ),
    class = "gaussianprocesses_convergence_warning"
  )
  expect_false(model$converged)

  expect_s3_class(
    model,
    "gaussianprocesses_heteroscedastic_model"
  )
  expect_s3_class(model$mean_model, "gaussianprocesses_model")
  expect_s3_class(model$noise_model, "gaussianprocesses_model")
  expect_length(model$training_noise_variance, length(x))
  expect_true(all(model$training_noise_variance > 0))
  expect_gte(model$iterations, 1L)
})


test_that("heteroscedastic simulation recovers higher-noise regions", {
  set.seed(1001)

  x <- seq(0, 1, length.out = 60)
  true_noise <- 0.015 + 0.30 * x^2
  latent <- sin(2 * pi * x)
  y <- latent +
    stats::rnorm(length(x), sd = sqrt(true_noise))

  model <- fit_heteroscedastic_gp(
    x,
    y,
    kernel = rbf_kernel(
      variance = 1,
      length_scale = 0.18
    ),
    noise_kernel = rbf_kernel(
      variance = 1,
      length_scale = 0.25
    ),
    initial_noise_variance = 0.08,
    noise_floor = 1e-4,
    damping = 0.6
  )
  expect_true(model$converged)

  prediction <- predict_heteroscedastic_gp(
    model,
    x
  )

  low_noise <- mean(prediction$observation_noise_variance[x <= 0.3])
  high_noise <- mean(prediction$observation_noise_variance[x >= 0.7])

  expect_gt(high_noise, low_noise)
})


test_that("heteroscedastic prediction separates uncertainty components", {
  set.seed(8)

  x <- seq(0, 1, length.out = 30)
  true_noise <- 0.02 + 0.15 * x
  y <- cos(2 * pi * x) +
    stats::rnorm(length(x), sd = sqrt(true_noise))

  model <- suppress_convergence_warning(
    fit_heteroscedastic_gp(
      x,
      y,
      kernel = matern32_kernel(
        variance = 1,
        length_scale = 0.25
      ),
      noise_kernel = rbf_kernel(
        variance = 1,
        length_scale = 0.4
      ),
      initial_noise_variance = 0.06,
      max_iterations = 3
    )
  )

  prediction <- predict_heteroscedastic_gp(
    model,
    c(0.2, 0.8)
  )

  expect_equal(
    prediction$observation_variance,
    prediction$latent_variance +
      prediction$observation_noise_variance,
    tolerance = 1e-12
  )
  expect_true(
    all(prediction$noise_log_variance_uncertainty >= 0)
  )
  expect_true(
    all(
      prediction$noise_variance_interval[, "lower"] <=
        prediction$noise_variance_interval[, "upper"]
    )
  )
})


test_that("heteroscedastic fitting validates approximation controls", {
  x <- c(0, 1, 2)
  y <- c(0, 1, 0)

  expect_error(
    fit_heteroscedastic_gp(
      x,
      y,
      kernel = rbf_kernel(),
      noise_kernel = rbf_kernel(),
      noise_floor = 0
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    fit_heteroscedastic_gp(
      x,
      y,
      kernel = rbf_kernel(),
      noise_kernel = rbf_kernel(),
      damping = 1.5
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    fit_heteroscedastic_gp(
      x[1:2],
      y[1:2],
      kernel = rbf_kernel(),
      noise_kernel = rbf_kernel()
    ),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("heteroscedastic objects have concise print methods", {
  set.seed(99)

  x <- seq(0, 1, length.out = 15)
  y <- sin(2 * pi * x) +
    stats::rnorm(length(x), sd = 0.15)

  model <- suppress_convergence_warning(
    fit_heteroscedastic_gp(
      x,
      y,
      kernel = rbf_kernel(length_scale = 0.2),
      noise_kernel = rbf_kernel(length_scale = 0.3),
      initial_noise_variance = 0.03,
      max_iterations = 2
    )
  )
  prediction <- predict_heteroscedastic_gp(
    model,
    c(0.25, 0.75)
  )

  model_output <- capture.output(print(model))
  prediction_output <- capture.output(print(prediction))

  expect_true(
    any(grepl("heteroscedastic", model_output, fixed = TRUE))
  )
  expect_true(
    any(grepl("noise variance range", model_output, fixed = TRUE))
  )
  expect_true(
    any(grepl("Heteroscedastic", prediction_output, fixed = TRUE))
  )
})


test_that("standardized noise targets are unbiased at the fixed point", {
  # A leave-one-out residual has variance v + sigma^2. With the current noise
  # equal to sigma^2, the standardized target estimates log(sigma^2); the
  # earlier rule estimates log(v + sigma^2). Bounds are five Monte Carlo
  # standard errors: sd(log chi^2_1) = sqrt(trigamma(1/2)) over 2e5 draws.
  set.seed(31)
  n_draws <- 2e5
  noise <- 0.3
  latent_variance <- 0.3
  loo <- list(
    residual = stats::rnorm(n_draws, sd = sqrt(latent_variance + noise)),
    variance = rep(latent_variance + noise, n_draws)
  )
  bound <- 5 * sqrt(trigamma(0.5) / n_draws)
  targets <- function(rule) {
    gaussianprocesses:::.noise_targets(
      loo,
      training_noise = rep(noise, n_draws),
      noise_floor = 1e-300,
      rule = rule
    )
  }

  expect_lt(abs(mean(targets("standardized")) - log(noise)), bound)
  expect_lt(
    abs(mean(targets("residual")) - log(latent_variance + noise)),
    bound
  )
})


test_that("the noise target rule is standardized and can be studied", {
  simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 40)
  fit <- function() {
    fit_heteroscedastic_gp(
      simulation$x,
      simulation$observed,
      kernel = matern52_kernel(length_scale = 0.2),
      noise_kernel = rbf_kernel(length_scale = 0.4)
    )
  }

  model <- fit()
  expect_identical(model$noise_target, "standardized")
  expect_true(model$converged)
  expect_lte(model$iterations, 50L)

  old <- options(gaussianprocesses.heteroscedastic_target = "residual")
  on.exit(options(old), add = TRUE)
  residual_model <- fit()
  expect_identical(residual_model$noise_target, "residual")

  # The standardized rule estimates less noise than the upward-biased one.
  expect_lt(
    mean(log(model$training_noise_variance)),
    mean(log(residual_model$training_noise_variance))
  )

  options(gaussianprocesses.heteroscedastic_target = "neither")
  expect_error(fit(), class = "gaussianprocesses_validation_error")
})


test_that("an unconverged noise estimate warns and says so", {
  simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 30)

  expect_warning(
    model <- fit_heteroscedastic_gp(
      simulation$x,
      simulation$observed,
      kernel = matern52_kernel(length_scale = 0.2),
      noise_kernel = rbf_kernel(length_scale = 0.4),
      max_iterations = 2
    ),
    "did not converge in 2 iterations",
    class = "gaussianprocesses_convergence_warning"
  )
  expect_false(model$converged)
  expect_identical(model$iterations, 2L)
})
