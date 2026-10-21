test_that("mathematical vignette examples reproduce core identities", {
  x <- c(-1, 0, 1)
  y <- c(1, -0.25, 0.5)
  kernel <- rbf_kernel(
    variance = 1.2,
    length_scale = 0.7
  )
  noise <- 0.1

  model <- fit_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = noise
  )

  covariance <- evaluate_kernel(kernel, x) +
    diag(noise, length(x))

  expect_equal(
    model$alpha,
    solve(covariance, y),
    tolerance = 1e-10
  )

  expected_lml <-
    -0.5 * drop(crossprod(y, solve(covariance, y))) -
    0.5 * as.numeric(
      determinant(covariance, logarithm = TRUE)$modulus
    ) -
    0.5 * length(y) * log(2 * pi)

  expect_equal(
    log_marginal_likelihood(model),
    expected_lml,
    tolerance = 1e-10
  )
})


test_that("kernel composition example matches explicit covariance algebra", {
  time <- 0:6

  trend <- linear_kernel(variance = 0.05)
  seasonal <- periodic_kernel(
    variance = 1,
    length_scale = 0.8,
    period = 12
  )
  local <- matern32_kernel(
    variance = 0.3,
    length_scale = 2
  )

  combined <- sum_kernel(
    trend,
    product_kernel(seasonal, local)
  )

  expect_equal(
    evaluate_kernel(combined, time),
    evaluate_kernel(trend, time) +
      evaluate_kernel(seasonal, time) *
      evaluate_kernel(local, time),
    tolerance = 1e-12
  )
})


test_that("time-series vignette workflow distinguishes interpolation and forecast", {
  time <- 0:12
  response <- sin(2 * pi * time / 6)

  model <- fit_time_series_gp(
    time,
    response,
    kernel = time_series_kernel(
      period = 6,
      trend_variance = 0.01,
      periodic_variance = 1,
      periodic_length_scale = 0.8,
      local_variance = 0.1,
      local_length_scale = 1.5
    ),
    noise_variance = 0.05
  )

  prediction <- forecast_gp(
    model,
    c(4.5, 13)
  )

  expect_equal(
    prediction$region,
    c("interpolation", "forecast")
  )
})


test_that("ARD vignette interpretation maps to scaled covariance", {
  x <- rbind(
    c(0, 0),
    c(1, 0),
    c(0, 1)
  )

  covariance <- evaluate_kernel(
    rbf_kernel(
      variance = 1,
      length_scale = c(0.4, 2)
    ),
    x
  )

  expect_lt(
    covariance[1, 2],
    covariance[1, 3]
  )
})
