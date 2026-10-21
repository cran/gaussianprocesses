# Proper scoring rules.

test_that("the CRPS matches numerical integration of its definition", {
  cases <- rbind(
    c(observed = 0.3, mean = 0, sd = 1),
    c(1.7, 0.2, 0.5),
    c(-4, 1, 2.5),
    c(10, 10, 0.01),
    c(0, 3, 0.4)
  )

  for (row in seq_len(nrow(cases))) {
    y <- cases[row, 1]
    mean <- cases[row, 2]
    sd <- cases[row, 3]
    distribution <- function(t) stats::pnorm(t, mean, sd)
    integral <- stats::integrate(
      function(t) distribution(t)^2,
      -Inf,
      y,
      rel.tol = 1e-12
    )$value + stats::integrate(
      function(t) (1 - distribution(t))^2,
      y,
      Inf,
      rel.tol = 1e-12
    )$value

    expect_equal(
      gp_scores(y, mean, sd)$pointwise$crps,
      integral,
      tolerance = 1e-8,
      info = row
    )
  }
})


test_that("scores of calibrated forecasts have their expected values", {
  set.seed(30)
  n <- 20000
  mean <- stats::rnorm(n)
  sd <- stats::runif(n, 0.5, 2)
  observed <- stats::rnorm(n, mean, sd)
  levels <- c(0.5, 0.8, 0.9, 0.95)

  scores <- gp_scores(observed, mean, sd, interval_levels = levels)
  summary <- scores$summary

  # Monte Carlo standard errors of each quantity.
  expect_lt(abs(summary[["pit_mean"]] - 0.5), 5 * sqrt(1 / 12 / n))
  expect_lt(
    abs(summary[["pit_variance"]] - 1 / 12),
    5 * sqrt(1 / 80 - 1 / 144) / sqrt(n)
  )
  expect_true(all(
    abs(scores$calibration$coverage - levels) <
      5 * sqrt(levels * (1 - levels) / n)
  ))

  # E[log score] = -log(2 pi) / 2 - E[log sd] - 1 / 2, with -z^2 / 2 having
  # standard deviation 1 / sqrt(2).
  expected_log_score <- -0.5 * log(2 * pi) - mean(log(sd)) - 0.5
  expect_lt(
    abs(summary[["log_score"]] - expected_log_score),
    5 / sqrt(2 * n)
  )
  expect_identical(nrow(scores$pointwise), as.integer(n))
  expect_identical(scores$calibration$level, levels)

  # Proper scores prefer the calibrated forecast to overconfident and
  # underconfident ones.
  for (factor in c(0.5, 2)) {
    miscalibrated <- gp_scores(observed, mean, factor * sd)$summary
    expect_lt(miscalibrated[["log_score"]], summary[["log_score"]])
    expect_gt(miscalibrated[["crps"]], summary[["crps"]])
  }
})


test_that("interval scores, SMSE, and MSLL follow their definitions", {
  observed <- c(0, 2, -3, 1)
  mean <- c(0.5, 0, 0, 1)
  sd <- c(1, 1, 2, 0.5)

  scores <- gp_scores(
    observed,
    mean,
    sd,
    interval_levels = 0.9,
    reference_mean = 0.2,
    reference_variance = 4
  )

  quantile <- stats::qnorm(0.95)
  lower <- mean - quantile * sd
  upper <- mean + quantile * sd
  expected_interval <- (upper - lower) +
    (2 / 0.1) * pmax(lower - observed, 0) +
    (2 / 0.1) * pmax(observed - upper, 0)

  expect_equal(scores$calibration$interval_score, mean(expected_interval))
  expect_equal(scores$calibration$coverage, 0.75)
  expect_equal(scores$calibration$width, mean(upper - lower))
  expect_equal(
    scores$summary[["smse"]],
    mean((observed - mean)^2) / stats::var(observed)
  )
  expect_equal(
    scores$summary[["msll"]],
    mean(
      -stats::dnorm(observed, mean, sd, log = TRUE) +
        stats::dnorm(observed, 0.2, 2, log = TRUE)
    )
  )
  expect_true(is.na(gp_scores(observed, mean, sd)$summary[["msll"]]))
  expect_true(is.na(gp_scores(1, 0, 1)$summary[["smse"]]))

  for (call in list(
    function() gp_scores(observed, mean[-1], sd),
    function() gp_scores(observed, mean, c(1, 1, 0, 1)),
    function() gp_scores(observed, mean, sd, interval_levels = 1),
    function() gp_scores(observed, mean, sd, interval_levels = c(0.5, 0.5)),
    function() gp_scores(observed, mean, sd, reference_mean = 0),
    function() gp_scores(observed, mean, sd, reference_mean = 0,
                         reference_variance = 0)
  )) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }
})


test_that("models are scored on held-out data and by leave-one-out", {
  set.seed(1)
  x <- seq(0, 6, length.out = 40)
  y <- sin(x) + stats::rnorm(40, sd = 0.1)
  training <- seq(1, 40, by = 2)
  model <- fit_gp(x[training], y[training], rbf_kernel(),
                  noise_variance = 0.01)

  holdout <- gp_holdout_scores(model, x[-training], y[-training])
  prediction <- predict_gp(model, x[-training])
  expect_equal(
    holdout,
    gp_scores(
      y[-training],
      prediction$mean,
      prediction$observation_sd,
      reference_mean = mean(y[training]),
      reference_variance = stats::var(y[training])
    )
  )
  expect_lt(holdout$summary[["msll"]], 0)

  loo <- gp_loo_scores(model, interval_levels = 0.9)
  leave_one_out <- loo_gp(model)
  expect_equal(loo$pointwise$mean, leave_one_out$mean)
  expect_equal(loo$pointwise$sd, leave_one_out$sd)
  expect_identical(loo$calibration$level, 0.9)

  sparse <- fit_sparse_gp(x[training], y[training], rbf_kernel(),
                          noise_variance = 0.01, n_inducing = 8)
  expect_s3_class(
    gp_holdout_scores(sparse, x[-training], y[-training]),
    "gaussianprocesses_scores"
  )

  # Derivative observations do not enter the MSLL reference.
  derivative <- fit_gp(
    c(x[training], 1, 2),
    c(y[training], cos(1), cos(2)),
    rbf_kernel(),
    noise_variance = 0.01,
    derivative = c(rep(0, 20), 1, 1)
  )
  derivative_scores <- gp_holdout_scores(
    derivative,
    x[-training],
    y[-training],
    observation_noise_variance = 0.01
  )
  derivative_prediction <- predict_gp(derivative, x[-training],
                                      observation_noise_variance = 0.01)
  expect_equal(
    derivative_scores$summary[["msll"]],
    gp_scores(
      y[-training],
      derivative_prediction$mean,
      derivative_prediction$observation_sd,
      reference_mean = mean(y[training]),
      reference_variance = stats::var(y[training])
    )$summary[["msll"]]
  )

  output <- utils::capture.output(print(holdout))
  expect_match(output[1], "Predictive scores for 20 observations")
  expect_true(any(grepl("coverage: 50%", output, fixed = TRUE)))
})
