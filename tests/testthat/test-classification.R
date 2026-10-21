# Binary classification: scores, probabilities, labels, and separable data.

classification_data <- function(n = 60, seed = 34) {
  set.seed(seed)
  x <- sort(stats::runif(n, 0, 6))

  list(x = x, y = stats::rbinom(n, 1, stats::plogis(2 * sin(x))))
}


test_that("classification scores follow their definitions", {
  observed <- c(1, 0, 1, 1, 0, 0, 1, 0)
  probability <- c(0.9, 0.2, 0.6, 0.6, 0.6, 0.05, 0.35, 0.5)

  scores <- gp_classification_scores(observed, probability, n_bins = 4)
  expect_s3_class(scores, "gaussianprocesses_classification_scores")

  log_density <- ifelse(observed == 1, log(probability), log(1 - probability))
  expect_equal(scores$summary[["log_loss"]], -mean(log_density))
  expect_equal(scores$summary[["brier_score"]],
               mean((probability - observed)^2))
  expect_identical(scores$summary[["n"]], 8)
  expect_identical(scores$summary[["base_rate"]], 0.5)
  expect_equal(scores$summary[["mean_probability"]], mean(probability))

  # The AUC counts every event-non-event pair, ties as one half.
  pairs <- outer(probability[observed == 1], probability[observed == 0], "-")
  expect_equal(scores$summary[["auc"]],
               mean((pairs > 0) + 0.5 * (pairs == 0)))

  # Bins [0, 0.25), [0.25, 0.5), [0.5, 0.75), [0.75, 1].
  expect_identical(scores$reliability$bin, 1:4)
  expect_identical(scores$reliability$n, c(2L, 1L, 4L, 1L))
  expect_equal(scores$reliability$mean_probability,
               c(0.125, 0.35, 0.575, 0.9))
  expect_equal(scores$reliability$observed_frequency, c(0, 1, 0.5, 1))
  expect_equal(
    scores$reliability$standard_error,
    sqrt(c(0.125 * 0.875 / 2, 0.35 * 0.65, 0.575 * 0.425 / 4, 0.9 * 0.1))
  )
  expect_equal(
    scores$summary[["ece"]],
    (2 * 0.125 + 1 * 0.65 + 4 * 0.075 + 1 * 0.1) / 8
  )
  expect_equal(scores$pointwise$log_density, log_density)

  # Supplied log densities are used for the log loss.
  supplied <- gp_classification_scores(observed, probability,
                                       log_density = rep(-1, 8))
  expect_identical(supplied$summary[["log_loss"]], 1)

  # Labels are coded as in bernoulli_likelihood().
  labels <- factor(ifelse(observed == 1, "yes", "no"), levels = c("no", "yes"))
  expect_identical(
    gp_classification_scores(labels, probability, n_bins = 4)$summary,
    scores$summary
  )
  expect_identical(
    gp_classification_scores(observed == 1, probability, n_bins = 4)$summary,
    scores$summary
  )
  expect_true(is.na(gp_classification_scores(c(1, 1), c(0.2, 0.7))$
                      summary[["auc"]]))

  output <- utils::capture.output(print(scores))
  expect_match(output[1], "Classification scores for 8 responses")
  expect_true(any(grepl("ROC AUC", output)))

  for (call in list(
    function() gp_classification_scores(observed, probability[-1]),
    function() gp_classification_scores(observed, probability + 0.2),
    function() gp_classification_scores(observed + 1, probability),
    function() gp_classification_scores(observed, probability, n_bins = 0),
    function() gp_classification_scores(observed, probability,
                                        log_density = rep(0.1, 8)),
    function() gp_classification_scores(observed, probability,
                                        log_density = rep(-1, 7))
  )) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }
})


test_that("calibrated forecasts score as expected", {
  set.seed(35)
  n <- 20000
  probability <- stats::runif(n)
  observed <- stats::rbinom(n, 1, probability)
  scores <- gp_classification_scores(observed, probability)

  # Every bin's observed frequency is within 5 standard errors of its mean
  # probability.
  reliability <- scores$reliability
  expect_lt(
    max(abs(reliability$observed_frequency - reliability$mean_probability) /
          reliability$standard_error),
    5
  )

  # E[(p - y)^2] = E[p (1 - p)] = 1/6, with standard deviation below 0.3 per
  # forecast; for p ~ U(0, 1) the AUC is 5/6.
  expect_lt(abs(scores$summary[["brier_score"]] - 1 / 6), 5 * 0.3 / sqrt(n))
  expect_lt(abs(scores$summary[["auc"]] - 5 / 6), 0.02)
  expect_lt(scores$summary[["ece"]], 0.02)

  # Proper scores prefer the calibrated forecasts.
  for (distorted in list(probability^2, sqrt(probability))) {
    worse <- gp_classification_scores(observed, distorted)$summary
    expect_gt(worse[["brier_score"]], scores$summary[["brier_score"]])
    expect_gt(worse[["log_loss"]], scores$summary[["log_loss"]])
    expect_gt(worse[["ece"]], scores$summary[["ece"]])
  }
})


test_that("held-out scores of latent models use the predictive distribution", {
  data <- classification_data()
  training <- seq(1, 60, by = 2)
  model <- fit_latent_gp(data$x[training], data$y[training],
                         rbf_kernel(4, 1), bernoulli_likelihood())
  held_out <- data$x[-training]
  prediction <- predict_latent_gp(model, held_out)

  scores <- gp_holdout_scores(model, held_out, data$y[-training], n_bins = 5)
  expect_equal(
    scores,
    gp_classification_scores(
      data$y[-training],
      prediction$probability,
      n_bins = 5,
      log_density = likelihood_predictive(
        bernoulli_likelihood(),
        prediction$mean,
        prediction$latent_variance,
        observed = data$y[-training]
      )$log_density
    )
  )
  expect_gt(scores$summary[["auc"]], 0.6)

  # Gaussian likelihoods give the scores of the Gaussian predictive.
  continuous <- sin(data$x)
  gaussian <- fit_latent_gp(data$x[training], continuous[training],
                            rbf_kernel(), gaussian_likelihood(0.01))
  exact <- fit_gp(data$x[training], continuous[training], rbf_kernel(),
                  noise_variance = 0.01)
  expect_equal(
    gp_holdout_scores(gaussian, held_out, continuous[-training]),
    gp_holdout_scores(exact, held_out, continuous[-training]),
    tolerance = 1e-8
  )

  for (call in list(
    function() gp_holdout_scores(model, held_out, data$y[1:3]),
    function() gp_holdout_scores(model, held_out, rep(2, 30))
  )) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }
})


test_that("class probabilities are exact for probit and converge for logit", {
  data <- classification_data(n = 100)
  probit <- fit_latent_gp(data$x, data$y, rbf_kernel(4, 1),
                          bernoulli_likelihood("probit"))
  prediction <- predict_latent_gp(probit, seq(-1, 7, by = 0.25))
  expect_equal(
    prediction$probability,
    .gaussian_expectation(function(f, index) stats::pnorm(f),
                          prediction$mean, prediction$latent_variance, 40),
    tolerance = 1e-10
  )

  # Logistic probabilities: 20 and 40 nodes agree to 1e-10 for latent
  # variances up to 1, and the default 40 nodes agree with 80 up to 3.
  mean <- rep(seq(-6, 6, by = 0.25), 3)
  logistic <- function(variance, order) {
    likelihood_predictive(bernoulli_likelihood(), mean,
                          rep(variance, length(mean)),
                          quadrature_order = order)$probability
  }

  for (variance in c(0.25, 0.5, 1)) {
    expect_lt(max(abs(logistic(variance, 20) - logistic(variance, 40))),
              1e-10, label = paste("variance", variance))
  }

  for (variance in c(1, 2, 3)) {
    expect_lt(max(abs(logistic(variance, 40) - logistic(variance, 80))),
              1e-10, label = paste("variance", variance))
  }

  # The documented errors of the plug-in and MacKay probabilities.
  grid <- seq(-10, 10, by = 0.05)
  plug_in_error <- function(variance) {
    max(abs(stats::plogis(grid) -
              logistic_on_grid(grid, variance)))
  }
  logistic_on_grid <- function(grid, variance) {
    likelihood_predictive(bernoulli_likelihood(), grid,
                          rep(variance, length(grid)),
                          quadrature_order = 80)$probability
  }
  mackay_error <- function(variance) {
    max(abs(stats::plogis(grid / sqrt(1 + pi * variance / 8)) -
              logistic_on_grid(grid, variance)))
  }

  expect_identical(round(plug_in_error(1), 2), 0.04)
  expect_identical(round(plug_in_error(4), 2), 0.11)
  expect_identical(round(mackay_error(1), 3), 0.005)
  expect_identical(round(mackay_error(4), 2), 0.01)
})


test_that("probability intervals and draws describe the class probability", {
  data <- classification_data()
  model <- fit_latent_gp(data$x, data$y, rbf_kernel(4, 1),
                         bernoulli_likelihood())
  x_new <- c(0.5, 2, 4.5, 9)
  prediction <- predict_latent_gp(model, x_new, interval_level = 0.9)

  expect_identical(prediction$probability_interval,
                   stats::plogis(prediction$latent_interval))
  expect_true(all(prediction$probability > prediction$probability_interval[,
                                                                          1]))
  expect_true(all(prediction$probability < prediction$probability_interval[,
                                                                          2]))

  draws <- sample_gp_posterior(model, x_new, n_draws = n_statistical_draws,
                               seed = 21)
  expect_identical(draws$probability, stats::plogis(draws$latent))

  # The draws average to the predictive probability, and their 5% and 95%
  # quantiles match the interval within Monte Carlo error.
  expect_lt(
    max(abs(rowMeans(draws$probability) - prediction$probability) /
          sqrt(apply(draws$probability, 1, stats::var) / n_statistical_draws)),
    5
  )
  coverage_lower <- rowMeans(
    draws$probability < prediction$probability_interval[, 1]
  )
  coverage_upper <- rowMeans(
    draws$probability > prediction$probability_interval[, 2]
  )
  standard_error <- sqrt(0.05 * 0.95 / n_statistical_draws)
  expect_lt(max(abs(c(coverage_lower, coverage_upper) - 0.05)),
            5 * standard_error)

  # Probit models transform with the normal distribution function.
  probit <- fit_latent_gp(data$x, data$y, rbf_kernel(4, 1),
                          bernoulli_likelihood("probit"))
  probit_draws <- sample_gp_posterior(probit, x_new, n_draws = 3, seed = 1)
  expect_identical(probit_draws$probability, stats::pnorm(probit_draws$latent))
  expect_null(predict_latent_gp(
    fit_latent_gp(data$x, rpois(60, 2), rbf_kernel(), poisson_likelihood()),
    x_new
  )$probability_interval)
})


test_that("labels map to the event coded 1", {
  data <- classification_data()
  labels <- factor(ifelse(data$y == 1, "malignant", "benign"),
                   levels = c("benign", "malignant"))
  kernel <- rbf_kernel(4, 1)
  x_new <- c(1, 3, 5)

  numeric_model <- fit_latent_gp(data$x, data$y, kernel,
                                 bernoulli_likelihood())
  factor_model <- fit_latent_gp(data$x, labels, kernel,
                                bernoulli_likelihood())
  logical_model <- fit_latent_gp(data$x, data$y == 1, kernel,
                                 bernoulli_likelihood())
  expect_identical(factor_model$response_levels, c("benign", "malignant"))
  expect_identical(predict_latent_gp(factor_model, x_new),
                   predict_latent_gp(numeric_model, x_new))
  expect_identical(predict_latent_gp(logical_model, x_new),
                   predict_latent_gp(numeric_model, x_new))

  # Reversing the levels makes the other class the event: under a zero mean
  # the latent function changes sign and the probabilities complement.
  reversed <- fit_latent_gp(
    data$x,
    factor(labels, levels = c("malignant", "benign")),
    kernel,
    bernoulli_likelihood()
  )
  expect_equal(reversed$mode, -numeric_model$mode, tolerance = 1e-12)
  expect_equal(predict_latent_gp(reversed, x_new)$probability,
               1 - predict_latent_gp(numeric_model, x_new)$probability,
               tolerance = 1e-12)

  for (y in list(
    factor(c(rep("a", 20), rep("b", 20), rep("c", 20))),
    c(data$y[-1], NA),
    c(data$y[-1], 2),
    as.character(data$y)
  )) {
    expect_error(fit_latent_gp(data$x, y, kernel, bernoulli_likelihood()),
                 class = "gaussianprocesses_validation_error")
  }
})


test_that("separable classes push the signal variance to its bound", {
  x <- 1:20
  y <- as.numeric(x > 10)

  expect_warning(
    bounded <- optimize_latent_gp(x, y, rbf_kernel(), bernoulli_likelihood(),
                                  n_starts = 1, upper = 100),
    class = "gaussianprocesses_bound_warning"
  )
  expect_identical(bounded$optimization$at_bound, "variance (upper)")
  expect_equal(bounded$kernel$parameters$variance, 100)

  # Without a binding bound the estimate is still very large.
  expect_no_warning(
    unbounded <- optimize_latent_gp(x, y, rbf_kernel(),
                                    bernoulli_likelihood(), n_starts = 1)
  )
  expect_identical(unbounded$optimization$at_bound, character(0))
  expect_gt(unbounded$kernel$parameters$variance, 1000)

  # Classes that overlap give a moderate variance.
  data <- classification_data(n = 100)
  overlapping <- optimize_latent_gp(data$x, data$y, rbf_kernel(),
                                    bernoulli_likelihood(), n_starts = 1)
  expect_lt(overlapping$kernel$parameters$variance, 100)
})
