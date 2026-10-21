# Poisson count models: exposure, estimated means, predictions, and scores.

count_data <- function(n = 60, seed = 41) {
  set.seed(seed)
  x <- sort(stats::runif(n, 0, 10))
  exposure <- stats::runif(n, 0.5, 2)

  list(
    x = x,
    exposure = exposure,
    y = stats::rpois(n, exposure * exp(1.5 + 0.8 * sin(x)))
  )
}


# F(k) by integrate(), split at the point where the Poisson distribution
# function falls, so that the quadrature does not miss it.
reference_count_distribution <- function(k, mean, variance, exposure = 1) {
  sd <- sqrt(variance)
  limits <- c(mean - 12 * sd, mean + 12 * sd)
  breaks <- sort(unique(pmin(pmax(
    c(limits, log((k + 1) / exposure) + c(-0.1, 0.1)),
    limits[1]
  ), limits[2])))

  sum(vapply(
    seq_len(length(breaks) - 1L),
    function(j) {
      stats::integrate(
        function(f) {
          stats::ppois(k, exposure * exp(f)) * stats::dnorm(f, mean, sd)
        },
        breaks[j],
        breaks[j + 1L],
        rel.tol = 1e-13,
        abs.tol = 1e-16,
        subdivisions = 2000L
      )$value
    },
    numeric(1)
  ))
}


test_that("exposure is an offset: doubling it and shifting the mean by -log 2 gives the same fit", {
  data <- count_data()
  kernel <- rbf_kernel(0.6, 1.2)
  x_new <- c(0.5, 4, 9.5)

  base <- fit_latent_gp(data$x, data$y, kernel, poisson_likelihood(),
                        mean = constant_mean(1.5), exposure = data$exposure)
  doubled <- fit_latent_gp(data$x, data$y, kernel, poisson_likelihood(),
                           mean = constant_mean(1.5 - log(2)),
                           exposure = 2 * data$exposure)

  expect_equal(doubled$mode, base$mode - log(2), tolerance = 1e-10)
  expect_equal(doubled$log_marginal_likelihood, base$log_marginal_likelihood,
               tolerance = 1e-12)

  # Predictions at matching exposures are identical.
  base_prediction <- predict_latent_gp(base, x_new, exposure = 1.5)
  doubled_prediction <- predict_latent_gp(doubled, x_new, exposure = 3)
  for (field in c("response_mean", "response_variance", "rate_interval",
                  "prediction_interval", "latent_variance")) {
    expect_equal(doubled_prediction[[field]], base_prediction[[field]],
                 tolerance = 1e-10, label = field)
  }

  # An estimated mean level absorbs the shift.
  estimated <- fit_latent_gp(data$x, data$y, kernel, poisson_likelihood(),
                             mean = constant_mean(estimate_coefficients()),
                             exposure = data$exposure)
  estimated_doubled <- fit_latent_gp(
    data$x, data$y, kernel, poisson_likelihood(),
    mean = constant_mean(estimate_coefficients()),
    exposure = 2 * data$exposure
  )
  expect_equal(coef(estimated_doubled), coef(estimated) - log(2),
               tolerance = 1e-6)
  expect_equal(estimated_doubled$log_marginal_likelihood,
               estimated$log_marginal_likelihood, tolerance = 1e-10)
})


test_that("estimated mean coefficients maximize the approximate marginal likelihood", {
  data <- count_data()
  kernel <- rbf_kernel(0.6, 1.2)

  model <- fit_latent_gp(data$x, data$y, kernel, poisson_likelihood(),
                         mean = linear_mean(), exposure = data$exposure)
  expect_named(coef(model), c("intercept", "x1"))
  expect_true(model$mean_fit$optimization$converged)
  expect_identical(attr(logLik(model), "df"), 2L)
  expect_output(print(model), "mean coefficients: intercept=")

  # The profile gradient omits the coefficients; their gradient is zero at
  # the estimate and matches finite differences elsewhere.
  expect_named(log_marginal_likelihood_gradient(model),
               c("variance", "length_scale"))
  expect_lt(max(abs(.laplace_gradient(model, include_mean = TRUE)[
    c("mean.intercept", "mean.x1")
  ])), 1e-4)

  for (likelihood in list(poisson_likelihood(), bernoulli_likelihood())) {
    y <- if (identical(likelihood$family, "poisson")) {
      data$y
    } else {
      as.numeric(data$y > stats::median(data$y))
    }
    exposure <- if (identical(likelihood$family, "poisson")) data$exposure
    coefficients <- c(intercept = 0.7, x1 = 0.05)
    at <- function(values) {
      .fit_latent_model(
        model$x, y, NULL, kernel, likelihood, model$mean,
        if (is.null(exposure)) rep(1, 60) else exposure,
        model$numerical_policy[1:5], 100L, 1e-14,
        coefficients = values
      )
    }
    analytical <- .laplace_gradient(at(coefficients), include_mean = TRUE)[
      c("mean.intercept", "mean.x1")
    ]
    numerical <- vapply(1:2, function(j) {
      step <- 1e-5
      up <- coefficients
      down <- coefficients
      up[j] <- up[j] + step
      down[j] <- down[j] - step
      (at(up)$log_marginal_likelihood - at(down)$log_marginal_likelihood) /
        (2 * step)
    }, numeric(1))
    expect_equal(unname(analytical), numerical, tolerance = 1e-6)
  }

  # With a Gaussian likelihood the estimate is the generalized least-squares
  # estimate of the exact model.
  continuous <- 2 + 0.3 * data$x + sin(data$x) +
    stats::rnorm(60, sd = 0.2)
  latent <- fit_latent_gp(data$x, continuous, rbf_kernel(1, 1),
                          gaussian_likelihood(0.04), mean = linear_mean())
  exact <- fit_gp(data$x, continuous, rbf_kernel(1, 1), noise_variance = 0.04,
                  mean = linear_mean(estimate_coefficients("ml")))
  expect_equal(coef(latent), coef(exact), tolerance = 1e-6)
  expect_equal(latent$log_marginal_likelihood, log_marginal_likelihood(exact),
               tolerance = 1e-10)

  # Under a Gaussian prior, coef() is the posterior mean of the exact model.
  prior <- linear_mean(coefficient_prior(c(1, 0), c(4, 1)))
  expect_equal(
    coef(fit_latent_gp(data$x, continuous, rbf_kernel(1, 1),
                       gaussian_likelihood(0.04), mean = prior)),
    coef(fit_gp(data$x, continuous, rbf_kernel(1, 1), noise_variance = 0.04,
                mean = prior)),
    tolerance = 1e-8
  )
  expect_identical(
    coef(fit_latent_gp(data$x, data$y, kernel, poisson_likelihood(),
                       mean = constant_mean(1.5))),
    c(intercept = 1.5)
  )

  # Joint optimization reaches the estimate of fit_latent_gp() at its optimum.
  optimized <- optimize_latent_gp(data$x, data$y, rbf_kernel(),
                                  poisson_likelihood(),
                                  mean = constant_mean(estimate_coefficients()),
                                  exposure = data$exposure, n_starts = 1)
  expect_named(optimized$optimization$optimized_parameters,
               c("variance", "length_scale", "mean.intercept"))
  expect_identical(attr(logLik(optimized), "df"), 3L)
  refit <- fit_latent_gp(data$x, data$y, optimized$kernel,
                         poisson_likelihood(),
                         mean = constant_mean(estimate_coefficients()),
                         exposure = data$exposure)
  expect_equal(coef(refit), coef(optimized), tolerance = 1e-4)
  expect_equal(refit$log_marginal_likelihood,
               optimized$log_marginal_likelihood, tolerance = 1e-8)

  expect_error(
    fit_latent_gp(data$x, data$y, kernel, poisson_likelihood(),
                  mean = constant_mean(estimate_coefficients("reml"))),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("the count distribution function is accurate for small and large rates", {
  cases <- expand.grid(mean = log(c(0.3, 5, 200, 5000)),
                       variance = c(0.001, 0.09, 1, 4))

  for (i in seq_len(nrow(cases))) {
    mean <- cases$mean[i]
    variance <- cases$variance[i]
    counts <- unique(pmax(0, round(exp(mean + c(-2, 0, 2) * sqrt(variance)))))
    expected <- vapply(counts, reference_count_distribution, numeric(1),
                       mean = mean, variance = variance)
    expect_lt(
      max(abs(
        .count_distribution(counts, rep(mean, length(counts)),
                            rep(variance, length(counts)),
                            rep(1, length(counts)), 40L) - expected
      )),
      1e-9,
      label = sprintf("mean %.2f, variance %g", mean, variance)
    )
  }

  expect_identical(.count_distribution(-1, 0, 1, 1, 40L), 0)
})


test_that("Poisson predictions give rate intervals, count intervals, and rate draws", {
  data <- count_data()
  model <- fit_latent_gp(data$x, data$y, rbf_kernel(0.6, 1.2),
                         poisson_likelihood(), mean = constant_mean(1.5),
                         exposure = data$exposure)
  x_new <- c(0.5, 3.3, 7, 13)
  exposure <- c(1, 2, 0.5, 1)
  prediction <- predict_latent_gp(model, x_new, interval_level = 0.9,
                                  exposure = exposure)

  expect_equal(prediction$rate_interval,
               exposure * exp(prediction$latent_interval))
  expect_equal(prediction$response_mean,
               exposure * exp(prediction$mean + prediction$latent_variance / 2))
  expect_null(prediction$probability)

  # The count intervals are the equal-tailed quantiles of the predictive
  # distribution, so they cover at least the level and are the shortest
  # such quantiles.
  for (i in seq_along(x_new)) {
    counts <- 0:400
    probabilities <- exp(likelihood_predictive(
      poisson_likelihood(),
      rep(prediction$mean[i], length(counts)),
      rep(prediction$latent_variance[i], length(counts)),
      observed = counts,
      exposure = exposure[i]
    )$log_density)
    distribution <- cumsum(probabilities)
    lower <- unname(prediction$prediction_interval[i, "lower"])
    upper <- unname(prediction$prediction_interval[i, "upper"])

    expect_identical(lower, as.numeric(min(counts[distribution >= 0.05])))
    expect_identical(upper, as.numeric(min(counts[distribution >= 0.95])))
    expect_gte(sum(probabilities[counts >= lower & counts <= upper]), 0.9)
  }

  # Gaussian likelihoods give the prediction interval of the exact model.
  continuous <- sin(data$x)
  expect_equal(
    predict_latent_gp(fit_latent_gp(data$x, continuous, rbf_kernel(),
                                    gaussian_likelihood(0.04)),
                      x_new)$prediction_interval,
    predict_gp(fit_gp(data$x, continuous, rbf_kernel(),
                      noise_variance = 0.04), x_new)$prediction_interval,
    tolerance = 1e-10
  )

  # Rate draws are the exponential of the latent draws times the exposure,
  # and their quantiles match the rate interval.
  draws <- sample_gp_posterior(model, x_new, n_draws = n_statistical_draws,
                               exposure = exposure, seed = 3)
  expect_equal(draws$rate, exposure * exp(draws$latent))
  outside <- c(
    rowMeans(draws$rate < prediction$rate_interval[, "lower"]),
    rowMeans(draws$rate > prediction$rate_interval[, "upper"])
  )
  expect_lt(max(abs(outside - 0.05)),
            5 * sqrt(0.05 * 0.95 / n_statistical_draws))
})


test_that("count scores follow their definitions", {
  data <- count_data()
  training <- seq(1, 60, by = 2)
  model <- fit_latent_gp(data$x[training], data$y[training],
                         rbf_kernel(0.6, 1.2), poisson_likelihood(),
                         mean = constant_mean(1.5),
                         exposure = data$exposure[training])
  held_out <- setdiff(seq_len(60), training)
  prediction <- predict_latent_gp(model, data$x[held_out],
                                  exposure = data$exposure[held_out])
  observed <- data$y[held_out]

  scores <- gp_count_scores(observed, prediction$mean,
                            prediction$latent_variance,
                            exposure = data$exposure[held_out],
                            interval_levels = c(0.5, 0.9), seed = 5)
  expect_s3_class(scores, "gaussianprocesses_count_scores")

  log_density <- likelihood_predictive(
    poisson_likelihood(), prediction$mean, prediction$latent_variance,
    observed = observed, exposure = data$exposure[held_out]
  )$log_density
  residual <- (observed - prediction$response_mean) /
    sqrt(prediction$response_variance)
  expect_equal(scores$summary[["log_score"]], mean(log_density))
  expect_equal(scores$summary[["dispersion"]], mean(residual^2))
  expect_equal(
    scores$summary[["dawid_sebastiani"]],
    mean(residual^2 + log(prediction$response_variance))
  )
  expect_equal(scores$calibration$level, c(0.5, 0.9))

  ninety <- predict_latent_gp(model, data$x[held_out], interval_level = 0.9,
                              exposure = data$exposure[held_out])
  expect_equal(
    scores$calibration$coverage[2],
    mean(observed >= ninety$prediction_interval[, "lower"] &
           observed <= ninety$prediction_interval[, "upper"])
  )

  # Randomized PIT values lie between F(y - 1) and F(y), reproducibly.
  below <- .count_distribution(observed - 1, prediction$mean,
                               prediction$latent_variance,
                               data$exposure[held_out], 40L)
  pit <- scores$pointwise$pit
  expect_true(all(pit >= below - 1e-12 & pit <= below + exp(log_density) +
                    1e-12))
  set.seed(1)
  before <- .Random.seed
  again <- gp_count_scores(observed, prediction$mean,
                           prediction$latent_variance,
                           exposure = data$exposure[held_out],
                           interval_levels = c(0.5, 0.9), seed = 5)
  expect_identical(.Random.seed, before)
  expect_identical(again, scores)

  # gp_holdout_scores() passes the prediction exposure through.
  expect_identical(
    gp_holdout_scores(model, data$x[held_out], observed,
                      interval_levels = c(0.5, 0.9), seed = 5,
                      exposure = data$exposure[held_out]),
    scores
  )

  output <- utils::capture.output(print(scores))
  expect_match(output[1], "Count scores for 30 observations")
  expect_true(any(grepl("dispersion", output)))

  for (call in list(
    function() gp_count_scores(observed + 0.5, prediction$mean,
                               prediction$latent_variance),
    function() gp_count_scores(observed, prediction$mean[-1],
                               prediction$latent_variance),
    function() gp_count_scores(observed, prediction$mean,
                               -prediction$latent_variance),
    function() gp_count_scores(observed, prediction$mean,
                               prediction$latent_variance, exposure = -1),
    function() gp_count_scores(observed, prediction$mean,
                               prediction$latent_variance,
                               interval_levels = 1),
    function() gp_holdout_scores(model, data$x[held_out], observed[-1])
  )) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }
})


test_that("randomized PIT values are uniform under the true model", {
  # One test count per data set, so that the PIT values are independent.
  set.seed(42)
  kernel <- rbf_kernel(0.5, 1.5)
  pit <- vapply(seq_len(150), function(replicate) {
    x <- stats::runif(31, 0, 10)
    f <- log(5) + drop(sample_gp_prior(x, kernel)$latent)
    exposure <- stats::runif(31, 0.5, 2)
    y <- stats::rpois(31, exposure * exp(f))
    model <- fit_latent_gp(x[-31], y[-31], kernel, poisson_likelihood(),
                           mean = constant_mean(log(5)),
                           exposure = exposure[-31])
    gp_holdout_scores(model, x[31], y[31], interval_levels = 0.9,
                      exposure = exposure[31])$pointwise$pit
  }, numeric(1))

  expect_gt(stats::ks.test(pit, "punif")$p.value, 0.01)
  expect_lt(abs(mean(pit) - 0.5), 5 * sqrt(1 / 12 / 150))
})


test_that("the dispersion check detects overdispersed counts", {
  set.seed(43)
  kernel <- rbf_kernel(0.5, 1.5)
  x <- stats::runif(260, 0, 10)
  f <- log(5) + drop(sample_gp_prior(x, kernel)$latent)
  training <- 1:60
  held_out <- 61:260
  dispersion <- function(rate) {
    y <- stats::rpois(260, rate)
    model <- fit_latent_gp(x[training], y[training], kernel,
                           poisson_likelihood(), mean = constant_mean(log(5)))
    gp_holdout_scores(model, x[held_out], y[held_out],
                      interval_levels = 0.9)$summary[["dispersion"]]
  }

  poisson <- dispersion(exp(f))
  overdispersed <- dispersion(exp(f) * stats::rgamma(260, 2, 2))
  expect_lt(abs(poisson - 1), 0.3)
  expect_gt(overdispersed, 2)
})
