# Likelihood specifications, derivatives, and predictive integrals.

likelihood_cases <- function() {
  list(
    list(likelihood = gaussian_likelihood(0.3), y = c(-1.2, 0.7)),
    list(likelihood = bernoulli_likelihood("logit"), y = c(0, 1)),
    list(likelihood = bernoulli_likelihood("probit"), y = c(0, 1)),
    list(likelihood = poisson_likelihood(), y = c(0, 3, 25)),
    list(likelihood = poisson_likelihood(), y = c(0, 4), exposure = 2.5)
  )
}


# Richardson-extrapolated central difference of g at f.
richardson_derivative <- function(g, f) {
  step <- 1e-3 * pmax(1, abs(f))
  central <- function(h) (g(f + h) - g(f - h)) / (2 * h)

  (4 * central(step / 2) - central(step)) / 3
}


test_that("likelihoods construct, print, and validate", {
  print_line <- function(likelihood) {
    utils::capture.output(print(likelihood))
  }

  expect_identical(
    print_line(gaussian_likelihood(0.1)),
    "GaussianLikelihood(variance=0.1)"
  )
  expect_identical(print_line(bernoulli_likelihood()),
                   "BernoulliLikelihood(link=logit)")
  expect_identical(print_line(bernoulli_likelihood("probit")),
                   "BernoulliLikelihood(link=probit)")
  expect_identical(print_line(poisson_likelihood()),
                   "PoissonLikelihood(link=log)")
  utils::capture.output(printed <- withVisible(print(poisson_likelihood())))
  expect_false(printed$visible)
  expect_s3_class(gaussian_likelihood(), "gaussianprocesses_likelihood")
  expect_identical(gaussian_likelihood()$parameters, list(variance = 1))

  for (likelihood in lapply(likelihood_cases(), `[[`, "likelihood")) {
    expect_true(.likelihood_entry(likelihood)$log_concave)
  }

  for (variance in list(0, -1, c(1, 2), NA_real_, "1")) {
    expect_error(
      gaussian_likelihood(variance),
      class = "gaussianprocesses_validation_error"
    )
  }

  expect_error(bernoulli_likelihood("cloglog"), "'arg'")
  expect_error(poisson_likelihood("identity"), "'arg'")

  # Specifications restored from files are validated before use.
  corrupted <- list(
    list(family = "bernoulli", link = "logit", parameters = list()),
    structure(
      list(family = "binomial", link = "logit", parameters = list()),
      class = "gaussianprocesses_likelihood"
    ),
    structure(
      list(family = "bernoulli", link = c("logit", "probit"),
           parameters = list()),
      class = "gaussianprocesses_likelihood"
    ),
    structure(
      list(family = "gaussian", link = "identity", parameters = list()),
      class = "gaussianprocesses_likelihood"
    ),
    structure(
      list(family = "gaussian", link = "identity",
           parameters = list(variance = -1)),
      class = "gaussianprocesses_likelihood"
    ),
    structure(
      list(family = "poisson", link = "log",
           parameters = list(dispersion = 2)),
      class = "gaussianprocesses_likelihood"
    )
  )

  for (likelihood in corrupted) {
    expect_error(
      evaluate_likelihood(likelihood, 1, 0),
      class = "gaussianprocesses_validation_error"
    )
  }
})


test_that("likelihood specifications survive saveRDS()", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))

  for (case in likelihood_cases()) {
    likelihood <- case$likelihood

    # Data only: no functions or environments to serialize.
    expect_false(any(unlist(rapply(
      unclass(likelihood),
      function(x) is.function(x) || is.environment(x),
      how = "unlist"
    ))))

    saveRDS(likelihood, path)
    restored <- readRDS(path)
    expect_identical(restored, likelihood)
    expect_identical(
      evaluate_likelihood(restored, case$y, seq_along(case$y) / 2,
                          exposure = case$exposure),
      evaluate_likelihood(likelihood, case$y, seq_along(case$y) / 2,
                          exposure = case$exposure)
    )
  }
})


test_that("log densities match their reference implementations", {
  f <- c(-40, -8, -1, 0, 0.5, 3, 40)
  ones <- rep(1, length(f))
  zeros <- rep(0, length(f))

  probit <- bernoulli_likelihood("probit")
  expect_identical(
    evaluate_likelihood(probit, ones, f)$log_density,
    stats::pnorm(f, log.p = TRUE)
  )
  expect_identical(
    evaluate_likelihood(probit, zeros, f)$log_density,
    stats::pnorm(-f, log.p = TRUE)
  )
  expect_equal(evaluate_likelihood(probit, 1, -40)$log_density, -804.608442,
               tolerance = 1e-9)

  # log sigma(f), computed without overflow on either side.
  log_logistic <- function(f) ifelse(f > 0, -log1p(exp(-f)), f - log1p(exp(f)))
  logit <- bernoulli_likelihood()
  expect_equal(evaluate_likelihood(logit, ones, f)$log_density,
               log_logistic(f), tolerance = 1e-14)
  expect_equal(evaluate_likelihood(logit, zeros, f)$log_density,
               log_logistic(-f), tolerance = 1e-14)

  counts <- c(0, 2, 5, 1, 0, 30, 7)
  exposure <- c(1, 0.5, 2, 3, 1, 1, 0.1)
  expect_equal(
    evaluate_likelihood(poisson_likelihood(), counts, f,
                        exposure = exposure)$log_density,
    stats::dpois(counts, exposure * exp(f), log = TRUE),
    tolerance = 1e-12
  )

  responses <- c(-1, 0, 2, 0.5, 0, 3, -0.2)
  expect_equal(
    evaluate_likelihood(gaussian_likelihood(0.4), responses, f)$log_density,
    stats::dnorm(responses, f, sqrt(0.4), log = TRUE),
    tolerance = 1e-14
  )
})


test_that("likelihood derivatives pass finite-difference checks to |f| = 40", {
  grid <- c(-40, -30, -20, -10, -5.5, -5, -4.5, -2, -0.3, 0, 0.3, 2, 4.5, 5,
            5.5, 10, 20, 30, 40)

  for (case in likelihood_cases()) {
    for (response in case$y) {
      y <- rep(response, length(grid))
      evaluate <- function(f) {
        evaluate_likelihood(case$likelihood, y, f, exposure = case$exposure)
      }
      values <- evaluate(grid)
      label <- paste(
        utils::capture.output(print(case$likelihood)),
        "y =",
        response
      )

      expect_true(all(is.finite(as.matrix(values))), label = label)

      for (pair in list(
        c("log_density", "first"),
        c("first", "second"),
        c("second", "third")
      )) {
        numerical <- richardson_derivative(
          function(f) evaluate(f)[[pair[1]]],
          grid
        )
        error <- abs(values[[pair[2]]] - numerical) /
          (abs(numerical) + 1e-8 * (1 + abs(values[[pair[1]]])))

        expect_lt(max(error), 1e-6, label = paste(label, pair[2]))
      }
    }
  }
})


test_that("Gaussian parameter derivatives match finite differences", {
  likelihood <- gaussian_likelihood(0.3)
  entry <- .likelihood_entry(likelihood)
  y <- c(-1.2, 0.7, 2)
  f <- c(0.1, 0.7, -1)
  derivatives <- entry$parameter_derivatives(y, f, likelihood$parameters, 1)

  at_log_variance <- function(log_variance, quantity) {
    parameters <- list(variance = exp(log_variance))
    switch(
      quantity,
      log_density = entry$log_density(y, f, parameters, 1),
      entry$derivatives(y, f, parameters, 1)[[quantity]]
    )
  }

  for (quantity in c("log_density", "first", "second")) {
    numerical <- vapply(
      seq_along(y),
      function(i) {
        richardson_derivative(
          function(theta) at_log_variance(theta, quantity)[i],
          log(0.3)
        )
      },
      numeric(1)
    )

    expect_equal(derivatives$variance[[quantity]], numerical,
                 tolerance = 1e-8, label = quantity)
  }
})


test_that("binary and count responses are coded and validated", {
  logit <- bernoulli_likelihood()
  f <- c(-1, 0.5, 2)

  coded <- evaluate_likelihood(logit, c(0, 1, 1), f)
  expect_identical(coded$y, c(0, 1, 1))
  expect_identical(
    evaluate_likelihood(logit, c(FALSE, TRUE, TRUE), f),
    coded
  )

  # The second factor level is the event, as in glm().
  expect_identical(
    evaluate_likelihood(
      logit,
      factor(c("no", "yes", "yes"), levels = c("no", "yes")),
      f
    ),
    coded
  )
  expect_identical(
    evaluate_likelihood(
      logit,
      factor(c("yes", "no", "no"), levels = c("yes", "no")),
      f
    ),
    coded
  )

  invalid <- list(
    function() evaluate_likelihood(logit, c(0, 2, 1), f),
    function() evaluate_likelihood(logit, c(TRUE, NA, FALSE), f),
    function() evaluate_likelihood(logit, factor(c("a", "b", "c")), f),
    function() evaluate_likelihood(logit, factor(c("a", "b", NA)), f),
    function() evaluate_likelihood(logit, c("0", "1", "1"), f),
    function() evaluate_likelihood(logit, c(0, 1), f),
    function() evaluate_likelihood(logit, c(0, 1, 1), c(0, NA, 1)),
    function() evaluate_likelihood(logit, c(0, 1, 1), f, exposure = 2),
    function() evaluate_likelihood(poisson_likelihood(), c(0, -1, 2), f),
    function() evaluate_likelihood(poisson_likelihood(), c(0, 1.5, 2), f),
    function() evaluate_likelihood(poisson_likelihood(), c(0, 1, 2), f,
                                   exposure = c(1, 2)),
    function() evaluate_likelihood(poisson_likelihood(), c(0, 1, 2), f,
                                   exposure = 0),
    function() evaluate_likelihood(gaussian_likelihood(), c(0, Inf, 1), f),
    function() evaluate_likelihood(gaussian_likelihood(), 1, 0,
                                   exposure = 1)
  )

  for (call in invalid) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }

  # Exposure is one value or one per response.
  expect_identical(
    evaluate_likelihood(poisson_likelihood(), c(0, 1, 2), f, exposure = 2),
    evaluate_likelihood(poisson_likelihood(), c(0, 1, 2), f,
                        exposure = c(2, 2, 2))
  )
})


test_that("predictive distributions match their closed forms", {
  mean <- c(-1.5, 0, 0.8, 2)
  variance <- c(0.2, 1, 0.5, 1.5)

  gaussian <- likelihood_predictive(gaussian_likelihood(0.3), mean, variance,
                                    observed = c(-1, 0.5, 0.8, 4))
  expect_equal(gaussian$mean, mean)
  expect_equal(gaussian$variance, variance + 0.3)
  expect_equal(
    gaussian$log_density,
    stats::dnorm(c(-1, 0.5, 0.8, 4), mean, sqrt(variance + 0.3), log = TRUE)
  )

  # The probit probability Phi(mu / sqrt(1 + sigma^2)) is E[Phi(f)].
  probit <- likelihood_predictive(bernoulli_likelihood("probit"), mean,
                                  variance, observed = c(1, 0, 1, 0))
  expect_equal(
    probit$probability,
    .gaussian_expectation(function(f, index) stats::pnorm(f), mean, variance,
                          40),
    tolerance = 1e-12
  )
  expect_equal(probit$mean, probit$probability)
  expect_equal(probit$variance, probit$probability * (1 - probit$probability))
  expect_equal(
    probit$log_density,
    log(c(probit$probability[1], 1 - probit$probability[2],
          probit$probability[3], 1 - probit$probability[4])),
    tolerance = 1e-12
  )

  # Poisson with a lognormal rate.
  exposure <- c(1, 2, 0.5, 3)
  poisson <- likelihood_predictive(poisson_likelihood(), mean, variance,
                                   exposure = exposure)
  rate_moment <- function(power) {
    .gaussian_expectation(
      function(f, index) (exposure[index] * exp(f))^power,
      mean,
      variance,
      40
    )
  }
  expect_equal(poisson$mean, rate_moment(1), tolerance = 1e-12)
  expect_equal(
    poisson$variance,
    rate_moment(1) + rate_moment(2) - rate_moment(1)^2,
    tolerance = 1e-10
  )
  expect_null(poisson$probability)
  expect_null(poisson$log_density)
})


test_that("predictive densities without a closed form match integration", {
  integrate_density <- function(density, mean, variance) {
    sd <- sqrt(variance)
    stats::integrate(
      function(f) density(f) * stats::dnorm(f, mean, sd),
      mean - 20 * sd,
      mean + 20 * sd,
      rel.tol = 1e-12,
      subdivisions = 1000L
    )$value
  }

  # The third count lies in the tail of its latent distribution, where
  # quadrature centred at the mean loses several digits.
  mean <- c(log(3), 0, 2, -4)
  variance <- c(0.2, 1, 0.5, 2)
  observed <- c(4, 0, 12, 3)
  exposure <- c(2, 1, 0.5, 1)
  poisson <- likelihood_predictive(poisson_likelihood(), mean, variance,
                                   observed = observed, exposure = exposure)
  expected <- vapply(
    seq_along(mean),
    function(i) {
      log(integrate_density(
        function(f) stats::dpois(observed[i], exposure[i] * exp(f)),
        mean[i],
        variance[i]
      ))
    },
    numeric(1)
  )
  expect_equal(poisson$log_density, expected, tolerance = 1e-9)

  # The predictive probabilities of all counts sum to one.
  counts <- 0:150
  total <- sum(exp(likelihood_predictive(
    poisson_likelihood(),
    rep(log(3), length(counts)),
    rep(0.2, length(counts)),
    observed = counts,
    exposure = 2
  )$log_density))
  expect_equal(total, 1, tolerance = 1e-10)

  mean <- c(-6, -1, 0, 1.5, 5)
  variance <- c(0.5, 2, 1, 0.3, 3)
  logit <- likelihood_predictive(bernoulli_likelihood(), mean, variance,
                                 observed = c(1, 0, 1, 1, 0))
  probability <- vapply(
    seq_along(mean),
    function(i) integrate_density(stats::plogis, mean[i], variance[i]),
    numeric(1)
  )
  expect_equal(logit$probability, probability, tolerance = 1e-8)
  expect_equal(
    logit$log_density,
    log(c(probability[1], 1 - probability[2], probability[3:4],
          1 - probability[5])),
    tolerance = 1e-8
  )
  expect_equal(logit$variance, logit$probability * (1 - logit$probability))

  # Zero variance gives the likelihood itself.
  degenerate <- likelihood_predictive(poisson_likelihood(), c(0, 1), c(0, 0),
                                      observed = c(2, 3))
  expect_equal(degenerate$log_density,
               stats::dpois(c(2, 3), exp(c(0, 1)), log = TRUE))
})


test_that("likelihood_predictive() validates its arguments", {
  logit <- bernoulli_likelihood()

  for (call in list(
    function() likelihood_predictive(logit, c(0, 1), c(1, -1)),
    function() likelihood_predictive(logit, c(0, 1), 1),
    function() likelihood_predictive(logit, c(0, 1), c(1, 1), observed = 1),
    function() likelihood_predictive(logit, c(0, 1), c(1, 1),
                                     observed = c(0, 2)),
    function() likelihood_predictive(logit, c(0, 1), c(1, 1),
                                     quadrature_order = 0),
    function() likelihood_predictive(logit, c(0, 1), c(1, 1), exposure = 2),
    function() likelihood_predictive("logit", 0, 1)
  )) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }
})
