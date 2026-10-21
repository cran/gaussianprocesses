# Hyperparameter uncertainty and multi-start agreement.

uncertainty_models <- function() {
  set.seed(31)
  x <- seq(0, 10, length.out = 50)
  y <- sin(x) + 0.3 * x + stats::rnorm(50, sd = 0.2)

  list(
    rbf = optimize_gp(x, y, rbf_kernel(), noise_variance = 0.1,
                      n_starts = 1),
    trend = optimize_gp(x, y, rbf_kernel(), noise_variance = 0.1,
                        mean = linear_mean(), n_starts = 1),
    changepoint = {
      kernel <- changepoint_kernel(rbf_kernel(1, 2), rbf_kernel(0.5, 0.5),
                                   location = 5, steepness = 2)
      simulated <- simulate_gp_data(seq(0, 10, length.out = 60), kernel,
                                    noise_variance = 0.01, seed = 1)
      optimize_gp(simulated$x, simulated$observed, kernel,
                  noise_variance = 0.01, n_starts = 1)
    }
  )
}

# Second differences of the log marginal likelihood in the coordinates.
likelihood_second_differences <- function(model, step = 1e-3) {
  setup <- gaussianprocesses:::.estimated_parameter_setup(model)
  eta <- setup$coordinates
  likelihood <- function(shift) {
    log_marginal_likelihood(
      gaussianprocesses:::.refit_at_coordinates(model, eta + shift, setup)
    )
  }
  unit <- function(index) {
    stats::setNames(replace(numeric(length(eta)), index, step), names(eta))
  }
  centre <- likelihood(0)
  hessian <- matrix(0, length(eta), length(eta),
                    dimnames = list(names(eta), names(eta)))

  for (i in seq_along(eta)) {
    for (j in seq_along(eta)) {
      hessian[i, j] <- if (i == j) {
        (likelihood(unit(i)) - 2 * centre + likelihood(-unit(i))) / step^2
      } else {
        (likelihood(unit(i) + unit(j)) - likelihood(unit(i) - unit(j)) -
          likelihood(unit(j) - unit(i)) + likelihood(-unit(i) - unit(j))) /
          (4 * step^2)
      }
    }
  }

  hessian
}


test_that("the information agrees with second differences of the likelihood", {
  # Slow; skipped on CRAN to keep its check within about a minute.
  skip_on_cran()
  for (name in names(uncertainty_models())) {
    model <- uncertainty_models()[[name]]
    uncertainty <- gp_hyperparameter_uncertainty(model)
    reference <- -likelihood_second_differences(model)

    expect_lte(
      max(abs(uncertainty$information - reference)) / max(abs(reference)),
      1e-4,
      label = name
    )
    # The estimate is a maximum: small gradient, positive definite.
    expect_lt(max(abs(uncertainty$gradient)), 1e-3, label = name)
    expect_true(all(uncertainty$eigenvalues > 0), label = name)
  }
})


test_that("standard errors, correlations, and intervals are consistent", {
  model <- uncertainty_models()$changepoint
  uncertainty <- gp_hyperparameter_uncertainty(model, interval_level = 0.9)
  parameters <- uncertainty$parameters

  expect_identical(
    parameters$parameter,
    c(names(kernel_parameters(model$kernel, flatten = TRUE)),
      "noise_variance")
  )
  expect_identical(
    parameters$coordinate,
    c("identity", rep("log", 6))
  )
  expect_equal(
    parameters$standard_error,
    unname(sqrt(diag(uncertainty$covariance)))
  )
  expect_equal(
    unname(uncertainty$covariance %*% uncertainty$information),
    diag(nrow(parameters)),
    tolerance = 1e-8
  )
  expect_equal(unname(diag(uncertainty$correlation)), rep(1, 7))

  # Intervals are symmetric on the coordinate scale and positive for log
  # coordinates.
  z <- stats::qnorm(0.95)
  location <- parameters[1, ]
  expect_equal(location$upper - location$estimate,
               z * location$standard_error)
  positive <- parameters[-1, ]
  expect_true(all(positive$lower > 0))
  expect_equal(log(positive$upper) - log(positive$estimate),
               z * positive$standard_error)

  output <- utils::capture.output(print(uncertainty))
  expect_match(output[1], "Hyperparameter uncertainty")
  expect_true(any(grepl("upper 90%", output, fixed = TRUE)))
})


test_that("the identifiability flag detects the Matern infill ridge", {
  fit <- function(range, length_scale) {
    x <- seq(0, range, length.out = 100)
    truth <- matern32_kernel(1, length_scale)
    simulated <- simulate_gp_data(x, truth, noise_variance = 1e-4, seed = 1)
    optimize_gp(x, simulated$observed, truth, noise_variance = 1e-4,
                optimize_noise = FALSE, n_starts = 1)
  }

  # Infill: 100 points within one length scale. Only variance /
  # length_scale^3 is identified, so the ridge direction in (log variance,
  # log length scale) is (3, 1) / sqrt(10).
  infill <- gp_hyperparameter_uncertainty(fit(1, 3))
  expect_true(infill$poorly_identified)
  expect_gt(infill$weakest_log_scale_sd, log(2))
  expect_equal(
    unname(infill$weakest_direction),
    c(3, 1) / sqrt(10),
    tolerance = 0.05
  )
  expect_true(any(grepl(
    "poorly identified along",
    utils::capture.output(print(infill))
  )))

  # Twenty length scales of data identify both.
  wide <- gp_hyperparameter_uncertainty(fit(100, 5))
  expect_false(wide$poorly_identified)
  expect_lt(wide$weakest_log_scale_sd, log(2))
  expect_false(any(grepl(
    "poorly identified",
    utils::capture.output(print(wide))
  )))

  # A redundant parameter makes the information singular.
  set.seed(2)
  x <- seq(0, 5, length.out = 30)
  redundant <- fit_gp(
    x,
    sin(x) + stats::rnorm(30, sd = 0.1),
    scale_kernel(rbf_kernel(), 2),
    noise_variance = 0.01
  )
  expect_true(gp_hyperparameter_uncertainty(redundant)$poorly_identified)
})


test_that("the profile likelihood peaks at the estimate", {
  model <- uncertainty_models()$rbf
  uncertainty <- gp_hyperparameter_uncertainty(model)
  standard_error <- uncertainty$parameters$standard_error[2]
  estimate <- kernel_parameters(model$kernel, flatten = TRUE)[["length_scale"]]

  profile <- gp_profile_likelihood(
    model,
    "length_scale",
    values = estimate * exp(c(-1, 0, 1) * standard_error)
  )

  expect_s3_class(profile, "gaussianprocesses_profile_likelihood")
  expect_equal(profile$deviance[2], 0, tolerance = 1e-8)
  expect_true(all(profile$converged))
  # One standard error away the deviance is about 1 when the likelihood is
  # nearly quadratic.
  expect_true(all(abs(profile$deviance[c(1, 3)] - 1) < 0.3))
  expect_identical(attr(profile, "parameter"), "length_scale")

  default_grid <- gp_profile_likelihood(model, "noise_variance", n_values = 3)
  expect_identical(nrow(default_grid), 3L)
  expect_lte(min(default_grid$deviance), 1e-8)

  single <- fit_gp(seq(0, 1, length.out = 10), sin(1:10), rbf_kernel(),
                   noise_variance = 0.1)
  expect_identical(
    nrow(gp_profile_likelihood(single, "variance", values = c(0.5, 1))),
    2L
  )

  expect_error(
    gp_profile_likelihood(model, "period"),
    "one of the estimated parameters",
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    gp_profile_likelihood(model, "length_scale", values = -1),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    gp_hyperparameter_uncertainty(model, step = 0),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("optimize_gp() and summary() report multi-start agreement", {
  truth <- changepoint_kernel(rbf_kernel(1, 2), rbf_kernel(1, 0.3),
                              location = 6, steepness = 5)
  x <- seq(0, 10, length.out = 80)
  simulated <- simulate_gp_data(x, truth, noise_variance = 0.01, seed = 3)

  model <- optimize_gp(
    x,
    simulated$observed,
    changepoint_kernel(rbf_kernel(), rbf_kernel(), location = 2,
                       steepness = 2),
    noise_variance = 0.05,
    n_starts = 4
  )
  optima <- model$optimization$optima

  # Starts 1 and 2 reach one optimum; 3 and 4 reach others.
  expect_identical(nrow(optima), 3L)
  expect_identical(sum(optima$n_starts), 4L)
  expect_identical(optima$gap[1], 0)
  expect_true(all(diff(optima$log_marginal_likelihood) < 0))
  expect_true(all(optima$gap[-1] > 1))
  expect_equal(optima$log_marginal_likelihood[1],
               log_marginal_likelihood(model))

  output <- utils::capture.output(print(summary(model)))
  expect_true(any(grepl(
    "starts: 4, reaching 3 distinct optima (log marginal likelihood gaps 0,",
    output,
    fixed = TRUE
  )))

  single <- uncertainty_models()$rbf
  expect_identical(nrow(single$optimization$optima), 1L)
  expect_true(any(grepl(
    "starts: 1, reaching 1 distinct optimum",
    utils::capture.output(print(summary(single))),
    fixed = TRUE
  )))

  # Starts with equal likelihoods but different parameters, or failed
  # starts, are told apart.
  results <- list(
    list(value = -10, par = c(a = 0, b = 0), convergence = 0L),
    list(value = -10, par = c(a = 1, b = 0), convergence = 0L),
    list(value = 1e100, par = c(a = 0, b = 0), convergence = 52L),
    list(value = -10 + 1e-4, par = c(a = 0, b = 0.001), convergence = 0L)
  )
  distinct <- gaussianprocesses:::.distinct_optima(results)
  expect_identical(distinct$starts, c("1, 4", "2"))
  expect_identical(distinct$gap, c(0, 0))
})
