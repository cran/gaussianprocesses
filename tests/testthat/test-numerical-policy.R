# Scale invariance under the relative jitter policy. Multiplying y by c and
# every variance parameter by c^2 multiplies the covariance by c^2; with
# jitter proportional to the covariance, the jittered covariance scales the
# same way, so in exact arithmetic
#   mean -> c mean, variance -> c^2 variance, draws -> c draws,
#   log marginal likelihood -> original - n log|c|, log-scale gradient fixed.
# inst/benchmarks/jitter-scaling.R compares both jitter policies.
#
# An exactly singular covariance sits on the boundary of positive
# definiteness, where rounding decides whether a factorization without
# jitter succeeds, and that decision differs between platforms and scales.
# The tests therefore use problems whose jitter does not depend on it: a
# singular covariance with jitter already in the first attempt, and a dense
# grid whose many negligible eigenvalues make every scale need the ladder.

relative_difference <- function(actual, expected) {
  max(abs(actual - expected)) / max(abs(expected))
}

# Duplicated inputs without noise, with relative jitter 1e-10 from the first
# attempt, so the same jitter applies at every scale on every platform.
invariance_x <- c(-1, 0, 0, 0.4, 1.1, 1.1, 2)
invariance_x_new <- c(-0.5, 0.3, 0.75, 2.5)

scaled_problem <- function(c) {
  y <- c * (sin(3 * invariance_x) + 0.2 * invariance_x)
  kernel <- rbf_kernel(1.3 * c^2, 0.8)
  model <- fit_gp(
    invariance_x,
    y,
    kernel,
    noise_variance = 0,
    initial_jitter = 1e-10
  )
  prediction <- predict_gp(model, invariance_x_new)

  list(
    model = model,
    log_marginal_likelihood = log_marginal_likelihood(model) +
      length(invariance_x) * log(c),
    gradient = log_marginal_likelihood_gradient(model),
    mean = prediction$mean / c,
    variance = prediction$latent_variance / c^2,
    prior_draws = sample_gp_prior(
      invariance_x,
      kernel,
      n_draws = 3,
      seed = 1,
      initial_jitter = 1e-10
    )$latent / c,
    posterior_draws = sample_gp_posterior(
      model,
      invariance_x_new,
      n_draws = 3,
      seed = 2,
      initial_jitter = 1e-10
    )$latent / c
  )
}


test_that("results do not depend on the units of y when jitter is applied", {
  reference <- scaled_problem(1)
  expect_gt(reference$model$jitter, 0)

  for (c in c(1e-4, 1e4)) {
    scaled <- scaled_problem(c)

    expect_equal(scaled$model$jitter, c^2 * reference$model$jitter)

    # Agreement is about 1e-11 or better on the platforms tested.
    for (field in c("mean", "variance", "prior_draws", "posterior_draws")) {
      expect_lt(
        relative_difference(scaled[[field]], reference[[field]]),
        1e-8,
        label = paste(field, "at c =", c)
      )
    }

    # The log determinant includes logs of tiny pivots of the near-singular
    # covariance, which carry the rounding of the inexact scaling by c^2;
    # agreement is about 2e-8 to 7e-8 on the platforms tested.
    for (field in c("log_marginal_likelihood", "gradient")) {
      expect_lt(
        relative_difference(scaled[[field]], reference[[field]]),
        1e-6,
        label = paste(field, "at c =", c)
      )
    }
  }
})


test_that("the jitter ladder chooses the same relative jitter at every scale", {
  # 40 nearly coincident inputs make C numerically rank deficient at any
  # scale, so every fit climbs the ladder.
  x <- seq(0, 1, length.out = 40)
  x_new <- c(-0.5, 0.3, 0.75, 2.5)
  fit_scaled <- function(c) {
    fit_gp(x, c * sin(3 * x), rbf_kernel(1.3 * c^2, 0.5), noise_variance = 0)
  }

  reference <- fit_scaled(1)
  reference_prediction <- predict_gp(reference, x_new)
  expect_gt(reference$jitter, 0)

  for (c in c(1e-4, 1e4)) {
    model <- fit_scaled(c)
    prediction <- predict_gp(model, x_new)

    expect_identical(model$cholesky_attempts, reference$cholesky_attempts)
    expect_equal(model$jitter, c^2 * reference$jitter)

    # The near-singular covariance amplifies the rounding of the inexact
    # scaling by c^2; the agreement is about 1e-8 on the platforms tested.
    expect_lt(
      relative_difference(prediction$mean / c, reference_prediction$mean),
      1e-7
    )
    expect_lt(
      relative_difference(
        prediction$latent_variance / c^2,
        reference_prediction$latent_variance
      ),
      1e-7
    )
  }
})


test_that("absolute jitter, the previous policy, depends on the units of y", {
  old <- options(gaussianprocesses.jitter_scaling = "absolute")
  on.exit(options(old), add = TRUE)

  reference <- scaled_problem(1)
  scaled <- scaled_problem(1e-4)

  # An absolute jitter of 1e-10 is 1e-2 of a covariance of scale 1e-8, which
  # changes the posterior materially.
  expect_gt(relative_difference(scaled$mean, reference$mean), 1e-2)
})


test_that("models, diagnostics, and samples report jitter consistently", {
  model <- scaled_problem(1e-4)$model
  diagnostics <- gp_numerical_diagnostics(model)

  expect_identical(diagnostics$numerical_jitter, model$jitter)
  expect_equal(
    diagnostics$covariance_scale,
    mean(diag(gaussianprocesses:::.model_observation_covariance(model)))
  )
  expect_equal(
    diagnostics$relative_jitter,
    model$jitter / diagnostics$covariance_scale
  )
  expect_equal(diagnostics$relative_jitter, 1e-10)
  expect_null(diagnostics$jitter_relative_to_covariance_scale)

  printed <- utils::capture.output(print(model))
  expect_true(
    any(grepl("1e-10 times the covariance scale", printed, fixed = TRUE))
  )
  expect_true(
    any(grepl(
      "1e-10 times the covariance scale",
      utils::capture.output(print(gp_diagnostics(model))),
      fixed = TRUE
    ))
  )

  draws <- sample_gp_prior(seq(0, 1, length.out = 100), rbf_kernel(4, 2))
  expect_gt(draws$jitter, 0)
  expect_equal(draws$covariance_scale, 4)
  expect_equal(draws$jitter, draws$relative_jitter * draws$covariance_scale)
  expect_true(
    any(grepl(
      "times the covariance scale 4",
      utils::capture.output(print(draws)),
      fixed = TRUE
    ))
  )
})
