test_that("seeded simulations are reproducible without changing caller RNG", {
  set.seed(123)
  before <- .Random.seed

  first <- simulate_gp_data(
    seq(0, 1, length.out = 10),
    kernel = rbf_kernel(length_scale = 0.2),
    noise_variance = 0.05,
    seed = 1729
  )

  after <- .Random.seed

  second <- simulate_gp_data(
    seq(0, 1, length.out = 10),
    kernel = rbf_kernel(length_scale = 0.2),
    noise_variance = 0.05,
    seed = 1729
  )

  expect_equal(first$latent, second$latent)
  expect_equal(first$observed, second$observed)
  expect_equal(before, after)
})


test_that("simulation respects supplied heteroscedastic noise vector", {
  x <- seq(0, 1, length.out = 12)
  noise <- 0.01 + 0.1 * x

  simulation <- simulate_gp_data(
    x,
    kernel = matern32_kernel(length_scale = 0.25),
    noise_variance = noise,
    seed = 8
  )

  expect_equal(simulation$noise_variance, noise)
  expect_equal(
    simulation$observed,
    simulation$latent + simulation$observation_error,
    tolerance = 1e-12
  )
})


test_that("standard scenarios expose intended structures", {
  smooth <- gp_simulation_scenario(
    "smooth_1d",
    n = 20,
    seed = 1
  )
  singular <- gp_simulation_scenario(
    "near_singular",
    n = 20,
    seed = 1
  )
  ard <- gp_simulation_scenario(
    "ard_2d",
    n = 20,
    seed = 1
  )
  hetero <- gp_simulation_scenario(
    "heteroscedastic_1d",
    n = 20,
    seed = 1
  )

  expect_equal(ncol(smooth$x), 1L)
  expect_equal(ncol(ard$x), 2L)
  expect_equal(
    length(ard$kernel$parameters$length_scale),
    2L
  )
  expect_lt(
    min(diff(sort(singular$x[, 1]))),
    1e-10
  )
  expect_gt(
    max(hetero$noise_variance),
    min(hetero$noise_variance)
  )
})


test_that("prediction comparison separates mean and variance errors", {
  reference <- list(
    mean = c(0, 1),
    latent_variance = c(0.2, 0.4)
  )
  candidate <- list(
    mean = c(0.1, 0.8),
    latent_variance = c(0.3, 0.35)
  )

  metrics <- compare_gp_predictions(
    reference,
    candidate
  )

  expect_equal(
    metrics$posterior_mean_rmse,
    sqrt(mean(c(0.1, -0.2)^2))
  )
  expect_equal(
    metrics$latent_variance_mae,
    mean(c(0.1, 0.05))
  )
})


test_that("benchmark suite keeps accuracy performance and stability separate", {
  # Slow; skipped on CRAN to keep its check within about a minute.
  skip_on_cran()
  suite <- run_gp_benchmark_suite(
    scenarios = c("smooth_1d", "ard_2d"),
    n = 12,
    n_inducing = 5,
    repeats = 1,
    seed = 42
  )

  expect_s3_class(
    suite,
    "gaussianprocesses_benchmark_suite"
  )
  expect_true(is.data.frame(suite$accuracy))
  expect_true(is.data.frame(suite$performance))
  expect_true(is.data.frame(suite$stability))

  expect_true(
    all(
      c(
        "latent_truth_rmse",
        "posterior_mean_rmse_vs_exact",
        "latent_variance_mae_vs_exact"
      ) %in% names(suite$accuracy)
    )
  )

  expect_true(
    all(
      c(
        "fit_elapsed_seconds",
        "predict_elapsed_seconds",
        "retained_model_bytes"
      ) %in% names(suite$performance)
    )
  )

  expect_true(
    all(
      c(
        "exact_numerical_jitter",
        "exact_reciprocal_condition_number"
      ) %in% names(suite$stability)
    )
  )
})


test_that("near-singular benchmark exposes conditioning separately", {
  suite <- run_gp_benchmark_suite(
    scenarios = "near_singular",
    n = 12,
    n_inducing = 5,
    repeats = 1,
    seed = 7
  )

  expect_equal(
    suite$stability$scenario,
    "near_singular"
  )
  expect_true(
    is.finite(
      suite$stability$exact_reciprocal_condition_number
    )
  )
  expect_gte(
    suite$stability$exact_numerical_jitter,
    0
  )
})


test_that("exact scaling benchmark is reproducible in structure", {
  result <- benchmark_exact_gp_scaling(
    sizes = c(8, 12),
    repeats = 1,
    seed = 99
  )

  expect_equal(result$n, c(8, 12))
  expect_true(all(result$retained_model_bytes > 0))
  expect_true(all(result$median_fit_elapsed_seconds >= 0))
  expect_true(all(result$reciprocal_condition_number > 0))
})


test_that("benchmark arguments are validated", {
  expect_error(
    gp_simulation_scenario(
      "smooth_1d",
      n = 3
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    run_gp_benchmark_suite(
      scenarios = "unknown",
      n = 10,
      n_inducing = 4
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    run_gp_benchmark_suite(
      n = 10,
      n_inducing = 11
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    benchmark_exact_gp_scaling(
      sizes = c(3, 10)
    ),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("seeds outside the integer range fail validation cleanly", {
  # Previously coerced to NA with warnings and failed inside set.seed().
  expect_warning(
    expect_error(
      simulate_gp_data(1:3, rbf_kernel(), seed = 3e9),
      class = "gaussianprocesses_validation_error"
    ),
    regexp = NA
  )

  # Derived per-scenario seeds must not overflow either.
  expect_warning(
    expect_error(
      run_gp_benchmark_suite(
        scenarios = c("smooth_1d", "ard_2d"),
        n = 8,
        n_inducing = 3,
        repeats = 1,
        seed = .Machine$integer.max
      ),
      class = "gaussianprocesses_validation_error"
    ),
    regexp = NA
  )

  expect_s3_class(
    simulate_gp_data(1:3, rbf_kernel(), seed = .Machine$integer.max),
    "gaussianprocesses_simulation"
  )
})


test_that("simulation and benchmark objects print concisely", {
  simulation <- gp_simulation_scenario(
    "smooth_1d",
    n = 10,
    seed = 1
  )
  suite <- run_gp_benchmark_suite(
    scenarios = "smooth_1d",
    n = 10,
    n_inducing = 4,
    repeats = 1,
    seed = 1
  )

  simulation_output <- capture.output(print(simulation))
  suite_output <- capture.output(print(suite))

  expect_true(
    any(grepl("Gaussian-process simulation", simulation_output, fixed = TRUE))
  )
  expect_true(
    any(grepl("benchmark suite", suite_output, fixed = TRUE))
  )
})
