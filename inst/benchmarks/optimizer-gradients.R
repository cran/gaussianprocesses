# Analytical versus numerical gradients in optimize_gp()
#
# Run manually with the package attached, for example
#   Rscript -e 'library(gaussianprocesses); source("inst/benchmarks/optimizer-gradients.R")'
# This file is intentionally excluded from CI.
#
# Each problem is simulated from a known kernel and optimized from the same
# starting kernel with both gradient methods and the default three starts.
# Elapsed times are medians over repeated runs on this machine; evaluation
# counts and log marginal likelihoods are deterministic.

repeats <- 5L

problems <- list(
  list(
    name = "1D RBF, n = 15",
    x = seq(-3, 3, length.out = 15),
    truth = rbf_kernel(1.2, 0.7),
    start = rbf_kernel(1, 1)
  ),
  list(
    name = "1D RBF, n = 60",
    x = seq(-3, 3, length.out = 60),
    truth = rbf_kernel(1.2, 0.7),
    start = rbf_kernel(1, 1)
  ),
  list(
    name = "2D ARD RBF, n = 60",
    x = cbind(seq(-2, 2, length.out = 60), cos(1:60)),
    truth = rbf_kernel(1, c(0.8, 1.5)),
    start = rbf_kernel(1, c(1, 1))
  ),
  list(
    name = "3D ARD Matern 5/2, n = 80",
    x = cbind(seq(-2, 2, length.out = 80), sin(1:80), cos(2 * (1:80))),
    truth = matern52_kernel(1.3, c(0.6, 1.2, 2.5)),
    start = matern52_kernel(1, c(1, 1, 1))
  ),
  list(
    name = "2D ARD RBF + periodic, n = 60",
    x = cbind(seq(-2, 2, length.out = 60), cos(1:60)),
    truth = sum_kernel(
      rbf_kernel(0.8, c(0.6, 2)),
      periodic_kernel(0.4, 1, period = 1.3)
    ),
    start = sum_kernel(
      rbf_kernel(1, c(1, 1)),
      periodic_kernel(0.5, 1, period = 1.2)
    )
  ),
  list(
    name = "2D ARD RBF, n = 300",
    x = cbind(seq(-2, 2, length.out = 300), cos(1:300)),
    truth = rbf_kernel(1, c(0.8, 1.5)),
    start = rbf_kernel(1, c(1, 1))
  )
)

rows <- list()

for (problem in problems) {
  y <- simulate_gp_data(
    problem$x,
    kernel = problem$truth,
    noise_variance = 0.02,
    seed = 7
  )$observed

  for (gradient in c("analytical", "numerical")) {
    elapsed <- numeric(repeats)

    for (repetition in seq_len(repeats)) {
      elapsed[[repetition]] <- system.time(
        model <- optimize_gp(
          problem$x,
          y,
          kernel = problem$start,
          noise_variance = 0.05,
          gradient = gradient
        )
      )[["elapsed"]]
    }

    optimization <- model$optimization

    rows[[length(rows) + 1L]] <- data.frame(
      problem = problem$name,
      parameters = length(optimization$optimized_parameters),
      gradient = gradient,
      function_evaluations = optimization$function_evaluations,
      gradient_evaluations = optimization$gradient_evaluations,
      model_fits = optimization$model_fits,
      median_seconds = stats::median(elapsed),
      log_marginal_likelihood = optimization$log_marginal_likelihood,
      converged = optimization$converged
    )
  }
}

results <- do.call(rbind, rows)
print(results, row.names = FALSE)
