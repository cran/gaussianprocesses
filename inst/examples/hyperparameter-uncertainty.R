# Do reported standard errors match the spread of the estimates?
#
# Simulate 200 data sets from a known Gaussian process, estimate its
# hyperparameters on each with optimize_gp(), and compare the standard errors
# that gp_hyperparameter_uncertainty() reports with the standard deviation of
# the estimates across data sets. Everything is on the optimizer's log scale.

library(gaussianprocesses)

truth <- rbf_kernel(variance = 1, length_scale = 0.8)
true_noise <- 0.04
x <- seq(0, 10, length.out = 60)
true_values <- log(c(variance = 1, length_scale = 0.8, noise_variance = 0.04))

replicates <- lapply(
  seq_len(200),
  function(seed) {
    simulated <- simulate_gp_data(
      x,
      truth,
      noise_variance = true_noise,
      seed = seed
    )
    model <- optimize_gp(
      x,
      simulated$observed,
      kernel = truth,
      noise_variance = true_noise,
      n_starts = 1
    )
    uncertainty <- gp_hyperparameter_uncertainty(model)

    list(
      estimate = uncertainty$parameters$coordinate_estimate,
      standard_error = uncertainty$parameters$standard_error,
      poorly_identified = uncertainty$poorly_identified
    )
  }
)

estimates <- do.call(rbind, lapply(replicates, `[[`, "estimate"))
standard_errors <- do.call(rbind, lapply(replicates, `[[`, "standard_error"))
colnames(estimates) <- colnames(standard_errors) <- names(true_values)

z <- stats::qnorm(0.975)
covered <- abs(sweep(estimates, 2, true_values)) <= z * standard_errors

study <- data.frame(
  parameter = names(true_values),
  spread_of_estimates = apply(estimates, 2, stats::sd),
  # The scaled median absolute deviation equals the standard deviation for
  # normal estimates but ignores their tails.
  robust_spread = apply(estimates, 2, stats::mad),
  median_standard_error = apply(standard_errors, 2, stats::median),
  coverage_95 = colMeans(covered),
  row.names = NULL
)
study

# Monte Carlo standard error of the coverage of 200 intervals.
sqrt(0.95 * 0.05 / 200)
sum(vapply(replicates, `[[`, logical(1), "poorly_identified"))
