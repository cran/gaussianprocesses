# Hyperparameter learning with FITC and VFE against the exact model
#
# Run manually after installing the package, for example
#   Rscript inst/benchmarks/sparse-vfe-fitc.R
# It takes a few minutes and is excluded from CI.
#
# Bauer, van der Wilk, and Rasmussen (2016) report that FITC's objective is
# not a bound on the marginal likelihood and can exceed it, that FITC tends
# to underestimate the noise variance, and that VFE overestimates the noise
# when the inducing points cover the data poorly but recovers the exact
# hyperparameters as inducing points are added. This script checks those
# claims with fixed inducing points (selected by farthest-point selection,
# not optimized).
#
# Each of 20 replicates draws 400 training and 500 test inputs uniformly on
# [-4, 4], responses sin(2x) plus noise of variance 0.1, and
# fits an RBF kernel and the noise variance by maximizing the exact log
# marginal likelihood, the FITC objective, and the VFE bound, with m = 6,
# 10, 20, and 40 inducing points and the same three starts. For each fit it
# records the noise estimate relative to the true noise variance, the
# objective minus the exact log marginal likelihood at the same
# hyperparameters, the exact log marginal likelihood at those
# hyperparameters minus the best exact value found in the replicate, the
# trace term, and the mean test log predictive density.
#
# The signal has a single length scale, so that the exact likelihood has a
# single optimum in these data; with a second, faster component the exact
# optimizer often stopped at a local optimum, which confounded the
# comparison. As a safeguard, the reference is the largest exact log
# marginal likelihood at any of the replicate's estimates, exact or sparse.

library(gaussianprocesses)

n_replicates <- 20L
inducing_counts <- c(6L, 10L, 20L, 40L)
noise <- 0.1
signal <- function(x) sin(2 * x)

rows <- list()

for (replicate in seq_len(n_replicates)) {
  set.seed(1000L + replicate)
  x <- stats::runif(400, -4, 4)
  y <- signal(x) + stats::rnorm(400, sd = sqrt(noise))
  x_test <- stats::runif(500, -4, 4)
  y_test <- signal(x_test) + stats::rnorm(500, sd = sqrt(noise))

  exact <- optimize_gp(x, y, rbf_kernel(), noise_variance = 0.5, n_starts = 3)
  exact_prediction <- predict_gp(exact, x_test)
  replicate_rows <- list(data.frame(
    replicate = replicate, method = "exact", m = NA_integer_,
    noise_ratio = exact$noise_variance / noise, objective_gap = 0,
    exact_value = log_marginal_likelihood(exact), trace_term = 0,
    test_log_density = mean(stats::dnorm(
      y_test, exact_prediction$mean, exact_prediction$observation_sd,
      log = TRUE
    ))
  ))

  for (m in inducing_counts) {
    for (method in c("fitc", "vfe")) {
      model <- optimize_sparse_gp(x, y, rbf_kernel(), noise_variance = 0.5,
                                  n_inducing = m, method = method,
                                  n_starts = 3)
      exact_at <- log_marginal_likelihood(
        fit_gp(x, y, model$kernel, noise_variance = model$noise_variance)
      )
      prediction <- predict_sparse_gp(model, x_test)
      replicate_rows[[length(replicate_rows) + 1L]] <- data.frame(
        replicate = replicate, method = method, m = m,
        noise_ratio = model$noise_variance / noise,
        objective_gap = model$log_marginal_likelihood - exact_at,
        exact_value = exact_at,
        trace_term = model$trace_term,
        test_log_density = mean(stats::dnorm(
          y_test, prediction$mean, prediction$observation_sd, log = TRUE
        ))
      )
    }
  }

  replicate_rows <- do.call(rbind, replicate_rows)
  replicate_rows$exact_gap <- replicate_rows$exact_value -
    max(replicate_rows$exact_value)
  # Paired with the exact fit of the same data.
  replicate_rows$noise_vs_exact <- replicate_rows$noise_ratio /
    replicate_rows$noise_ratio[1L]
  replicate_rows$density_vs_exact <- replicate_rows$test_log_density -
    replicate_rows$test_log_density[1L]
  rows[[length(rows) + 1L]] <- replicate_rows
  message("replicate ", replicate)
}

results <- do.call(rbind, rows)

# Means over replicates with Monte Carlo standard errors.
summary <- do.call(rbind, lapply(
  split(results, list(results$method, results$m), drop = TRUE),
  function(group) {
    columns <- c("noise_ratio", "noise_vs_exact", "objective_gap",
                 "exact_gap", "trace_term", "test_log_density",
                 "density_vs_exact")
    row <- data.frame(method = group$method[1], m = group$m[1],
                      positive_gaps = sum(group$objective_gap > 1e-8),
                      largest_gap = max(group$objective_gap))
    for (column in columns) {
      row[[column]] <- mean(group[[column]])
      row[[paste0(column, "_se")]] <- stats::sd(group[[column]]) /
        sqrt(nrow(group))
    }
    row
  }
))
summary <- summary[order(summary$method, summary$m), ]
rownames(summary) <- NULL
print(summary, digits = 3)

# The exact model's own results.
exact_rows <- results[results$method == "exact", ]
for (column in c("noise_ratio", "exact_gap", "test_log_density")) {
  cat("exact", column, ":", round(mean(exact_rows[[column]]), 4), "(",
      round(stats::sd(exact_rows[[column]]) / sqrt(n_replicates), 4), ")\n")
}
