# Bias and convergence of the heteroscedastic noise estimator
#
# Run manually with the package attached, for example
#   Rscript -e 'library(gaussianprocesses); source("inst/benchmarks/heteroscedastic-noise.R")'
# This file is intentionally excluded from CI. Set `replicates` before
# sourcing to change the number of simulated data sets (default 50).
#
# Data come from gp_simulation_scenario("heteroscedastic_1d"), whose noise
# variance sigma^2(x) = 0.015 + 0.25 x^2 is known. Each data set is fitted
# with the true mean kernel and an RBF noise kernel, allowing up to 200
# iterations, under the three ways of forming the noise GP's targets (see
# .noise_targets() in R/gp-heteroscedastic.R):
#   residual      log squared leave-one-out residuals, corrected only for the
#                 mean of log chi^2_1 (the estimator before this study);
#   standardized  the same, minus log((v_i + s_i^2) / s_i^2) for the
#                 leave-one-out latent variance v_i and current noise s_i^2;
#   subtracted    log(max(e_i^2 - v_i, floor)), corrected for log chi^2_1.
# For each fit the error is log(estimated / true noise variance) at the
# training inputs. The tables give, over replicates, the mean of the median
# error (the bias), the mean root-mean-square error, their Monte Carlo
# standard errors, the share of fits that converged within 200, 6 (the
# former default), 30, and 50 iterations, and the iterations needed.

if (!exists("replicates")) {
  replicates <- 50L
}

rules <- c("residual", "standardized", "subtracted")
sizes <- c(40L, 100L, 200L)
old_option <- getOption("gaussianprocesses.heteroscedastic_target")

fit_replicate <- function(n, seed) {
  simulation <- gp_simulation_scenario("heteroscedastic_1d", n = n, seed = seed)
  model <- withCallingHandlers(
    fit_heteroscedastic_gp(
      simulation$x,
      simulation$observed,
      kernel = simulation$kernel,
      noise_kernel = rbf_kernel(length_scale = 0.4),
      max_iterations = 200
    ),
    gaussianprocesses_convergence_warning = function(warning) {
      invokeRestart("muffleWarning")
    }
  )
  error <- log(model$training_noise_variance / simulation$noise_variance)

  c(
    bias = stats::median(error),
    rmse = sqrt(mean(error^2)),
    iterations = model$iterations,
    converged = model$converged
  )
}

rows <- list()

for (rule in rules) {
  options(gaussianprocesses.heteroscedastic_target = rule)

  for (n in sizes) {
    results <- t(vapply(
      seq_len(replicates),
      function(seed) fit_replicate(n, seed),
      numeric(4)
    ))

    rows[[length(rows) + 1L]] <- data.frame(
      rule = rule,
      n = n,
      bias = mean(results[, "bias"]),
      bias_mcse = stats::sd(results[, "bias"]) / sqrt(replicates),
      rmse = mean(results[, "rmse"]),
      rmse_mcse = stats::sd(results[, "rmse"]) / sqrt(replicates),
      converged = mean(results[, "converged"]),
      # The iterations do not depend on the cap, so the share converged
      # within k iterations follows from the uncapped runs.
      within_6 = mean(results[, "converged"] & results[, "iterations"] <= 6),
      within_30 = mean(results[, "converged"] & results[, "iterations"] <= 30),
      within_50 = mean(results[, "converged"] & results[, "iterations"] <= 50),
      iterations_median = stats::median(results[, "iterations"]),
      iterations_90 = unname(stats::quantile(results[, "iterations"], 0.9)),
      iterations_max = max(results[, "iterations"])
    )
  }
}

options(gaussianprocesses.heteroscedastic_target = old_option)

study <- do.call(rbind, rows)
numeric_columns <- c("bias", "bias_mcse", "rmse", "rmse_mcse")
study[numeric_columns] <- lapply(study[numeric_columns], round, 3)
cat(replicates, "replicates per row\n")
print(study, row.names = FALSE)
