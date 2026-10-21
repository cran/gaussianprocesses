# Time and peak memory of marginal exact-GP prediction
#
# Run manually after installing the package, for example
#   Rscript inst/benchmarks/prediction-memory.R
# This file is intentionally excluded from CI and needs the callr package.
# Run it on an otherwise idle machine: timings vary with load.
#
# Each measurement runs in a fresh R process. Peak memory is R's allocation
# high-water mark from gc() during the call, above what was in use before it.
# It covers vectors and matrices allocated by R, which dominate here; BLAS
# workspace is not included. A fresh process matters: an earlier large
# computation raises R's garbage-collection threshold, and later peaks would
# then include garbage that a fresh session collects.
#
# The "previous m x m" path reproduces what predict_gp() computed for
# marginal predictions before the blocked path: it formed the m x m prior and
# posterior covariances and took their diagonals. It is run only for m up to
# 8000, because at m = 50000 those matrices alone would need about 40 GB.

run_case <- function(m, path) {
  library(gaussianprocesses)

  previous_marginal_prediction <- function(model, x) {
    x <- matrix(x, ncol = 1L)
    cross_covariance <- evaluate_kernel(model$kernel, model$x, x)
    prior_covariance <- evaluate_kernel(model$kernel, x)
    posterior_mean <- evaluate_mean(model$mean, x) +
      drop(crossprod(cross_covariance, model$alpha))
    projected <- forwardsolve(t(model$cholesky), cross_covariance)
    latent_covariance <- prior_covariance - crossprod(projected)
    latent_covariance <- 0.5 * (latent_covariance + t(latent_covariance))

    list(
      mean = posterior_mean,
      latent_variance = pmax(diag(latent_covariance), 0)
    )
  }

  x <- seq(0, 10, length.out = 500L)
  model <- fit_gp(
    x,
    sin(x) + 0.3 * cos(3 * x),
    kernel = rbf_kernel(variance = 1, length_scale = 0.8),
    noise_variance = 0.01
  )
  x_new <- seq(-1, 11, length.out = m)
  compute <- switch(
    path,
    blocked = function() predict_gp(model, x_new)[c("mean", "latent_variance")],
    previous = function() previous_marginal_prediction(model, x_new)
  )

  baseline <- sum(gc(reset = TRUE)[, 2L])
  seconds <- system.time(value <- compute())[["elapsed"]]
  peak_mb <- sum(gc()[, 6L]) - baseline

  list(value = value, seconds = seconds, peak_mb = peak_mb)
}


relative_difference <- function(actual, expected) {
  max(abs(actual - expected)) / max(abs(expected))
}


n <- 500L
block_size <- gaussianprocesses:::.prediction_block_size(
  NULL,
  n_training = n,
  n_predictions = .Machine$integer.max
)
rows <- list()

for (m in c(1000L, 8000L, 50000L)) {
  blocked <- callr::r(run_case, list(m = m, path = "blocked"))
  row <- data.frame(
    m = m,
    path = "blocked marginal",
    block_size = min(block_size, m),
    seconds = blocked$seconds,
    peak_mb = blocked$peak_mb,
    mean_difference = NA_real_,
    variance_difference = NA_real_
  )

  if (m <= 8000L) {
    previous <- callr::r(run_case, list(m = m, path = "previous"))
    row$mean_difference <- relative_difference(
      blocked$value$mean,
      previous$value$mean
    )
    row$variance_difference <- relative_difference(
      blocked$value$latent_variance,
      previous$value$latent_variance
    )
    row <- rbind(
      row,
      data.frame(
        m = m,
        path = "previous m x m",
        block_size = NA_integer_,
        seconds = previous$seconds,
        peak_mb = previous$peak_mb,
        mean_difference = NA_real_,
        variance_difference = NA_real_
      )
    )
  }

  rows[[length(rows) + 1L]] <- row
}

results <- do.call(rbind, rows)
results$seconds <- round(results$seconds, 2)
results$peak_mb <- round(results$peak_mb, 1)
results$mean_difference <- signif(results$mean_difference, 2)
results$variance_difference <- signif(results$variance_difference, 2)
cat(n, "training inputs; differences are relative to the previous path\n")
print(results, row.names = FALSE)
