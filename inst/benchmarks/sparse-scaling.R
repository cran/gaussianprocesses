# Time and peak memory of sparse Gaussian-process regression (FITC and VFE)
#
# Run manually after installing the package, for example
#   Rscript inst/benchmarks/sparse-scaling.R
# This file is intentionally excluded from CI and needs the callr package.
# Run it on an otherwise idle machine: timings vary with load.
#
# For n training inputs and m inducing points it measures, each in a fresh R
# process,
#   fit       fit_sparse_gp() at fixed hyperparameters;
#   gradient  one evaluation of the analytical gradient of the objective
#             with respect to the kernel parameters and the noise variance,
#             the cost per optimizer step on top of a fit;
#   predict   marginal prediction at 1000 new inputs.
# Peak memory is R's allocation high-water mark from gc() during the call,
# above what was in use before it, as in prediction-memory.R; BLAS workspace
# is not included. The retained size is that of the fitted model, which
# keeps the m x n cross-covariance. Each measurement is repeated three times
# in fresh processes, and the median time is reported: single timings of
# the largest cases varied by a factor of two on a desktop machine.

run_case <- function(n, m, method, step) {
  library(gaussianprocesses)

  set.seed(1)
  x <- cbind(stats::runif(n, 0, 10), stats::runif(n, 0, 10))
  y <- sin(x[, 1]) * cos(x[, 2]) + stats::rnorm(n, sd = 0.1)
  kernel <- rbf_kernel(variance = 1, length_scale = c(1.5, 1.5))
  fit <- function() {
    fit_sparse_gp(x, y, kernel, noise_variance = 0.01, n_inducing = m,
                  method = method)
  }
  model <- if (!identical(step, "fit")) fit()
  x_new <- cbind(stats::runif(1000, 0, 10), stats::runif(1000, 0, 10))
  compute <- switch(
    step,
    fit = fit,
    gradient = function() {
      gaussianprocesses:::.sparse_gradient(model, include_noise = TRUE)
    },
    predict = function() predict_sparse_gp(model, x_new)
  )

  baseline <- sum(gc(reset = TRUE)[, 2L])
  seconds <- system.time(value <- compute())[["elapsed"]]
  peak_mb <- sum(gc()[, 6L]) - baseline

  list(
    seconds = seconds,
    peak_mb = peak_mb,
    model_mb = if (identical(step, "fit")) {
      as.numeric(utils::object.size(value)) / 2^20
    } else {
      NA_real_
    }
  )
}


rows <- list()

for (n in c(5000L, 10000L, 20000L, 50000L)) {
  for (m in c(100L, 200L, 500L)) {
    for (method in c("fitc", "vfe")) {
      for (step in c("fit", "gradient", "predict")) {
        repeats <- lapply(seq_len(3L), function(repeat_index) {
          callr::r(run_case, list(n = n, m = m, method = method, step = step))
        })
        seconds <- vapply(repeats, function(result) result$seconds,
                          numeric(1))
        rows[[length(rows) + 1L]] <- data.frame(
          n = n,
          m = m,
          method = method,
          step = step,
          seconds = round(stats::median(seconds), 2),
          peak_mb = round(max(vapply(repeats, function(result) result$peak_mb,
                                     numeric(1))), 1),
          model_mb = round(repeats[[1L]]$model_mb, 1)
        )
        message(n, " ", m, " ", method, " ", step, ": ",
                paste(round(seconds, 2), collapse = ", "), " s")
      }
    }
  }
}

results <- do.call(rbind, rows)
print(results, row.names = FALSE)
