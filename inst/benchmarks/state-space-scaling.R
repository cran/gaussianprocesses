# Time of state-space (Kalman) inference against exact inference
#
# Run manually after installing the package, for example
#   Rscript inst/benchmarks/state-space-scaling.R
# This file is intentionally excluded from CI and needs the callr package.
# Run it on an otherwise idle machine: timings vary with load.
#
# For n irregularly spaced training times it measures, each in a fresh R
# process,
#   likelihood  one Kalman filter pass (the log marginal likelihood only),
#               the cost of one objective value in optimize_time_series_gp();
#   fit         fit_time_series_gp(method = "state_space"): filter and
#               Rauch-Tung-Striebel smoother;
#   forecast    forecast_gp() at 1000 new times, which filters and smooths
#               the training and new times together;
#   exact       fit_time_series_gp(method = "exact"), for n up to 5000.
# The kernels are Matern-5/2 (a state of dimension 3) and the sum of a
# Matern-5/2 and a Matern-1/2 kernel (dimension 4). Each measurement is
# repeated three times in fresh processes, after an untimed fit on 50
# points, and the median is reported. The log-log slope of time against n
# is the least-squares slope over all n; linear scaling gives a slope near
# 1, and Cholesky-based exact inference a slope near 3 once the O(n^3)
# factorization dominates.

run_case <- function(n, kernel_name, step) {
  library(gaussianprocesses)

  # Exponential gaps give irregular, strictly increasing times; sorted
  # uniform draws can tie, because runif() has 32-bit resolution.
  set.seed(1)
  time <- cumsum(stats::rexp(n, rate = 10))
  y <- sin(time) + stats::rnorm(n, sd = 0.3)
  kernel <- switch(
    kernel_name,
    matern52 = matern52_kernel(length_scale = 2),
    sum = sum_kernel(matern52_kernel(length_scale = 2),
                     matern12_kernel(variance = 0.1, length_scale = 0.5))
  )
  future <- seq(0, max(time) + 5, length.out = 1000)
  model <- if (identical(step, "forecast")) {
    fit_time_series_gp(time, y, kernel, noise_variance = 0.09,
                       method = "state_space")
  }
  # An untimed fit on 50 points first, so that lazy loading and the kernel
  # registry's initialization are not part of the timed call.
  invisible(fit_time_series_gp(time[1:50], y[1:50], kernel,
                               noise_variance = 0.09,
                               method = if (identical(step, "exact")) {
                                 "exact"
                               } else {
                                 "state_space"
                               }))
  compute <- switch(
    step,
    likelihood = function() {
      gaussianprocesses:::.fit_state_space_model(time, y, kernel, 0.09,
                                                 smooth = FALSE)
    },
    fit = function() {
      fit_time_series_gp(time, y, kernel, noise_variance = 0.09,
                         method = "state_space")
    },
    forecast = function() forecast_gp(model, future),
    exact = function() {
      fit_time_series_gp(time, y, kernel, noise_variance = 0.09)
    }
  )

  system.time(compute())[["elapsed"]]
}


sizes <- c(1000L, 2000L, 5000L, 10000L, 20000L, 50000L, 100000L)
rows <- list()

for (kernel_name in c("matern52", "sum")) {
  for (step in c("likelihood", "fit", "forecast", "exact")) {
    for (n in sizes) {
      if (identical(step, "exact") && n > 5000L) {
        next
      }

      seconds <- vapply(seq_len(3L), function(repeat_index) {
        callr::r(run_case, list(n = n, kernel_name = kernel_name,
                                step = step))
      }, numeric(1))
      rows[[length(rows) + 1L]] <- data.frame(
        kernel = kernel_name,
        step = step,
        n = n,
        seconds = stats::median(seconds)
      )
      message(kernel_name, " ", step, " ", n, ": ",
              paste(round(seconds, 3), collapse = ", "), " s")
    }
  }
}

results <- do.call(rbind, rows)
print(results, row.names = FALSE)

slopes <- do.call(rbind, lapply(
  split(results, list(results$kernel, results$step), drop = TRUE),
  function(group) {
    data.frame(
      kernel = group$kernel[[1L]],
      step = group$step[[1L]],
      n_range = paste(range(group$n), collapse = "-"),
      slope = round(unname(stats::coef(
        stats::lm(log(seconds) ~ log(n), data = group)
      )[[2L]]), 2),
      microseconds_per_point = round(1e6 * max(group$seconds / group$n), 1)
    )
  }
))
print(slopes, row.names = FALSE)
