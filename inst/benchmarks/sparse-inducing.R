# Inducing-point selection and optimization for sparse Gaussian processes
#
# Run manually after installing the package, for example
#   Rscript inst/benchmarks/sparse-inducing.R
# It takes several minutes and is excluded from CI.
#
# 1. Gap closed by optimizing Z. For each simulation scenario of
#    gp_simulation_scenario() with 400 observations, the exact log marginal
#    likelihood is maximized, and the VFE bound is maximized with m = 5, 10,
#    and 20 inducing points that are fixed at farthest-point selection, fixed
#    at greedy variance selection, or optimized jointly from the greedy
#    selection. The gap is the exact maximum minus the bound's maximum, and
#    the fraction closed is one minus the gap with optimized inducing points
#    over the gap with fixed ones.
# 2. FITC and VFE with free inducing points. On 20 replicates of the design
#    of sparse-vfe-fitc.R (400 inputs on [-4, 4], sin(2x), noise variance
#    0.1), both are optimized with m = 10 inducing points, fixed or free, and
#    the noise estimate and the objective minus the exact log marginal
#    likelihood at the estimate are recorded.
# 3. Greedy variance selection costs O(n m^2): the median of three timings
#    for n from 5 000 to 40 000 and m from 50 to 800, with a length scale
#    short enough that 800 points do not explain the prior variance (with a
#    longer one the selection stops early and the timings mean nothing).

library(gaussianprocesses)

quietly <- function(expression) {
  withCallingHandlers(
    expression,
    gaussianprocesses_inducing_collapse_warning = function(condition) {
      invokeRestart("muffleWarning")
    }
  )
}

# 1. Gap closed by optimizing Z on the simulation scenarios. -----------------
scenarios <- c("smooth_1d", "near_singular", "ard_2d", "heteroscedastic_1d")
gap_rows <- list()

for (scenario in scenarios) {
  data <- gp_simulation_scenario(scenario, n = 400, seed = 1)
  # Every parameter starts at 1, away from the simulated values.
  start_values <- kernel_parameters(data$kernel, flatten = TRUE)
  start_values[] <- 1
  start_kernel <- update_kernel_parameters(data$kernel, start_values)
  exact <- optimize_gp(data$x, data$observed, start_kernel,
                       noise_variance = 0.1, n_starts = 3)
  exact_maximum <- log_marginal_likelihood(exact)

  for (m in c(5L, 10L, 20L)) {
    bound <- function(selection, optimize_inducing) {
      model <- quietly(optimize_sparse_gp(
        data$x, data$observed, start_kernel, noise_variance = 0.1,
        n_inducing = m, selection = selection, method = "vfe",
        optimize_inducing = optimize_inducing, n_starts = 3
      ))
      model$log_marginal_likelihood
    }
    farthest <- bound("farthest", FALSE)
    variance <- bound("variance", FALSE)
    optimized <- bound("variance", TRUE)

    gap_rows[[length(gap_rows) + 1L]] <- data.frame(
      scenario = scenario,
      m = m,
      gap_farthest = exact_maximum - farthest,
      gap_variance = exact_maximum - variance,
      gap_optimized = exact_maximum - optimized,
      closed_from_farthest = 1 - (exact_maximum - optimized) /
        (exact_maximum - farthest),
      closed_from_variance = 1 - (exact_maximum - optimized) /
        (exact_maximum - variance)
    )
    message(scenario, " m = ", m)
  }
}

gaps <- do.call(rbind, gap_rows)
print(gaps, digits = 3)

# 2. FITC and VFE with free inducing points. --------------------------------
free_rows <- list()

for (replicate in seq_len(20L)) {
  set.seed(1000L + replicate)
  x <- stats::runif(400, -4, 4)
  y <- sin(2 * x) + stats::rnorm(400, sd = sqrt(0.1))
  exact <- optimize_gp(x, y, rbf_kernel(), noise_variance = 0.5, n_starts = 3)

  for (method in c("fitc", "vfe")) {
    for (free in c(FALSE, TRUE)) {
      model <- quietly(optimize_sparse_gp(
        x, y, rbf_kernel(), noise_variance = 0.5, n_inducing = 10,
        selection = "variance", method = method, optimize_inducing = free,
        n_starts = 3
      ))
      exact_at <- log_marginal_likelihood(
        fit_gp(x, y, model$kernel, noise_variance = model$noise_variance)
      )
      free_rows[[length(free_rows) + 1L]] <- data.frame(
        replicate = replicate,
        method = method,
        inducing = if (free) "optimized" else "fixed",
        noise_vs_exact = model$noise_variance / exact$noise_variance,
        objective_gap = model$log_marginal_likelihood - exact_at,
        trace_term = model$trace_term,
        collapsed_pairs = if (is.null(model$inducing_collapse)) {
          0L
        } else {
          nrow(model$inducing_collapse)
        }
      )
    }
  }

  message("replicate ", replicate)
}

free <- do.call(rbind, free_rows)
free_summary <- do.call(rbind, lapply(
  split(free, list(free$method, free$inducing), drop = TRUE),
  function(group) {
    data.frame(
      method = group$method[1],
      inducing = group$inducing[1],
      noise_vs_exact = mean(group$noise_vs_exact),
      noise_vs_exact_se = stats::sd(group$noise_vs_exact) / sqrt(nrow(group)),
      objective_gap = mean(group$objective_gap),
      objective_gap_se = stats::sd(group$objective_gap) / sqrt(nrow(group)),
      above_exact = sum(group$objective_gap > 1e-8),
      trace_term = stats::median(group$trace_term),
      replicates_with_collapse = sum(group$collapsed_pairs > 0)
    )
  }
))
print(free_summary, digits = 3)

# 3. Cost of greedy variance selection. -------------------------------------
timing_rows <- list()
kernel <- rbf_kernel(variance = 1, length_scale = c(0.05, 0.05))

for (n in c(5000L, 10000L, 20000L, 40000L)) {
  set.seed(n)
  x <- cbind(stats::runif(n), stats::runif(n))

  for (m in c(50L, 100L, 200L, 400L, 800L)) {
    seconds <- vapply(seq_len(3L), function(repeat_index) {
      system.time(
        select_inducing_points(x, m, method = "variance", kernel = kernel)
      )[["elapsed"]]
    }, numeric(1))
    timing_rows[[length(timing_rows) + 1L]] <- data.frame(
      n = n, m = m, seconds = stats::median(seconds)
    )
  }
}

timings <- do.call(rbind, timing_rows)
print(timings)

# Log-log slopes of time in n (at the largest m) and in m (at the largest n).
slope <- function(values, sizes) unname(stats::coef(stats::lm(log(values) ~ log(sizes)))[2])
largest_m <- timings[timings$m == max(timings$m), ]
largest_n <- timings[timings$n == max(timings$n), ]
cat("slope in n at m =", max(timings$m), ":",
    round(slope(largest_m$seconds, largest_m$n), 2), "\n")
cat("slope in m at n =", max(timings$n), ":",
    round(slope(largest_n$seconds, largest_n$m), 2), "\n")
