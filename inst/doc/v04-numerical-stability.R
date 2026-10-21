## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(gaussianprocesses)

## -----------------------------------------------------------------------------
model <- fit_gp(
  x = c(0, 0),
  y = c(1, 1),
  kernel = rbf_kernel(),
  noise_variance = 0,
  fallback_jitter = 1e-10
)

model$noise_variance
model$jitter
model$cholesky_attempts

## -----------------------------------------------------------------------------
# 40 nearly coincident inputs without noise make C numerically singular, so
# every fit needs jitter.
x <- seq(0, 1, length.out = 40)

fit_in_units <- function(c) {
  fit_gp(
    x,
    c * sin(3 * x),
    kernel = rbf_kernel(variance = c^2, length_scale = 0.5),
    noise_variance = 0
  )
}

unit <- fit_in_units(1)
small <- fit_in_units(1e-4)

# The absolute jitter scales with C; the relative jitter is the same.
c(unit = unit$jitter, small = small$jitter)
c(
  unit = gp_numerical_diagnostics(unit)$relative_jitter,
  small = gp_numerical_diagnostics(small)$relative_jitter
)

# Relative difference of the posterior means after undoing the scaling.
x_new <- c(0.25, 0.5, 1.5)
unit_mean <- predict_gp(unit, x_new)$mean
small_mean <- predict_gp(small, x_new)$mean / 1e-4
max(abs(small_mean - unit_mean)) / max(abs(unit_mean))

## -----------------------------------------------------------------------------
simulation <- gp_simulation_scenario(
  "near_singular",
  n = 30,
  seed = 1729
)

near_singular_model <- fit_gp(
  simulation$x,
  simulation$observed,
  kernel = simulation$kernel,
  noise_variance = simulation$noise_variance
)

gp_numerical_diagnostics(
  near_singular_model
)

## -----------------------------------------------------------------------------
suite <- run_gp_benchmark_suite(
  scenarios = "near_singular",
  n = 20,
  n_inducing = 6,
  repeats = 1,
  seed = 42
)

suite$stability
suite$performance
suite$accuracy

