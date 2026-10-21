# Standard local benchmark suite
#
# This script is intentionally excluded from CI performance gating.

suite <- run_gp_benchmark_suite(
  n = 120,
  n_inducing = 25,
  repeats = 3,
  seed = 1729
)

cat("\nAccuracy\n")
print(suite$accuracy)

cat("\nPerformance\n")
print(suite$performance)

cat("\nNumerical stability\n")
print(suite$stability)

cat("\nExact-GP scaling\n")
print(
  benchmark_exact_gp_scaling(
    sizes = c(50, 100, 200, 400),
    repeats = 3,
    seed = 1729
  )
)
