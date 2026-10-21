# Sparse GP benchmark
#
# Run manually. This file is intentionally not part of CI because exact GP
# scaling can make larger benchmark cases expensive.

set.seed(1729)

n <- 500
x <- sort(stats::runif(n, -4, 4))
latent <- sin(x) + 0.2 * cos(3 * x)
y <- latent + stats::rnorm(n, sd = 0.15)

kernel <- matern32_kernel(
  variance = 1,
  length_scale = 0.8
)

result <- benchmark_sparse_gp(
  x,
  y,
  kernel = kernel,
  n_inducing = 40,
  prediction_x = seq(-4, 4, length.out = 300),
  noise_variance = 0.15^2,
  selection = "quantile"
)

print(result$comparison)
