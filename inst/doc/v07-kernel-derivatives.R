## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(gaussianprocesses)

## -----------------------------------------------------------------------------
kernel <- scale_kernel(
  sum_kernel(
    product_kernel(
      rbf_kernel(variance = 0.9, length_scale = c(0.5, 1.5)),
      periodic_kernel(variance = 1.2, length_scale = 0.7, period = 1.6)
    ),
    white_noise_kernel(variance = 0.05)
  ),
  scale = 0.7
)

x <- cbind(
  seq(-1, 1, length.out = 7),
  cos(1:7)
)

gradient <- kernel_gradient(kernel, x)

names(gradient)

## -----------------------------------------------------------------------------
all.equal(gradient[["scale"]], evaluate_kernel(kernel, x))

## -----------------------------------------------------------------------------
x_train <- seq(-2, 2, length.out = 20)
model <- fit_gp(
  x_train,
  sin(2 * x_train),
  kernel = rbf_kernel(variance = 1, length_scale = 0.6),
  noise_variance = 0.02
)

log_marginal_likelihood_gradient(model)

## -----------------------------------------------------------------------------
finite_difference <- function(kernel, x, step) {
  parameters <- kernel_parameters(kernel, flatten = TRUE)
  derivatives <- lapply(seq_along(parameters), function(j) {
    covariance_at <- function(direction) {
      shifted <- parameters[j] * exp(direction * step)
      evaluate_kernel(update_kernel_parameters(kernel, shifted), x)
    }
    (covariance_at(1) - covariance_at(-1)) / (2 * step)
  })
  names(derivatives) <- names(parameters)
  derivatives
}

discrepancy <- function(step) {
  numerical <- finite_difference(kernel, x, step)
  max(mapply(function(a, n) max(abs(a - n)), gradient, numerical))
}

data.frame(
  step = c(1e-3, 1e-4, 1e-5),
  max_abs_discrepancy = sapply(c(1e-3, 1e-4, 1e-5), discrepancy)
)

## -----------------------------------------------------------------------------
kernel <- product_kernel(
  matern52_kernel(variance = 1.2, length_scale = c(0.6, 1.4)),
  linear_kernel(variance = 0.5)
)
x <- cbind(c(-0.5, 0, 0.8), c(1, 0.3, -0.2))

gradient <- kernel_input_gradient(kernel, x)

step <- 1e-6
shift <- cbind(step, rep(0, nrow(x)))
numerical <- (evaluate_kernel(kernel, x + shift, x) -
  evaluate_kernel(kernel, x - shift, x)) / (2 * step)

max(abs(gradient[[1]] - numerical))

## -----------------------------------------------------------------------------
x_train <- seq(0, 2 * pi, length.out = 30)
model <- fit_gp(
  x_train,
  sin(x_train),
  kernel = rbf_kernel(length_scale = 1.2),
  noise_variance = 1e-4
)

x_new <- seq(0, 2 * pi, length.out = 7)
gradient <- predict_gradient_gp(model, x_new)

data.frame(
  x = round(x_new, 3),
  estimate = gradient$mean[, 1],
  sd = gradient$sd[, 1],
  derivative = cos(x_new)
)

## -----------------------------------------------------------------------------
rough <- fit_gp(x_train, sin(x_train), matern12_kernel(), noise_variance = 1e-4)

tryCatch(
  predict_gradient_gp(rough, 1),
  gaussianprocesses_smoothness_error = conditionMessage
)

## -----------------------------------------------------------------------------
x_observed <- c(-2, -1, 0, 1, 2)
kernel <- rbf_kernel(length_scale = 0.8)

values_only <- fit_gp(
  x_observed,
  sin(x_observed),
  kernel = kernel,
  noise_variance = 1e-4
)
with_gradients <- fit_gp(
  c(x_observed, x_observed),
  c(sin(x_observed), cos(x_observed)),
  kernel = kernel,
  noise_variance = 1e-4,
  derivative = rep(c(0, 1), each = 5)
)
with_gradients

## -----------------------------------------------------------------------------
joint <- crossprod(with_gradients$cholesky)
min(eigen(joint, symmetric = TRUE, only.values = TRUE)$values)

## ----fig.alt = "Posterior standard deviation of the latent function over x from -2.5 to 2.5. With values only, it rises to about 0.23 between the five observations; with values and gradients it stays below about 0.02 between them and is lower everywhere, also beyond the data.", fig.width = 7, fig.height = 4----
grid <- seq(-2.5, 2.5, length.out = 201)
before <- predict_gp(values_only, grid)
after <- predict_gp(with_gradients, grid)

plot(
  grid,
  before$latent_sd,
  type = "l",
  lty = 2,
  ylim = c(0, max(before$latent_sd)),
  xlab = "x",
  ylab = "posterior standard deviation"
)
lines(grid, after$latent_sd)
points(x_observed, rep(0, 5), pch = 19)
legend(
  "top",
  legend = c("values", "values and gradients"),
  lty = c(2, 1),
  bty = "n"
)

## -----------------------------------------------------------------------------
between <- abs(grid - round(grid)) > 0.4 & abs(grid) < 2
c(
  values = max(before$latent_sd[between]),
  values_and_gradients = max(after$latent_sd[between])
)

## -----------------------------------------------------------------------------
predict_gradient_gp(with_gradients, x_observed)$mean[, 1] - cos(x_observed)

## -----------------------------------------------------------------------------
gp_numerical_diagnostics(with_gradients)$observation_types

