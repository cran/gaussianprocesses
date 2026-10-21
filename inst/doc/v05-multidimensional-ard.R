## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)
library(gaussianprocesses)

## -----------------------------------------------------------------------------
x <- rbind(
  c(0, 0),
  c(1, 0),
  c(0, 1)
)

isotropic <- rbf_kernel(
  variance = 1,
  length_scale = 1
)

ard <- rbf_kernel(
  variance = 1,
  length_scale = c(0.4, 2)
)

evaluate_kernel(isotropic, x)
evaluate_kernel(ard, x)

## -----------------------------------------------------------------------------
kernel_parameters(
  ard,
  flatten = TRUE
)

## -----------------------------------------------------------------------------
ard_updated <- update_kernel_parameters(
  ard,
  c("length_scale[2]" = 3)
)

kernel_parameters(
  ard_updated,
  flatten = TRUE
)

## -----------------------------------------------------------------------------
simulation <- gp_simulation_scenario(
  "ard_2d",
  n = 40,
  seed = 1729
)

head(simulation$x)
simulation$kernel

## -----------------------------------------------------------------------------
model <- fit_gp(
  simulation$x,
  simulation$observed,
  kernel = simulation$kernel,
  noise_variance = simulation$noise_variance
)

prediction <- predict_gp(
  model,
  simulation$x,
  observation_noise_variance = 0
)

sqrt(
  mean(
    (prediction$mean - simulation$latent)^2
  )
)

## ----eval=FALSE---------------------------------------------------------------
# estimated <- optimize_gp(
#   simulation$x,
#   simulation$observed,
#   kernel = rbf_kernel(
#     variance = 1,
#     length_scale = c(1, 1)
#   ),
#   noise_variance = 0.04,
#   optimize_noise = FALSE
# )
# 
# kernel_parameters(
#   estimated$kernel,
#   flatten = TRUE
# )

## -----------------------------------------------------------------------------
additive <- sum_kernel(
  select_dimensions(rbf_kernel(), 1),
  select_dimensions(rbf_kernel(length_scale = c(1, 1)), 2:3)
)

additive

## -----------------------------------------------------------------------------
names(kernel_parameters(additive, flatten = TRUE))

## -----------------------------------------------------------------------------
set.seed(2025)

f1 <- function(x1) sin(2 * x1)
f2 <- function(x2, x3) 1.5 * exp(-(x2^2 + x3^2) / 2)
latent <- function(x) f1(x[, 1]) + f2(x[, 2], x[, 3])

x_train <- matrix(runif(3 * 100, -2, 2), ncol = 3)
y_train <- latent(x_train) + rnorm(100, sd = 0.1)
x_test <- matrix(runif(3 * 500, -2, 2), ncol = 3)

## -----------------------------------------------------------------------------
additive_model <- optimize_gp(
  x_train,
  y_train,
  kernel = additive,
  noise_variance = 0.05
)
full_model <- optimize_gp(
  x_train,
  y_train,
  kernel = rbf_kernel(length_scale = c(1, 1, 1)),
  noise_variance = 0.05
)

test_rmse <- function(model) {
  sqrt(mean((predict_gp(model, x_test)$mean - latent(x_test))^2))
}

data.frame(
  kernel = c("additive", "ARD RBF on all columns"),
  log_marginal_likelihood = c(
    log_marginal_likelihood(additive_model),
    log_marginal_likelihood(full_model)
  ),
  test_rmse = c(test_rmse(additive_model), test_rmse(full_model))
)

## -----------------------------------------------------------------------------
component_mean <- function(model, index, x) {
  component <- model$kernel$children[[index]]
  drop(evaluate_kernel(component, x, model$x) %*% model$alpha)
}

grid <- seq(-2, 2, length.out = 101)
# Component 1 ignores the second and third columns.
estimate <- component_mean(additive_model, 1, cbind(grid, 0, 0))

## ----fig.alt = "The centred posterior mean of the first additive component as a solid line over x1 from -2 to 2, and the centred true function sin(2 x1) as a dashed line. The two curves nearly coincide."----
centre <- function(values) values - mean(values)

plot(
  grid,
  centre(f1(grid)),
  type = "l",
  lty = 2,
  xlab = expression(x[1]),
  ylab = "centred component 1"
)
lines(grid, centre(estimate))
legend(
  "topleft",
  legend = c("posterior mean", "true function"),
  lty = c(1, 2),
  bty = "n"
)

max(abs(centre(estimate) - centre(f1(grid))))

