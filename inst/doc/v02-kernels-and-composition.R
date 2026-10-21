## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(gaussianprocesses)

## -----------------------------------------------------------------------------
x <- seq(-2, 2, length.out = 5)

kernel_rbf(
  x,
  variance = 2,
  length_scale = 0.8
)

specification <- rbf_kernel(
  variance = 2,
  length_scale = 0.8
)

evaluate_kernel(specification, x)

## -----------------------------------------------------------------------------
rbf_kernel(
  variance = 1.5,
  length_scale = 0.7
)

## -----------------------------------------------------------------------------
rational_quadratic_kernel(
  variance = 1,
  length_scale = 0.5,
  alpha = 1.2
)

## -----------------------------------------------------------------------------
periodic_kernel(
  variance = 1,
  length_scale = 0.8,
  period = 12
)

## -----------------------------------------------------------------------------
trend <- linear_kernel(variance = 0.05)

seasonal <- periodic_kernel(
  variance = 1,
  length_scale = 0.8,
  period = 12
)

local <- matern32_kernel(
  variance = 0.3,
  length_scale = 2
)

combined <- sum_kernel(
  trend,
  product_kernel(seasonal, local)
)

combined

## -----------------------------------------------------------------------------
time <- 0:6

K_combined <- evaluate_kernel(combined, time)

K_manual <-
  evaluate_kernel(trend, time) +
  evaluate_kernel(seasonal, time) *
  evaluate_kernel(local, time)

max(abs(K_combined - K_manual))

## -----------------------------------------------------------------------------
kernel_parameters(
  combined,
  flatten = TRUE
)

## -----------------------------------------------------------------------------
updated <- update_kernel_parameters(
  combined,
  c("kernel2.kernel1.period" = 24)
)

kernel_parameters(
  updated,
  flatten = TRUE
)

