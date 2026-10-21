## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(gaussianprocesses)

## -----------------------------------------------------------------------------
x <- seq(-1.5, 1.5, length.out = 12)
y <- sin(2 * x) + 0.1 * x

initial_kernel <- rbf_kernel(
  variance = 0.5,
  length_scale = 1.5
)

initial_model <- fit_gp(
  x,
  y,
  kernel = initial_kernel,
  noise_variance = 0.05
)

log_marginal_likelihood(initial_model)

## -----------------------------------------------------------------------------
optimized <- optimize_gp(
  x,
  y,
  kernel = initial_kernel,
  noise_variance = 0.05,
  optimize_noise = TRUE,
  n_starts = 2,
  control = list(maxit = 30)
)

kernel_parameters(
  optimized$kernel,
  flatten = TRUE
)

optimized$noise_variance
optimized$optimization$log_marginal_likelihood

## -----------------------------------------------------------------------------
optimized$optimization[c(
  "method",
  "converged",
  "best_start",
  "counts"
)]

## -----------------------------------------------------------------------------
x_new <- seq(-2, 2, length.out = 25)

prediction <- predict_gp(
  optimized,
  x_new,
  interval_level = 0.9
)

head(
  data.frame(
    x = x_new,
    mean = prediction$mean,
    latent_variance = prediction$latent_variance,
    observation_variance = prediction$observation_variance
  )
)

## -----------------------------------------------------------------------------
loo <- loo_gp(optimized)

head(loo)

## -----------------------------------------------------------------------------
gp_calibration_diagnostics(
  optimized,
  interval_level = 0.9
)

## -----------------------------------------------------------------------------
set.seed(26)
x_trend <- sort(runif(30, 0, 10))
truth <- function(x) 2 + 0.8 * x + sin(1.3 * x)
y_trend <- truth(x_trend) + rnorm(30, sd = 0.2)

fit_with_mean <- function(mean) {
  optimize_gp(
    x_trend,
    y_trend,
    kernel = rbf_kernel(variance = 1, length_scale = 1),
    noise_variance = 0.05,
    mean = mean
  )
}

trend_models <- list(
  "zero mean" = fit_with_mean(zero_mean()),
  "estimated, ML" = fit_with_mean(linear_mean()),
  "estimated, REML" = fit_with_mean(
    linear_mean(estimate_coefficients("reml"))
  ),
  "vague prior" = fit_with_mean(linear_mean(coefficient_prior()))
)

data.frame(
  likelihood = sapply(trend_models, function(model) {
    model$optimization$likelihood
  }),
  variance = sapply(trend_models, function(model) {
    model$kernel$parameters$variance
  }),
  length_scale = sapply(trend_models, function(model) {
    model$kernel$parameters$length_scale
  })
)

## -----------------------------------------------------------------------------
estimated <- trend_models[["estimated, REML"]]

cbind(
  estimate = coef(estimated),
  standard_error = sqrt(diag(vcov(estimated)))
)

## -----------------------------------------------------------------------------
t(sapply(trend_models, function(model) {
  prediction <- predict_gp(model, 15)
  c(mean = prediction$mean, latent_sd = prediction$latent_sd)
}))

## ----trend-bands, fig.alt = "Latent 95% intervals over x from 0 to 15 for the REML plug-in model (dashed lines) and the vague-prior model (solid lines), with the data as points on 0 to 10 and the true function as a grey line. The two intervals coincide over the data and the vague-prior interval becomes wider beyond x = 10.", fig.width = 7, fig.height = 4----
x_grid <- seq(0, 15, length.out = 200)
plug_in <- predict_gp(trend_models[["estimated, REML"]], x_grid)
marginal <- predict_gp(trend_models[["vague prior"]], x_grid)

plot(
  x_trend,
  y_trend,
  xlim = c(0, 15),
  ylim = range(marginal$latent_interval, y_trend),
  xlab = "x",
  ylab = "f(x)"
)
lines(x_grid, truth(x_grid), col = "grey60", lwd = 2)
lines(x_grid, marginal$mean)
matlines(x_grid, marginal$latent_interval, lty = 1, col = "black")
matlines(x_grid, plug_in$latent_interval, lty = 2, col = "black")
legend(
  "topleft",
  legend = c("vague prior", "REML plug-in", "true function"),
  lty = c(1, 2, 1),
  col = c("black", "black", "grey60"),
  lwd = c(1, 1, 2),
  bty = "n"
)

