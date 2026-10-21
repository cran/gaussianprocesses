# Spectral mixture versus a structured composite kernel on Mauna Loa CO2
#
# Monthly atmospheric CO2 concentrations (datasets::co2) from 1985 to 1994
# train two models, and 1995 to 1997 are held out. Both estimate a linear
# trend by generalized least squares (linear_mean()) and their kernel
# hyperparameters and noise by maximum likelihood:
#
# - a composite kernel built from prior knowledge: a smooth trend
#   deviation, an annual cycle that drifts slowly, and short-term
#   variation;
# - a spectral mixture with three components, initialized from the
#   spectrum of the detrended training series.
#
# Each model is scored by the log predictive density of the 36 held-out
# months, the sum of the log densities of the observations under the
# predictive distribution including noise, and by the RMSE of the
# predictive mean.

library(gaussianprocesses)

years <- as.numeric(stats::time(datasets::co2))
level <- as.numeric(datasets::co2)
training <- years >= 1985 & years < 1995
held_out <- years >= 1995

composite <- sum_kernel(
  rbf_kernel(variance = 1, length_scale = 10),
  product_kernel(
    periodic_kernel(variance = 4, length_scale = 1, period = 1),
    rbf_kernel(variance = 1, length_scale = 50)
  ),
  rbf_kernel(variance = 0.5, length_scale = 1)
)
composite_model <- optimize_gp(
  years[training],
  level[training],
  kernel = composite,
  noise_variance = 0.05,
  mean = linear_mean(),
  n_starts = 3
)

detrended <- stats::residuals(stats::lm(level[training] ~ years[training]))
initial <- initialize_spectral_mixture(
  years[training],
  detrended,
  n_components = 3
)
initial

spectral_model <- optimize_gp(
  years[training],
  level[training],
  kernel = initial,
  noise_variance = 0.05,
  mean = linear_mean(),
  n_starts = 3
)
spectral_model$kernel

score <- function(model) {
  prediction <- predict_gp(model, years[held_out])

  c(
    log_marginal_likelihood = log_marginal_likelihood(model),
    held_out_log_predictive_density = sum(
      stats::dnorm(
        level[held_out],
        mean = prediction$mean,
        sd = prediction$observation_sd,
        log = TRUE
      )
    ),
    held_out_rmse = sqrt(mean((prediction$mean - level[held_out])^2))
  )
}

scores <- rbind(
  composite = score(composite_model),
  spectral_mixture = score(spectral_model)
)
scores

grid <- seq(1985, 1998, by = 1 / 24)
composite_path <- predict_gp(composite_model, grid)
spectral_path <- predict_gp(spectral_model, grid)

plot(
  years[training | held_out],
  level[training | held_out],
  pch = ifelse(held_out[training | held_out], 1, 20),
  cex = 0.6,
  xlab = "year",
  ylab = "CO2 (ppm)"
)
abline(v = 1995, col = "grey60", lty = 3)
lines(grid, composite_path$mean, lty = 2)
lines(grid, spectral_path$mean)
legend(
  "topleft",
  legend = c("spectral mixture", "composite", "held out"),
  lty = c(1, 2, NA),
  pch = c(NA, NA, 1),
  bty = "n"
)
