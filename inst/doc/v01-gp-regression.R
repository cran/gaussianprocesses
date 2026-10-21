## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)
library(gaussianprocesses)

## -----------------------------------------------------------------------------
x <- c(-1.5, -0.5, 0.5, 1.5)
y <- c(0.9, -0.1, 0.2, 1.0)

kernel <- rbf_kernel(
  variance = 1.2,
  length_scale = 0.7
)

model <- fit_gp(
  x,
  y,
  kernel = kernel,
  noise_variance = 0.05
)

model

## -----------------------------------------------------------------------------
K <- evaluate_kernel(model$kernel, model$x)
C <- K + diag(model$training_noise_variance)
C

model$cholesky
all.equal(crossprod(model$cholesky), C + diag(model$jitter, nrow(C)))

model$alpha

## -----------------------------------------------------------------------------
x_new <- seq(-2, 2, length.out = 80)

prediction <- predict_gp(
  model,
  x_new,
  interval_level = 0.95
)

head(
  data.frame(
    x = x_new,
    mean = prediction$mean,
    latent_sd = prediction$latent_sd,
    observation_sd = prediction$observation_sd
  )
)

## ----fig.alt = "Posterior mean of the fitted Gaussian process for x from -2 to 2, with a hatched 95% prediction band and the four training observations as points."----
plot(
  x_new,
  prediction$mean,
  type = "l",
  xlab = "x",
  ylab = "response",
  ylim = range(prediction$prediction_interval)
)
polygon(
  c(x_new, rev(x_new)),
  c(
    prediction$prediction_interval[, "lower"],
    rev(prediction$prediction_interval[, "upper"])
  ),
  border = NA,
  density = 15
)
lines(x_new, prediction$mean)
points(x, y)

## ----fig.alt = "Two panels. Left: five smooth prior function draws spread around zero across x from -2 to 2. Right: five posterior function draws that pass close to the four training observations and fan out between and beyond them.", fig.width = 7, fig.height = 3.5----
prior_draws <- sample_gp_prior(x_new, kernel, n_draws = 5, seed = 1)
posterior_draws <- sample_gp_posterior(model, x_new, n_draws = 5, seed = 1)

old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
matplot(
  x_new, prior_draws$latent,
  type = "l", lty = 1, col = "grey40",
  xlab = "x", ylab = "f(x)", main = "Prior draws"
)
matplot(
  x_new, posterior_draws$latent,
  type = "l", lty = 1, col = "grey40",
  xlab = "x", ylab = "f(x)", main = "Posterior draws"
)
points(x, y, pch = 19)
par(old_par)

## -----------------------------------------------------------------------------
c(prior = prior_draws$jitter, posterior = posterior_draws$jitter)

## ----fig.alt = "One latent posterior function draw as a line, with the corresponding noisy observation draws as small points scattered around it, and the training observations as larger points."----
noisy <- sample_gp_posterior(
  model,
  x_new,
  n_draws = 1,
  type = "observation",
  seed = 2
)

plot(
  x_new, noisy$observation[, 1],
  pch = 20, cex = 0.6, col = "grey50",
  xlab = "x", ylab = "response"
)
lines(x_new, noisy$latent[, 1])
points(x, y, pch = 19)

## -----------------------------------------------------------------------------
many <- sample_gp_posterior(
  model,
  c(-1, 0.5),
  n_draws = 5000,
  type = "observation",
  seed = 3
)

round(
  rbind(
    latent_sample = apply(many$latent, 1, var),
    latent_theory = diag(many$latent_covariance),
    observation_sample = apply(many$observation, 1, var),
    observation_theory = diag(many$latent_covariance) +
      many$observation_noise_variance
  ),
  4
)

## -----------------------------------------------------------------------------
known_noise <- c(0.02, 0.04, 0.08, 0.16)

heterogeneous_fit <- fit_gp(
  x,
  y,
  kernel = kernel,
  noise_variance = known_noise
)

heterogeneous_fit$training_noise_variance

