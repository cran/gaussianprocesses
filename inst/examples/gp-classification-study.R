# Calibration of Gaussian-process classification on data from a known GP
# classifier.
#
# Each replicate draws a latent function f from a Gaussian process with an
# RBF kernel (variance 4, length scale 1) at 280 uniform inputs on [0, 10],
# and binary responses y ~ Bernoulli(plogis(f)). The first 80 are for
# training, the other 200 for testing. Four forecasts of the test responses
# are scored:
#   oracle     plogis(f), the true probabilities;
#   known      the predictive probability with the true hyperparameters;
#   plug-in    plogis() of the latent mean with the true hyperparameters,
#              which ignores the latent uncertainty;
#   estimated  the predictive probability with hyperparameters estimated
#              by optimize_latent_gp().
# Scores are averaged over replicates, with Monte Carlo standard errors, and
# the test forecasts of all replicates are pooled for reliability curves.

library(gaussianprocesses)

n_replicates <- 300
n_training <- 80
n_test <- 200
kernel <- rbf_kernel(variance = 4, length_scale = 1)
methods <- c("oracle", "known", "plug-in", "estimated")

set.seed(20261007)
scores <- list()
forecasts <- list()
estimates <- matrix(NA_real_, n_replicates, 2,
                    dimnames = list(NULL, c("variance", "length_scale")))
bound_warnings <- 0L

for (replicate in seq_len(n_replicates)) {
  x <- stats::runif(n_training + n_test, 0, 10)
  f <- drop(sample_gp_prior(x, kernel)$latent)
  y <- stats::rbinom(length(x), 1, stats::plogis(f))
  training <- seq_len(n_training)
  test <- n_training + seq_len(n_test)

  known <- fit_latent_gp(x[training], y[training], kernel,
                         bernoulli_likelihood())
  known_prediction <- predict_latent_gp(known, x[test])

  estimated <- withCallingHandlers(
    optimize_latent_gp(x[training], y[training], rbf_kernel(),
                       bernoulli_likelihood(), n_starts = 1),
    gaussianprocesses_bound_warning = function(condition) {
      bound_warnings <<- bound_warnings + 1L
      invokeRestart("muffleWarning")
    }
  )
  estimates[replicate, ] <- kernel_parameters(estimated$kernel,
                                              flatten = TRUE)

  probability <- list(
    oracle = stats::plogis(f[test]),
    known = known_prediction$probability,
    `plug-in` = stats::plogis(known_prediction$mean),
    estimated = predict_latent_gp(estimated, x[test])$probability
  )

  for (method in methods) {
    result <- gp_classification_scores(y[test], probability[[method]])
    scores[[length(scores) + 1L]] <- data.frame(
      replicate = replicate,
      method = method,
      brier_score = result$summary[["brier_score"]],
      log_loss = result$summary[["log_loss"]],
      ece = result$summary[["ece"]]
    )
    forecasts[[length(forecasts) + 1L]] <- data.frame(
      method = method,
      probability = probability[[method]],
      observed = y[test]
    )
  }
}

scores <- do.call(rbind, scores)
forecasts <- do.call(rbind, forecasts)

# Means over replicates with Monte Carlo standard errors.
summarize <- function(values) {
  c(mean = mean(values), se = stats::sd(values) / sqrt(length(values)))
}
study <- do.call(rbind, lapply(methods, function(method) {
  rows <- scores[scores$method == method, ]
  brier <- summarize(rows$brier_score)
  log_loss <- summarize(rows$log_loss)
  ece <- summarize(rows$ece)
  data.frame(
    method = method,
    brier_score = brier[["mean"]],
    brier_se = brier[["se"]],
    log_loss = log_loss[["mean"]],
    log_loss_se = log_loss[["se"]],
    ece = ece[["mean"]],
    ece_se = ece[["se"]]
  )
}))

# Paired differences within replicates, whose standard errors are much
# smaller than those of the scores themselves.
paired <- function(method, baseline) {
  difference <- function(score) {
    values <- scores[[score]][scores$method == method] -
      scores[[score]][scores$method == baseline]
    c(mean(values), stats::sd(values) / sqrt(n_replicates))
  }
  brier <- difference("brier_score")
  log_loss <- difference("log_loss")

  data.frame(
    comparison = paste(method, "-", baseline),
    brier_score = brier[[1L]],
    brier_se = brier[[2L]],
    log_loss = log_loss[[1L]],
    log_loss_se = log_loss[[2L]]
  )
}
differences <- rbind(
  paired("known", "oracle"),
  paired("plug-in", "known"),
  paired("estimated", "known")
)

# Pooled reliability curves: test forecasts of every replicate in ten
# equal-width bins.
reliability <- do.call(rbind, lapply(methods, function(method) {
  rows <- forecasts[forecasts$method == method, ]
  bins <- gp_classification_scores(rows$observed, rows$probability)$reliability
  cbind(method = method, bins)
}))
reliability$z <- (reliability$observed_frequency -
                    reliability$mean_probability) / reliability$standard_error

# Calibration of the pooled forecasts: the ECE of all test forecasts
# together, and the largest bin deviation in standard errors.
pooled <- do.call(rbind, lapply(methods, function(method) {
  rows <- forecasts[forecasts$method == method, ]
  bins <- reliability[reliability$method == method, ]
  data.frame(
    method = method,
    ece = gp_classification_scores(rows$observed, rows$probability)$
      summary[["ece"]],
    largest_deviation = max(abs(bins$z))
  )
}))

print(study, digits = 3)
print(differences, digits = 3)
print(pooled, digits = 3)
cat("median estimates:", round(apply(estimates, 2, stats::median), 2),
    "\nbound warnings:", bound_warnings, "\n")
