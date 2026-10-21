# Coverage and calibration of Poisson Gaussian-process models on counts from
# a known log-intensity process.
#
# Each replicate draws a log intensity f from a Gaussian process with mean
# log(5) and an RBF kernel (variance 0.5, length scale 1.5) at 160 uniform
# inputs on [0, 10], exposures uniform on [0.5, 2], and counts
# y ~ Poisson(E exp(f)). The first 60 are for training, the other 100 for
# testing. The model is fitted with the true hyperparameters and mean
# ("known") and with hyperparameters and mean level estimated by
# optimize_latent_gp() ("estimated"). For each fit the script records, at
# the test inputs,
#   the coverage of the rate intervals E exp(mu -+ z sigma) of the true
#     rate E exp(f) at levels 0.8 and 0.95;
#   the coverage and mean width of the count prediction intervals;
#   the randomized PIT values of the test counts and their dispersion.
# A last set of replicates adds overdispersion, multiplying each rate by an
# independent gamma variable with mean 1 and variance 0.5, which the Poisson
# model does not describe.

library(gaussianprocesses)

n_replicates <- 150
n_overdispersed <- 40
n_training <- 60
n_test <- 100
levels <- c(0.8, 0.95)
kernel <- rbf_kernel(variance = 0.5, length_scale = 1.5)
log_level <- log(5)

simulate_counts <- function(overdispersion = 0) {
  x <- stats::runif(n_training + n_test, 0, 10)
  f <- log_level + drop(sample_gp_prior(x, kernel)$latent)
  exposure <- stats::runif(length(x), 0.5, 2)
  rate <- exposure * exp(f)

  if (overdispersion > 0) {
    rate <- rate * stats::rgamma(length(x), 1 / overdispersion,
                                 1 / overdispersion)
  }

  list(x = x, f = f, exposure = exposure, y = stats::rpois(length(x), rate))
}

# The rate intervals at both levels come from one prediction: they are the
# exponential of the latent intervals, scaled by the exposure.
evaluate <- function(model, data, test) {
  prediction <- predict_latent_gp(model, data$x[test],
                                  exposure = data$exposure[test])
  true_rate <- data$exposure[test] * exp(data$f[test])
  rate_coverage <- vapply(levels, function(level) {
    half_width <- stats::qnorm((1 + level) / 2) * prediction$latent_sd
    lower <- data$exposure[test] * exp(prediction$mean - half_width)
    upper <- data$exposure[test] * exp(prediction$mean + half_width)
    mean(true_rate >= lower & true_rate <= upper)
  }, numeric(1))
  scores <- gp_holdout_scores(model, data$x[test], data$y[test],
                              interval_levels = levels,
                              exposure = data$exposure[test])

  list(
    intervals = data.frame(
      level = levels,
      rate_coverage = rate_coverage,
      count_coverage = scores$calibration$coverage,
      count_width = scores$calibration$width
    ),
    pit = scores$pointwise$pit,
    dispersion = scores$summary[["dispersion"]]
  )
}

set.seed(20261008)
results <- list()
pit <- list()
dispersion <- list()
estimates <- matrix(NA_real_, n_replicates, 3,
                    dimnames = list(NULL, c("variance", "length_scale",
                                            "mean_level")))
training <- seq_len(n_training)
test <- n_training + seq_len(n_test)

for (replicate in seq_len(n_replicates + n_overdispersed)) {
  overdispersed <- replicate > n_replicates
  data <- simulate_counts(if (overdispersed) 0.5 else 0)
  fits <- list(known = fit_latent_gp(
    data$x[training], data$y[training], kernel, poisson_likelihood(),
    mean = constant_mean(log_level), exposure = data$exposure[training]
  ))

  if (!overdispersed) {
    fits$estimated <- optimize_latent_gp(
      data$x[training], data$y[training], rbf_kernel(), poisson_likelihood(),
      mean = constant_mean(estimate_coefficients()),
      exposure = data$exposure[training], n_starts = 1
    )
    estimates[replicate, ] <- c(
      kernel_parameters(fits$estimated$kernel, flatten = TRUE),
      coef(fits$estimated)
    )
  }

  for (fit in names(fits)) {
    label <- if (overdispersed) "overdispersed" else fit
    evaluation <- evaluate(fits[[fit]], data, test)
    results[[length(results) + 1L]] <- cbind(
      replicate = replicate,
      fit = label,
      evaluation$intervals
    )
    pit[[label]] <- rbind(pit[[label]], evaluation$pit)
    dispersion[[label]] <- c(dispersion[[label]], evaluation$dispersion)
  }
}

results <- do.call(rbind, results)

# Coverage averaged over replicates, with Monte Carlo standard errors from
# the spread of the per-replicate coverages.
coverage <- do.call(rbind, lapply(
  split(results, list(results$fit, results$level), drop = TRUE),
  function(rows) {
    data.frame(
      fit = rows$fit[1],
      level = rows$level[1],
      rate_coverage = mean(rows$rate_coverage),
      rate_coverage_se = stats::sd(rows$rate_coverage) / sqrt(nrow(rows)),
      count_coverage = mean(rows$count_coverage),
      count_coverage_se = stats::sd(rows$count_coverage) / sqrt(nrow(rows)),
      count_width = mean(rows$count_width)
    )
  }
))
coverage <- coverage[order(coverage$fit, coverage$level), ]
rownames(coverage) <- NULL

# PIT values of one data set share its latent function, so they are not
# independent. The Kolmogorov-Smirnov test uses the first test count of each
# data set, which are; the mean and variance pool all test counts.
calibration <- data.frame(
  fit = names(pit),
  pit_mean = vapply(pit, mean, numeric(1)),
  pit_variance = vapply(pit, function(values) stats::var(c(values)),
                        numeric(1)),
  ks_p_value = vapply(
    pit,
    function(values) {
      suppressWarnings(stats::ks.test(values[, 1], "punif"))$p.value
    },
    numeric(1)
  ),
  dispersion = vapply(dispersion, mean, numeric(1))
)
rownames(calibration) <- NULL

print(coverage, digits = 3)
print(calibration, digits = 3)
print(round(apply(estimates, 2, stats::median), 2))
