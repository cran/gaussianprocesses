# Does the jitter policy make results invariant to the units of y?
#
# Run manually with the package attached, for example
#   Rscript -e 'library(gaussianprocesses); source("inst/benchmarks/jitter-scaling.R")'
# This file is intentionally excluded from CI.
#
# Multiplying y by c and every variance parameter (kernel variance and noise
# variance) by c^2 multiplies the covariance C by c^2. Then, in exact
# arithmetic,
#   posterior mean       -> c   * original
#   posterior variance   -> c^2 * original
#   log marginal lik.    -> original - n log|c|
#   gradient in log parameters and seeded draws / c -> unchanged,
# provided the jittered covariance C + delta I also scales by c^2, that is,
# provided delta scales with C. The study compares the two jitter policies
# on problems that need jitter, for c in {1e-4, 1, 1e4}. Each entry is the
# largest relative difference from the c = 1 result after undoing the
# scaling. The policy is switched with the internal option
# gaussianprocesses.jitter_scaling.

cases <- list(
  "duplicated inputs, no noise" = list(
    x = c(-1, 0, 0, 0.4, 1.1, 1.1, 2),
    kernel = function(variance) rbf_kernel(variance, 0.8),
    noise = 0
  ),
  "dense grid, no noise" = list(
    x = seq(0, 1, length.out = 40),
    kernel = function(variance) rbf_kernel(variance, 0.5),
    noise = 0
  ),
  "dense grid, tiny noise" = list(
    x = seq(0, 1, length.out = 60),
    kernel = function(variance) matern52_kernel(variance, 0.6),
    noise = 1e-12
  )
)
scales <- c(1e-4, 1, 1e4)
x_new <- c(-0.5, 0.3, 0.75, 2.5)


relative_difference <- function(actual, expected) {
  max(abs(actual - expected)) / max(abs(expected))
}


scaled_results <- function(case, c) {
  x <- case$x
  y <- c * (sin(3 * x) + 0.2 * x)
  kernel <- case$kernel(1.3 * c^2)

  model <- tryCatch(
    fit_gp(x, y, kernel, noise_variance = case$noise * c^2),
    error = function(error) NULL
  )

  if (is.null(model)) {
    return(NULL)
  }

  prediction <- predict_gp(model, x_new)
  prior <- tryCatch(
    sample_gp_prior(x, kernel, n_draws = 3, seed = 1)$latent,
    error = function(error) NULL
  )
  posterior <- tryCatch(
    sample_gp_posterior(model, x_new, n_draws = 3, seed = 2)$latent,
    error = function(error) NULL
  )

  list(
    jitter = model$jitter,
    attempts = model$cholesky_attempts,
    log_marginal_likelihood = log_marginal_likelihood(model) +
      length(x) * log(c),
    gradient = log_marginal_likelihood_gradient(model),
    mean = prediction$mean / c,
    variance = prediction$latent_variance / c^2,
    prior_draws = if (!is.null(prior)) prior / c,
    posterior_draws = if (!is.null(posterior)) posterior / c
  )
}


rows <- list()
old_option <- getOption("gaussianprocesses.jitter_scaling")

for (policy in c("absolute", "relative")) {
  options(gaussianprocesses.jitter_scaling = policy)

  for (case_name in names(cases)) {
    reference <- scaled_results(cases[[case_name]], 1)

    for (c in scales) {
      result <- scaled_results(cases[[case_name]], c)
      compare <- function(field) {
        if (is.null(result) || is.null(result[[field]])) {
          return(NA_real_)
        }
        relative_difference(result[[field]], reference[[field]])
      }

      rows[[length(rows) + 1L]] <- data.frame(
        policy = policy,
        case = case_name,
        c = c,
        fitted = !is.null(result),
        attempts = if (is.null(result)) NA_integer_ else result$attempts,
        jitter = if (is.null(result)) NA_real_ else result$jitter,
        lml = compare("log_marginal_likelihood"),
        gradient = compare("gradient"),
        mean = compare("mean"),
        variance = compare("variance"),
        prior_draws = compare("prior_draws"),
        posterior_draws = compare("posterior_draws")
      )
    }
  }
}

options(gaussianprocesses.jitter_scaling = old_option)

results <- do.call(rbind, rows)
numeric_columns <- c(
  "jitter", "lml", "gradient", "mean", "variance", "prior_draws",
  "posterior_draws"
)
results[numeric_columns] <- lapply(results[numeric_columns], signif, 2)
print(results, row.names = FALSE)
