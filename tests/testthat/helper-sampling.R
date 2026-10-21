# Moment checks for Gaussian draws, shared by the sampling tests.
#
# Statistical checks use fixed seeds, so they are deterministic, and bounds of
# five Monte Carlo standard errors, so they stay valid if R's random number
# generator ever changes.
# CRAN asks for checks of about a minute, and R CMD check runs the tests
# without NOT_CRAN, so they use 2000 draws there; testthat::test_local(),
# devtools::test(), and CI use 20000. The five-standard-error bounds hold for
# any number of draws.
n_statistical_draws <- if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  20000L
} else {
  2000L
}

expect_gaussian_moments <- function(draws, mean, covariance, label) {
  n_draws <- ncol(draws)
  sample_covariance <- stats::cov(t(draws))

  # Sample means have standard error sqrt(variance / N).
  expect_lt(
    max(abs(rowMeans(draws) - mean) / sqrt(diag(covariance) / n_draws)),
    5,
    label = paste(label, "mean")
  )
  # Sample variances have relative standard error sqrt(2 / N).
  expect_lt(
    max(abs(diag(sample_covariance) / diag(covariance) - 1)),
    5 * sqrt(2 / n_draws),
    label = paste(label, "variance")
  )
  # Sample correlations have standard error at most 1 / sqrt(N).
  expect_lt(
    max(abs(stats::cov2cor(sample_covariance) - stats::cov2cor(covariance))),
    5 / sqrt(n_draws),
    label = paste(label, "correlation")
  )
}
