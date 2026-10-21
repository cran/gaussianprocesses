# API stability tiers and deprecations -----------------------------------------
#
# Every export is stable, experimental, or deprecated (#39). The policy is in
# ?gaussianprocesses-package and CONTRIBUTING.md. Each export's help page
# states its tier in a "Stability" section, made from the templates in
# man-roxygen/, and the pkgdown reference index groups the topics by tier.
# tests/testthat/test-api-stability.R checks that this file, the help pages,
# and the reference index agree. Exports not listed below are stable.

.experimental_exports <- c(
  # Approximate inference: sparse approximations, latent models and their
  # likelihoods and scores, and the heteroscedastic plug-in method.
  "fit_sparse_gp", "optimize_sparse_gp", "predict_sparse_gp",
  "select_inducing_points",
  "fit_latent_gp", "optimize_latent_gp", "predict_latent_gp",
  "gaussian_likelihood", "bernoulli_likelihood", "poisson_likelihood",
  "evaluate_likelihood", "likelihood_predictive",
  "gp_classification_scores", "gp_count_scores",
  "fit_heteroscedastic_gp", "predict_heteroscedastic_gp",
  # Benchmark and simulation-scenario helpers.
  "benchmark_exact_gp_scaling", "benchmark_sparse_gp",
  "run_gp_benchmark_suite", "compare_gp_predictions",
  "gp_simulation_scenario",
  # Recent kernels and methods. fit_gp(derivative = ),
  # fit_time_series_gp(method = "state_space"), and
  # time_series_kernel(changepoint = ) are experimental parts of stable
  # functions, stated on their help pages.
  "optimize_time_series_gp", "spectral_mixture_kernel",
  "initialize_spectral_mixture", "changepoint_kernel",
  "predict_gradient_gp", "kernel_input_gradient",
  # Hyperparameter uncertainty and diagnostics.
  "gp_hyperparameter_uncertainty", "gp_profile_likelihood",
  "gp_numerical_diagnostics", "gp_calibration_diagnostics"
)

# Deprecated names: none in 1.0.0. A stable name deprecated later is listed
# here and on ?gaussianprocesses-deprecated, warns for at least one minor
# release, and is removed in the next major release (CONTRIBUTING.md).
.deprecated_exports <- stats::setNames(character(0), character(0))

.deprecated_arguments <- data.frame(
  fun = character(0),
  argument = character(0),
  replacement = character(0),
  stringsAsFactors = FALSE
)

# Result fields kept as deprecated copies, which cannot warn when read.
.deprecated_fields <- data.frame(
  fun = character(0),
  field = character(0),
  replacement = character(0),
  stringsAsFactors = FALSE
)

# Names deprecated after 0.1.0 and removed in 1.0.0 without a release of
# warnings, by the decision in #40. The tests check that none returns.
.removed_in_1.0.0 <- list(
  exports = c(fit_gp_time = "fit_time_series_gp"),
  arguments = c(
    tolerance = "symmetry_tolerance or, in fit_heteroscedastic_gp(), convergence_tolerance",
    numerical_tolerance = "symmetry_tolerance"
  ),
  fields = c(noise_variance = "observation_noise_variance")
)


.api_tier <- function(name) {
  ifelse(
    name %in% names(.deprecated_exports),
    "deprecated",
    ifelse(name %in% .experimental_exports, "experimental", "stable")
  )
}
