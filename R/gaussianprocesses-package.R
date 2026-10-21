#' Gaussian Processes from First Principles
#'
#' `gaussianprocesses` provides mathematically explicit implementations of
#' Gaussian-process models with an emphasis on numerical stability,
#' uncertainty quantification, and reproducibility. The mathematical
#' vignettes map each model equation to the corresponding package objects.
#'
#' @section Stability and versioning:
#' Every exported function belongs to one of three tiers, which its help page
#' states in a "Stability" section and by which the reference index on the
#' package website is grouped.
#'
#' - **Stable.** From version 1.0.0, a stable interface changes
#'   incompatibly only in a major release. A stable function, argument, or
#'   result field that is replaced is first deprecated: it keeps working
#'   with a warning for at least one minor release, and is removed in the
#'   next major release.
#' - **Experimental.** An experimental interface may change in a minor
#'   release, and every change is listed in NEWS. The experimental
#'   interfaces are the sparse approximations, latent models with their
#'   likelihoods and scores, the heteroscedastic method, state-space
#'   inference, the spectral-mixture and changepoint kernels, derivative
#'   observations and gradient prediction, hyperparameter uncertainty, the
#'   numerical and calibration diagnostics, and the benchmark and
#'   simulation-scenario helpers.
#' - **Deprecated.** A deprecated name warns when used and will be removed
#'   in the next major release. Version 1.0.0 has no deprecated names.
#'
#' Results and options of a stable function that concern an experimental
#' model class, kernel, or argument follow the experimental tier: for
#' example, [predict_gp()] is stable, but its behaviour for a state-space
#' model is as experimental as state-space inference.
#'
#' Version numbers follow semantic versioning from 1.0.0. A major release
#' may change stable interfaces incompatibly; a minor release adds features
#' and may change experimental interfaces; a patch release only fixes
#' defects compatibly.
#'
#' Fitted models record the format version they were saved with (see
#' [gaussianprocesses_model]). Every 1.x release reads models saved by any
#' earlier release, back to 0.1.0.
#'
#' @keywords internal
"_PACKAGE"
