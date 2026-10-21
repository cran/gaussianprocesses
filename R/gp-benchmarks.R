# Reproducible Gaussian-process benchmark suite ------------------------------

.elapsed_median <- function(code, repeats) {
  repeats <- .assert_positive_integer(repeats, "repeats")

  timings <- numeric(repeats)

  for (index in seq_len(repeats)) {
    timings[[index]] <- system.time(
      code()
    )[["elapsed"]]
  }

  stats::median(timings)
}


#' Compare two GP posterior predictions numerically
#'
#' @param reference Reference prediction object with `mean` and
#'   `latent_variance` fields.
#' @param candidate Candidate prediction object with matching lengths.
#'
#' @return A named list of numerical accuracy differences.
#' @examples
#' x <- seq(-3, 3, length.out = 100)
#' y <- sin(x)
#' kernel <- rbf_kernel(length_scale = 0.8)
#'
#' exact <- predict_gp(fit_gp(x, y, kernel, noise_variance = 0.01), x)
#' sparse <- predict_sparse_gp(
#'   fit_sparse_gp(x, y, kernel, noise_variance = 0.01, n_inducing = 10),
#'   x
#' )
#'
#' compare_gp_predictions(exact, sparse)
#' @template experimental
#' @export
compare_gp_predictions <- function(
    reference,
    candidate) {
  required <- c("mean", "latent_variance")

  if (!all(required %in% names(reference)) ||
      !all(required %in% names(candidate))) {
    .gp_abort(
      paste0(
        "Both predictions must contain 'mean' and ",
        "'latent_variance'."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (length(reference$mean) != length(candidate$mean) ||
      length(reference$latent_variance) !=
        length(candidate$latent_variance)) {
    .gp_abort(
      "Reference and candidate predictions must have matching lengths.",
      "gaussianprocesses_validation_error"
    )
  }

  mean_difference <- candidate$mean - reference$mean
  variance_difference <-
    candidate$latent_variance -
    reference$latent_variance

  list(
    posterior_mean_rmse =
      sqrt(mean(mean_difference^2)),
    posterior_mean_max_abs_error =
      max(abs(mean_difference)),
    latent_variance_mae =
      mean(abs(variance_difference)),
    latent_variance_max_abs_error =
      max(abs(variance_difference))
  )
}


.simulation_truth_metrics <- function(
    prediction,
    latent_truth) {
  if (length(prediction$mean) != length(latent_truth)) {
    .gp_abort(
      "Prediction and latent truth must have matching lengths.",
      "gaussianprocesses_validation_error"
    )
  }

  residual <- prediction$mean - latent_truth

  list(
    latent_truth_rmse = sqrt(mean(residual^2)),
    latent_truth_mae = mean(abs(residual))
  )
}


.benchmark_single_scenario <- function(
    simulation,
    n_inducing,
    repeats,
    selection) {
  x <- simulation$x
  y <- simulation$observed
  kernel <- simulation$kernel
  noise <- simulation$noise_variance

  exact_fit_time <- .elapsed_median(
    function() {
      fit_gp(
        x,
        y,
        kernel = kernel,
        noise_variance = noise
      )
    },
    repeats = repeats
  )

  exact_model <- fit_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = noise
  )

  exact_predict_time <- .elapsed_median(
    function() {
      predict_gp(
        exact_model,
        x,
        observation_noise_variance = 0
      )
    },
    repeats = repeats
  )

  exact_prediction <- predict_gp(
    exact_model,
    x,
    observation_noise_variance = 0
  )

  sparse_fit_time <- .elapsed_median(
    function() {
      fit_sparse_gp(
        x,
        y,
        kernel = kernel,
        noise_variance = noise,
        n_inducing = n_inducing,
        selection = selection
      )
    },
    repeats = repeats
  )

  sparse_model <- fit_sparse_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = noise,
    n_inducing = n_inducing,
    selection = selection
  )

  sparse_predict_time <- .elapsed_median(
    function() {
      predict_sparse_gp(
        sparse_model,
        x,
        observation_noise_variance = 0
      )
    },
    repeats = repeats
  )

  sparse_prediction <- predict_sparse_gp(
    sparse_model,
    x,
    observation_noise_variance = 0
  )

  approximation_accuracy <- compare_gp_predictions(
    exact_prediction,
    sparse_prediction
  )
  exact_truth <- .simulation_truth_metrics(
    exact_prediction,
    simulation$latent
  )
  sparse_truth <- .simulation_truth_metrics(
    sparse_prediction,
    simulation$latent
  )

  numerical <- gp_numerical_diagnostics(exact_model)

  accuracy <- data.frame(
    model = c("exact", "FITC"),
    latent_truth_rmse = c(
      exact_truth$latent_truth_rmse,
      sparse_truth$latent_truth_rmse
    ),
    latent_truth_mae = c(
      exact_truth$latent_truth_mae,
      sparse_truth$latent_truth_mae
    ),
    posterior_mean_rmse_vs_exact = c(
      0,
      approximation_accuracy$posterior_mean_rmse
    ),
    posterior_mean_max_abs_error_vs_exact = c(
      0,
      approximation_accuracy$posterior_mean_max_abs_error
    ),
    latent_variance_mae_vs_exact = c(
      0,
      approximation_accuracy$latent_variance_mae
    ),
    latent_variance_max_abs_error_vs_exact = c(
      0,
      approximation_accuracy$latent_variance_max_abs_error
    ),
    check.names = FALSE
  )

  performance <- data.frame(
    model = c("exact", "FITC"),
    fit_elapsed_seconds = c(
      exact_fit_time,
      sparse_fit_time
    ),
    predict_elapsed_seconds = c(
      exact_predict_time,
      sparse_predict_time
    ),
    retained_model_bytes = c(
      as.numeric(utils::object.size(exact_model)),
      as.numeric(utils::object.size(sparse_model))
    ),
    n_training = nrow(x),
    n_inducing = c(nrow(x), sparse_model$n_inducing),
    check.names = FALSE
  )

  stability <- data.frame(
    exact_numerical_jitter = numerical$numerical_jitter,
    exact_reciprocal_condition_number =
      numerical$reciprocal_condition_number,
    exact_cholesky_attempts = numerical$cholesky_attempts,
    sparse_kuu_jitter = sparse_model$kuu_jitter,
    sparse_a_jitter = sparse_model$a_jitter,
    sparse_fitc_floor_count = sparse_model$fitc_floor_count,
    check.names = FALSE
  )

  list(
    accuracy = accuracy,
    performance = performance,
    stability = stability
  )
}


#' Run the standard GP benchmark suite
#'
#' Runs seeded synthetic scenarios through exact and FITC regression. Accuracy,
#' performance, and numerical stability are returned as separate tables so
#' faster execution is never conflated with better numerical agreement.
#'
#' @param scenarios Character vector of standard scenario names.
#' @param n Number of observations per scenario.
#' @param n_inducing Number of FITC inducing points.
#' @param repeats Number of timing repetitions. Medians are reported.
#' @param seed Base non-negative integer seed. Scenario-specific offsets are
#'   applied deterministically.
#'
#' @details
#' Timing results use median elapsed wall-clock time and are intended for
#' local comparison, not CI gating. Retained memory is measured with
#' [utils::object.size()] and is not a peak-allocation measurement.
#'
#' The exact Gaussian process is used as the numerical reference when
#' assessing FITC approximation error. Synthetic scenarios also retain their
#' latent draw so both exact and sparse posterior means can be compared
#' against the same known simulation truth.
#'
#' @return An object of class `gaussianprocesses_benchmark_suite` containing
#'   separate `accuracy`, `performance`, and `stability` data frames.
#' @examples
#' # A deliberately small run; see inst/benchmarks/run-suite.R for a full one.
#' suite <- run_gp_benchmark_suite(
#'   scenarios = "smooth_1d",
#'   n = 30,
#'   n_inducing = 8,
#'   repeats = 1
#' )
#' suite$accuracy
#' @template experimental
#' @export
run_gp_benchmark_suite <- function(
    scenarios = c(
      "smooth_1d",
      "near_singular",
      "ard_2d",
      "heteroscedastic_1d"
    ),
    n = 60L,
    n_inducing = 15L,
    repeats = 3L,
    seed = 1729L) {
  allowed <- c(
    "smooth_1d",
    "near_singular",
    "ard_2d",
    "heteroscedastic_1d"
  )

  if (!is.character(scenarios) ||
      length(scenarios) == 0L ||
      any(!scenarios %in% allowed)) {
    .gp_abort(
      paste0(
        "'scenarios' must contain only: ",
        paste(allowed, collapse = ", "),
        "."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  n <- .assert_positive_integer(n, "n")
  n_inducing <- .assert_positive_integer(
    n_inducing,
    "n_inducing"
  )
  repeats <- .assert_positive_integer(repeats, "repeats")
  .validate_seed(seed)
  .offset_seed(seed, length(scenarios) - 1L)

  if (n_inducing > n) {
    .gp_abort(
      "'n_inducing' cannot exceed 'n'.",
      "gaussianprocesses_validation_error"
    )
  }

  accuracy_blocks <- vector("list", length(scenarios))
  performance_blocks <- vector("list", length(scenarios))
  stability_blocks <- vector("list", length(scenarios))

  for (index in seq_along(scenarios)) {
    scenario <- scenarios[[index]]
    simulation <- gp_simulation_scenario(
      scenario,
      n = n,
      seed = .offset_seed(seed, index - 1L)
    )

    selection <- if (ncol(simulation$x) == 1L) {
      "quantile"
    } else {
      "farthest"
    }

    benchmark <- .benchmark_single_scenario(
      simulation,
      n_inducing = n_inducing,
      repeats = repeats,
      selection = selection
    )

    benchmark$accuracy$scenario <- scenario
    benchmark$performance$scenario <- scenario
    benchmark$stability$scenario <- scenario

    accuracy_blocks[[index]] <- benchmark$accuracy
    performance_blocks[[index]] <- benchmark$performance
    stability_blocks[[index]] <- benchmark$stability
  }

  accuracy <- do.call(rbind, accuracy_blocks)
  performance <- do.call(rbind, performance_blocks)
  stability <- do.call(rbind, stability_blocks)

  rownames(accuracy) <- NULL
  rownames(performance) <- NULL
  rownames(stability) <- NULL

  accuracy <- accuracy[
    c(
      "scenario",
      setdiff(names(accuracy), "scenario")
    )
  ]
  performance <- performance[
    c(
      "scenario",
      setdiff(names(performance), "scenario")
    )
  ]
  stability <- stability[
    c(
      "scenario",
      setdiff(names(stability), "scenario")
    )
  ]

  structure(
    list(
      accuracy = accuracy,
      performance = performance,
      stability = stability,
      settings = list(
        n = n,
        n_inducing = n_inducing,
        repeats = repeats,
        seed = seed,
        scenarios = scenarios
      )
    ),
    class = "gaussianprocesses_benchmark_suite"
  )
}


#' Benchmark exact-GP scaling
#'
#' Measures exact GP fit time and retained model size over increasing training
#' sizes using seeded RBF simulations. This function is intended for local
#' empirical scaling studies rather than CI performance gating.
#'
#' @param sizes Positive integer vector of training sizes.
#' @param repeats Number of timing repetitions per size.
#' @param seed Base non-negative integer seed.
#'
#' @return A data frame with one row per requested size.
#' @examples
#' # A deliberately small run; larger sizes show the cubic cost.
#' benchmark_exact_gp_scaling(sizes = c(20, 40), repeats = 1)
#' @template experimental
#' @export
benchmark_exact_gp_scaling <- function(
    sizes = c(50L, 100L, 200L, 400L),
    repeats = 3L,
    seed = 1729L) {
  if (!is.numeric(sizes) ||
      !is.null(dim(sizes)) ||
      length(sizes) == 0L ||
      any(!is.finite(sizes)) ||
      any(sizes < 4) ||
      any(sizes != floor(sizes))) {
    .gp_abort(
      "'sizes' must contain integers greater than or equal to 4.",
      "gaussianprocesses_validation_error"
    )
  }

  repeats <- .assert_positive_integer(repeats, "repeats")
  .validate_seed(seed)
  .offset_seed(seed, length(sizes) - 1L)

  blocks <- vector("list", length(sizes))

  for (index in seq_along(sizes)) {
    n <- as.integer(sizes[[index]])
    simulation <- gp_simulation_scenario(
      "smooth_1d",
      n = n,
      seed = .offset_seed(seed, index - 1L)
    )

    elapsed <- .elapsed_median(
      function() {
        fit_gp(
          simulation$x,
          simulation$observed,
          kernel = simulation$kernel,
          noise_variance = simulation$noise_variance
        )
      },
      repeats = repeats
    )

    model <- fit_gp(
      simulation$x,
      simulation$observed,
      kernel = simulation$kernel,
      noise_variance = simulation$noise_variance
    )
    numerical <- gp_numerical_diagnostics(model)

    blocks[[index]] <- data.frame(
      n = n,
      median_fit_elapsed_seconds = elapsed,
      retained_model_bytes =
        as.numeric(utils::object.size(model)),
      numerical_jitter = numerical$numerical_jitter,
      reciprocal_condition_number =
        numerical$reciprocal_condition_number,
      check.names = FALSE
    )
  }

  do.call(rbind, blocks)
}


#' @export
print.gaussianprocesses_benchmark_suite <- function(x, ...) {
  cat("Gaussian-process benchmark suite\n")
  cat(
    "  scenarios: ",
    length(x$settings$scenarios),
    "\n",
    sep = ""
  )
  cat(
    "  observations per scenario: ",
    x$settings$n,
    "\n",
    sep = ""
  )
  cat(
    "  FITC inducing points: ",
    x$settings$n_inducing,
    "\n",
    sep = ""
  )
  cat(
    "  timing repetitions: ",
    x$settings$repeats,
    "\n",
    sep = ""
  )
  invisible(x)
}
