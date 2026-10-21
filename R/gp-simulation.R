# Gaussian-process simulation utilities --------------------------------------

.validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(invisible(NULL))
  }

  .assert_numeric_scalar(seed, "seed", lower = 0)

  if (seed != floor(seed) || seed > .Machine$integer.max) {
    .gp_abort(
      sprintf(
        "'seed' must be NULL or an integer between 0 and %d.",
        .Machine$integer.max
      ),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(as.integer(seed))
}


.offset_seed <- function(seed, offset) {
  # Derived seeds are validated before use so that seed + offset cannot
  # overflow R's integer range part-way through a benchmark.
  .validate_seed(as.numeric(seed) + offset)
}


.with_local_seed <- function(seed, code) {
  .validate_seed(seed)

  if (is.null(seed)) {
    return(code())
  }

  had_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }

  on.exit({
    if (had_seed) {
      assign(
        ".Random.seed",
        old_seed,
        envir = .GlobalEnv
      )
    } else if (exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed))
  code()
}


#' Simulate observations from a Gaussian process
#'
#' Draws a latent Gaussian-process realization at supplied inputs and then adds
#' independent Gaussian observation noise. The kernel, mean, noise variances,
#' and numerical jitter used for simulation are retained in the result.
#'
#' @param x Numeric vector or matrix of simulation inputs.
#' @param kernel Gaussian-process kernel specification.
#' @param mean Gaussian-process mean specification. Its coefficients must
#'   be fixed or have a Gaussian [coefficient_prior()], whose uncertainty the
#'   simulated latent function includes.
#' @param noise_variance Non-negative scalar or one non-negative value per
#'   observation.
#' @param seed Optional non-negative integer seed. The caller's global random
#'   state is restored after simulation.
#' @param initial_jitter Non-negative jitter tried first, relative to the
#'   covariance scale (the mean of its diagonal).
#' @param fallback_jitter Positive relative jitter tried after the first
#'   failed factorization.
#' @param jitter_multiplier Multiplicative jitter escalation factor.
#' @param max_attempts Maximum Cholesky attempts.
#' @param symmetry_tolerance Relative tolerance for checking that the
#'   covariance matrix is symmetric.
#'
#' @return An object of class `gaussianprocesses_simulation`.
#' @examples
#' simulation <- simulate_gp_data(
#'   seq(0, 1, length.out = 30),
#'   kernel = rbf_kernel(length_scale = 0.2),
#'   noise_variance = 0.01,
#'   seed = 1
#' )
#' simulation
#'
#' head(cbind(latent = simulation$latent, observed = simulation$observed))
#' @template stable
#' @export
simulate_gp_data <- function(
    x,
    kernel,
    mean = zero_mean(),
    noise_variance = 0,
    seed = NULL,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps)) {
  x_matrix <- .as_input_matrix(x, name = "x")
  .assert_kernel(kernel)
  .assert_kernel_dimensions(kernel, ncol(x_matrix))
  .assert_mean_function(mean)
  .validate_seed(seed)

  noise <- .resolve_training_noise_variance(
    noise_variance,
    n_observations = nrow(x_matrix)
  )

  prior_mean <- evaluate_mean(mean, x_matrix)
  covariance <- evaluate_kernel(kernel, x_matrix)
  coefficient_covariance <- .mean_prior_covariance(mean, x_matrix)

  # Coefficients with a Gaussian prior N(b, B) add h(x)' B h(x').
  if (!is.null(coefficient_covariance)) {
    covariance <- covariance + coefficient_covariance
  }

  decomposition <- .stable_cholesky(
    covariance,
    initial_jitter = initial_jitter,
    fallback_jitter = fallback_jitter,
    jitter_multiplier = jitter_multiplier,
    max_attempts = max_attempts,
    tolerance = symmetry_tolerance
  )

  draw <- .with_local_seed(
    seed,
    function() {
      latent <- prior_mean +
        drop(
          t(decomposition$factor) %*%
            stats::rnorm(nrow(x_matrix))
        )

      observation_error <- if (all(noise == 0)) {
        rep(0, nrow(x_matrix))
      } else {
        stats::rnorm(
          nrow(x_matrix),
          sd = sqrt(noise)
        )
      }

      list(
        latent = latent,
        observation_error = observation_error,
        observed = latent + observation_error
      )
    }
  )

  structure(
    list(
      x = x_matrix,
      latent = draw$latent,
      observed = draw$observed,
      observation_error = draw$observation_error,
      mean = mean,
      mean_values = prior_mean,
      kernel = kernel,
      noise_variance = noise,
      seed = seed,
      simulation_jitter = decomposition$jitter,
      simulation_cholesky_attempts = decomposition$attempts
    ),
    class = "gaussianprocesses_simulation"
  )
}


#' Generate standard GP simulation scenarios
#'
#' Creates deterministic input designs and seeded observations for common
#' validation cases used throughout the benchmark suite.
#'
#' @param scenario One of `"smooth_1d"`, `"near_singular"`,
#'   `"ard_2d"`, or `"heteroscedastic_1d"`.
#' @param n Number of observations.
#' @param seed Non-negative integer seed.
#'
#' @return A `gaussianprocesses_simulation` with a `scenario` element.
#' @examples
#' simulation <- gp_simulation_scenario("ard_2d", n = 40)
#' simulation
#' simulation$kernel
#' @template experimental
#' @export
gp_simulation_scenario <- function(
    scenario = c(
      "smooth_1d",
      "near_singular",
      "ard_2d",
      "heteroscedastic_1d"
    ),
    n = 60L,
    seed = 1729L) {
  scenario <- match.arg(scenario)
  n <- .assert_positive_integer(n, "n")
  .validate_seed(seed)

  if (n < 4L) {
    .gp_abort(
      "'n' must be at least 4 for standard simulation scenarios.",
      "gaussianprocesses_validation_error"
    )
  }

  specification <- switch(
    scenario,
    smooth_1d = {
      x <- seq(-3, 3, length.out = n)

      list(
        x = x,
        kernel = rbf_kernel(
          variance = 1.25,
          length_scale = 0.8
        ),
        noise = 0.05
      )
    },
    near_singular = {
      base <- seq(-2, 2, length.out = n - 1L)
      insertion <- base[[ceiling(length(base) / 2L)]]
      x <- sort(c(base, insertion + 1e-12))

      list(
        x = x,
        kernel = rbf_kernel(
          variance = 1,
          length_scale = 1.5
        ),
        noise = 1e-12
      )
    },
    ard_2d = {
      phase <- seq(0, 2 * pi, length.out = n)
      x <- cbind(
        x1 = seq(-2, 2, length.out = n),
        x2 = cos(phase)
      )

      list(
        x = x,
        kernel = matern32_kernel(
          variance = 1.4,
          length_scale = c(0.45, 2.2)
        ),
        noise = 0.04
      )
    },
    heteroscedastic_1d = {
      x <- seq(0, 1, length.out = n)
      noise <- 0.015 + 0.25 * x^2

      list(
        x = x,
        kernel = matern52_kernel(
          variance = 1,
          length_scale = 0.2
        ),
        noise = noise
      )
    }
  )

  simulation <- simulate_gp_data(
    specification$x,
    kernel = specification$kernel,
    noise_variance = specification$noise,
    seed = seed
  )
  simulation$scenario <- scenario
  simulation
}


#' @export
print.gaussianprocesses_simulation <- function(x, ...) {
  cat("Gaussian-process simulation\n")
  cat("  observations: ", nrow(x$x), "\n", sep = "")
  cat("  input dimensions: ", ncol(x$x), "\n", sep = "")

  if (!is.null(x$scenario)) {
    cat("  scenario: ", x$scenario, "\n", sep = "")
  }

  cat(
    "  noise variance range: [",
    format(signif(min(x$noise_variance), 6), trim = TRUE),
    ", ",
    format(signif(max(x$noise_variance), 6), trim = TRUE),
    "]\n",
    sep = ""
  )
  cat(
    "  simulation jitter: ",
    format(signif(x$simulation_jitter, 6), trim = TRUE),
    "\n",
    sep = ""
  )
  invisible(x)
}
