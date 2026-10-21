# Numerical linear-algebra helpers -------------------------------------------

# Scale of a covariance matrix: the mean of its diagonal, the average prior
# variance. Jitter is a multiple of this scale, so a jittered factorization
# changes the matrix by the same relative amount whatever the units of the
# data. A matrix with no positive finite diagonal has scale 1.
.diagonal_scale <- function(diagonal) {
  scale <- mean(abs(diagonal))

  if (!is.finite(scale) || scale <= 0) {
    return(1)
  }

  scale
}


.covariance_scale <- function(x) {
  .diagonal_scale(diag(x))
}


# Jitter as printed by models, diagnostics, and samples. Objects created
# before the covariance scale was recorded print the absolute jitter only.
.format_jitter <- function(jitter, covariance_scale) {
  if (jitter == 0 || is.null(covariance_scale)) {
    return(format(signif(jitter, 3), trim = TRUE))
  }

  sprintf(
    "%s (%s times the covariance scale %s)",
    format(signif(jitter, 3), trim = TRUE),
    format(signif(jitter / covariance_scale, 3), trim = TRUE),
    format(signif(covariance_scale, 3), trim = TRUE)
  )
}


# How jitter is scaled: "relative" multiplies the jitter ladder by
# .covariance_scale(); "absolute" uses the ladder as it is, as earlier
# versions did. The option exists to reproduce earlier results and for the
# study in inst/benchmarks/jitter-scaling.R; it is not part of the documented
# API.
.jitter_scaling <- function() {
  scaling <- getOption("gaussianprocesses.jitter_scaling", "relative")

  if (!identical(scaling, "relative") && !identical(scaling, "absolute")) {
    .gp_abort(
      paste0(
        "Option 'gaussianprocesses.jitter_scaling' must be \"relative\" ",
        "or \"absolute\"."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  scaling
}


# Cholesky factor R of x + delta I with R'R = x + delta I. The jitter delta is
# 0 if x factorizes; otherwise it climbs the ladder initial_jitter (if
# positive), fallback_jitter, fallback_jitter * jitter_multiplier, ..., with
# every value multiplied by the covariance scale under relative scaling.
.stable_cholesky <- function(
    x,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    tolerance = .default_tolerance()) {
  .assert_symmetric_matrix(x, tolerance = tolerance)
  .assert_numeric_scalar(initial_jitter, "initial_jitter", lower = 0)
  .assert_numeric_scalar(
    fallback_jitter,
    "fallback_jitter",
    lower = 0,
    lower_open = TRUE
  )
  .assert_numeric_scalar(
    jitter_multiplier,
    "jitter_multiplier",
    lower = 1,
    lower_open = TRUE
  )
  max_attempts <- .assert_positive_integer(max_attempts, "max_attempts")

  scaling <- .jitter_scaling()
  scale <- if (identical(scaling, "relative")) .covariance_scale(x) else 1
  relative_jitter <- initial_jitter

  for (attempt in seq_len(max_attempts)) {
    candidate <- x

    if (relative_jitter > 0) {
      diag(candidate) <- diag(candidate) + relative_jitter * scale
    }

    factor <- tryCatch(
      chol(candidate),
      error = function(error) NULL
    )

    if (!is.null(factor)) {
      return(
        structure(
          list(
            factor = factor,
            jitter = relative_jitter * scale,
            relative_jitter = relative_jitter,
            covariance_scale = scale,
            jitter_scaling = scaling,
            attempts = attempt
          ),
          class = "gaussianprocesses_cholesky"
        )
      )
    }

    relative_jitter <- if (relative_jitter == 0) {
      fallback_jitter
    } else {
      relative_jitter * jitter_multiplier
    }
  }

  last_relative_jitter <- if (max_attempts == 1L) {
    initial_jitter
  } else {
    relative_jitter / jitter_multiplier
  }

  .gp_abort(
    sprintf(
      paste0(
        "Cholesky factorization failed after %d attempts. ",
        "The last attempted diagonal jitter was %.3e (%.3e times the ",
        "covariance scale %.3e)."
      ),
      max_attempts,
      last_relative_jitter * scale,
      last_relative_jitter,
      scale
    ),
    "gaussianprocesses_numerical_error"
  )
}
