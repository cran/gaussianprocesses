test_that("stable Cholesky does not add jitter when unnecessary", {
  x <- matrix(c(2, 1, 1, 2), nrow = 2)

  result <- gaussianprocesses:::.stable_cholesky(x)

  expect_s3_class(result, "gaussianprocesses_cholesky")
  expect_equal(result$jitter, 0)
  expect_equal(result$relative_jitter, 0)
  expect_equal(result$covariance_scale, 2)
  expect_equal(result$attempts, 1L)
  expect_equal(crossprod(result$factor), x, tolerance = 1e-12)
})


test_that("stable Cholesky adds explicit jitter for a singular matrix", {
  x <- matrix(1, nrow = 2, ncol = 2)

  result <- gaussianprocesses:::.stable_cholesky(
    x,
    fallback_jitter = 1e-10,
    max_attempts = 3
  )

  expected <- x + diag(result$jitter, nrow(x))

  expect_gt(result$jitter, 0)
  expect_gt(result$attempts, 1L)
  expect_equal(
    crossprod(result$factor),
    expected,
    tolerance = 1e-12
  )
})


test_that("jitter is relative to the mean diagonal of the matrix", {
  # Slightly indefinite (smallest eigenvalue about -5e-14 relative), so it
  # fails to factorize at every scale without jitter.
  indefinite <- matrix(c(1, 1, 0, 1, 1 - 1e-13, 0, 0, 0, 4), nrow = 3)

  for (c2 in c(1e-8, 1, 1e8)) {
    result <- gaussianprocesses:::.stable_cholesky(c2 * indefinite)

    # The same relative jitter is needed at every scale, so the absolute
    # jitter scales with the matrix.
    expect_identical(result$jitter_scaling, "relative")
    expect_equal(result$covariance_scale, c2 * 2)
    expect_identical(result$relative_jitter, 1e-10)
    expect_equal(result$jitter, c2 * 2e-10)
    expect_equal(
      crossprod(result$factor),
      c2 * indefinite + diag(result$jitter, 3),
      tolerance = 1e-12
    )
  }

  # A matrix without a positive diagonal falls back to scale 1.
  expect_identical(gaussianprocesses:::.covariance_scale(matrix(0, 2, 2)), 1)
})


test_that("absolute jitter scaling remains available internally", {
  singular <- 1e8 * matrix(1, nrow = 2, ncol = 2)
  old <- options(gaussianprocesses.jitter_scaling = "absolute")
  on.exit(options(old), add = TRUE)

  result <- gaussianprocesses:::.stable_cholesky(singular, max_attempts = 12)

  expect_identical(result$jitter_scaling, "absolute")
  expect_identical(result$covariance_scale, 1)
  expect_identical(result$jitter, result$relative_jitter)

  options(gaussianprocesses.jitter_scaling = "neither")
  expect_error(
    gaussianprocesses:::.stable_cholesky(diag(2)),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("stable Cholesky reports unrecoverable failure", {
  x <- diag(c(1, -1))

  expect_error(
    gaussianprocesses:::.stable_cholesky(
      x,
      fallback_jitter = 1e-8,
      jitter_multiplier = 10,
      max_attempts = 4
    ),
    "1.000e-06 times the covariance scale 1.000e\\+00",
    class = "gaussianprocesses_numerical_error"
  )
})


test_that("jitter is formatted with its relative size", {
  format_jitter <- gaussianprocesses:::.format_jitter

  expect_identical(format_jitter(0, 2), "0")
  expect_identical(
    format_jitter(2e-10, 2),
    "2e-10 (1e-10 times the covariance scale 2)"
  )
  expect_identical(format_jitter(3e-10, NULL), "3e-10")
})
