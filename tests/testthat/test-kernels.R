test_that("RBF kernel has the expected one-dimensional form", {
  x <- c(0, 1, 2)
  covariance <- kernel_rbf(x, variance = 2, length_scale = 1)

  expected <- 2 * exp(
    -0.5 * outer(x, x, FUN = function(a, b) (a - b)^2)
  )

  expect_equal(covariance, expected, tolerance = 1e-12)
  expect_equal(diag(covariance), rep(2, length(x)))
  expect_equal(covariance, t(covariance), tolerance = 1e-12)
})


test_that("kernel cross-covariances have consistent dimensions", {
  x <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  y <- matrix(c(0.5, 0.5, 3, 3), ncol = 2, byrow = TRUE)

  kernels <- list(
    function(a, b) kernel_rbf(a, b),
    function(a, b) kernel_matern12(a, b),
    function(a, b) kernel_matern32(a, b),
    function(a, b) kernel_matern52(a, b),
    function(a, b) kernel_rational_quadratic(a, b),
    function(a, b) kernel_periodic(a, b),
    function(a, b) kernel_linear(a, b),
    function(a, b) kernel_white_noise(a, b)
  )

  for (kernel in kernels) {
    expect_equal(dim(kernel(x, y)), c(3L, 2L))
  }
})


test_that("stationary kernels put marginal variance on the diagonal", {
  x <- matrix(
    c(
      0, 0,
      1, 2,
      3, 4
    ),
    ncol = 2,
    byrow = TRUE
  )
  variance <- 1.7

  kernels <- list(
    kernel_rbf,
    kernel_matern12,
    kernel_matern32,
    kernel_matern52,
    kernel_rational_quadratic,
    kernel_periodic
  )

  for (kernel in kernels) {
    covariance <- kernel(x, variance = variance)
    expect_equal(diag(covariance), rep(variance, nrow(x)), tolerance = 1e-12)
    expect_equal(covariance, t(covariance), tolerance = 1e-12)
  }
})


test_that("Matérn one-half is exponential in distance", {
  x <- c(0, 2)
  covariance <- kernel_matern12(x, variance = 3, length_scale = 2)

  expect_equal(covariance[1, 2], 3 * exp(-1), tolerance = 1e-12)
})


test_that("rational quadratic approaches RBF for large alpha", {
  x <- c(0, 0.5, 1.5)

  rq <- kernel_rational_quadratic(
    x,
    variance = 1.3,
    length_scale = 0.8,
    alpha = 1e7
  )
  rbf <- kernel_rbf(x, variance = 1.3, length_scale = 0.8)

  expect_equal(rq, rbf, tolerance = 1e-7)
})


test_that("periodic kernel repeats at the requested period in one dimension", {
  x <- c(0, 2, 4)
  covariance <- kernel_periodic(
    x,
    variance = 2,
    length_scale = 0.7,
    period = 2
  )

  expect_equal(covariance, matrix(2, 3, 3), tolerance = 1e-12)
})


test_that("linear kernel matches scaled inner products", {
  x <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
  y <- matrix(c(5, 6), ncol = 2)

  covariance <- kernel_linear(x, y, variance = 2)

  expect_equal(covariance, 2 * tcrossprod(x, y))
})


test_that("white-noise kernel identifies coincident observations", {
  x <- c(0, 1, 1, 2)
  covariance <- kernel_white_noise(x, variance = 0.25)

  expect_equal(diag(covariance), rep(0.25, length(x)))
  expect_equal(covariance[2, 3], 0.25)
  expect_equal(covariance[1, 4], 0)
})


test_that("white-noise kernel is exact on coincident multidimensional rows", {
  # The expanded form ||x||^2 + ||y||^2 - 2 x'y left round-off on several of
  # these diagonal entries, so white noise returned zero instead of variance.
  x <- cbind(3 * sin(1:20), 3 * cos(1:20), (1:20) / 7)

  covariance <- kernel_white_noise(x, variance = 0.4)

  expect_equal(diag(covariance), rep(0.4, nrow(x)))
  expect_equal(covariance, diag(0.4, nrow(x)))
  expect_equal(
    diag(gaussianprocesses:::.pairwise_squared_distance(x, x)),
    rep(0, nrow(x))
  )
  expect_equal(
    kernel_white_noise(x, x[c(3, 5), ], variance = 0.4)[c(3, 5), ],
    diag(0.4, 2)
  )
})


test_that("stationary kernels are translation invariant far from the origin", {
  x <- cbind(c(0, 0.5, 1.25), c(0.2, -0.1, 0.4))
  offset <- matrix(1e8, nrow = nrow(x), ncol = ncol(x))
  length_scale <- c(0.7, 0.3)

  evaluators <- list(
    kernel_rbf,
    kernel_matern12,
    kernel_matern32,
    kernel_matern52,
    kernel_rational_quadratic
  )

  # Adding 1e8 rounds inputs to about 1.5e-8, so only that much is lost.
  for (evaluate in evaluators) {
    expect_equal(
      evaluate(x + offset, length_scale = length_scale),
      evaluate(x, length_scale = length_scale),
      tolerance = 1e-6
    )
  }

  expect_equal(
    kernel_rbf(1e8 + c(0, 0.5))[1, 2],
    exp(-0.125),
    tolerance = 1e-6
  )
})


test_that("kernel parameters are validated consistently", {
  x <- c(0, 1)

  expect_error(
    kernel_rbf(x, variance = 0),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    kernel_matern32(x, length_scale = -1),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    kernel_rational_quadratic(x, alpha = 0),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    kernel_periodic(x, period = 0),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("kernel inputs require matching feature dimensions", {
  x <- matrix(1:6, ncol = 2)
  y <- matrix(1:9, ncol = 3)

  expect_error(
    kernel_rbf(x, y),
    class = "gaussianprocesses_validation_error"
  )
})
