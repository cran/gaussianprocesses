test_that("zero mean evaluates to zero for every observation", {
  x <- matrix(1:6, ncol = 2)

  expect_equal(
    evaluate_mean(zero_mean(), x),
    rep(0, nrow(x))
  )
})


test_that("constant mean evaluates consistently", {
  x <- c(0, 1, 2, 3)

  expect_equal(
    evaluate_mean(constant_mean(2.5), x),
    rep(2.5, length(x))
  )
})


test_that("constant mean validates its value", {
  expect_error(
    constant_mean(NA_real_),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    constant_mean(Inf),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("mean specifications are validated and printed", {
  expect_error(
    fit_gp(c(0, 1, 2), c(0, 1, 0), rbf_kernel(), noise_variance = 0.1, mean = 1),
    "'mean' must be a Gaussian-process mean specification",
    class = "gaussianprocesses_validation_error"
  )

  unknown <- structure(
    list(type = "quadratic", parameters = list()),
    class = "gaussianprocesses_mean"
  )
  expect_error(
    evaluate_mean(unknown, 1:3),
    "Unknown mean-function type 'quadratic'",
    class = "gaussianprocesses_validation_error"
  )
  expect_identical(gaussianprocesses:::.format_mean(unknown), "quadratic")

  expect_identical(
    gaussianprocesses:::.format_mean(constant_mean(2.5)),
    "ConstantMean(value=2.5)"
  )
  model <- fit_gp(
    c(0, 1, 2),
    c(2, 3, 2),
    rbf_kernel(),
    noise_variance = 0.1,
    mean = constant_mean(2)
  )
  expect_true(
    any(grepl("ConstantMean(value=2)", capture.output(print(model)), fixed = TRUE))
  )
})


test_that("parametric means evaluate their basis with fixed coefficients", {
  x <- cbind(c(0, 1, 2), c(-1, 0, 3))

  expect_equal(
    evaluate_mean(linear_mean(c(1, 2, -1)), x),
    1 + 2 * x[, 1] - x[, 2]
  )
  expect_equal(
    evaluate_mean(constant_mean(coefficient_prior(3, 1)), x),
    rep(3, 3)
  )

  quadratic <- polynomial_mean(2, c(1, 0.5, 2), center = 1, scale = 2)
  z <- (c(0, 1, 2) - 1) / 2
  expect_equal(evaluate_mean(quadratic, c(0, 1, 2)), 1 + 0.5 * z + 2 * z^2)

  two_inputs <- polynomial_mean(
    1,
    c(1, 2, 3),
    center = c(0, 1),
    scale = c(1, 2)
  )
  expect_equal(
    evaluate_mean(two_inputs, x),
    1 + 2 * x[, 1] + 3 * (x[, 2] - 1) / 2
  )
  expect_equal(
    evaluate_mean(polynomial_mean(0, 4, center = 0, scale = 1), x),
    rep(4, 3)
  )

  seasonal <- basis_mean(function(x) sin(x[, 1]), coefficients = 2)
  expect_equal(evaluate_mean(seasonal, c(0, pi / 2)), c(0, 2))

  # A Gaussian prior evaluates to its prior mean.
  expect_equal(
    evaluate_mean(linear_mean(coefficient_prior(c(1, 2, 0), 1)), x),
    1 + 2 * x[, 1]
  )
})


test_that("basis names follow the basis columns", {
  x <- cbind(c(0, 1, 2), c(-1, 0, 3))
  basis <- function(mean) gaussianprocesses:::.mean_basis(mean, x)

  expect_identical(colnames(basis(linear_mean())), c("intercept", "x1", "x2"))
  expect_identical(
    colnames(basis(polynomial_mean(3, center = 0, scale = 1))),
    c("intercept", "x1", "x1^2", "x1^3", "x2", "x2^2", "x2^3")
  )
  expect_identical(
    colnames(basis(polynomial_mean(1, center = 0, scale = 1))),
    c("intercept", "x1", "x2")
  )
  expect_identical(
    colnames(basis(polynomial_mean(0, center = 0, scale = 1))),
    "intercept"
  )
  expect_identical(
    colnames(basis(basis_mean(function(x) cbind(a = 1, b = x[, 2])))),
    c("a", "b")
  )
  expect_identical(
    colnames(basis(basis_mean(function(x) cbind(1, x[, 2])))),
    c("h1", "h2")
  )
})


test_that("polynomial means are standardized with the training inputs", {
  x <- cbind(seq(10, 20, length.out = 15), -1 + 2 * sin(1:15))
  y <- 0.1 * x[, 1]^2 + x[, 2]
  model <- fit_gp(x, y, rbf_kernel(length_scale = c(3, 2)),
                  noise_variance = 0.01, mean = polynomial_mean(2))
  ranges <- range(x[, 2])

  expect_equal(model$mean$center, c(15, mean(ranges)))
  expect_equal(model$mean$scale, c(5, diff(ranges) / 2))
  basis <- model$mean_fit$basis
  expect_true(all(abs(basis[, -1]) <= 1 + 1e-12))

  # Supplied values are kept, and inputs without range keep scale 1.
  constant_column <- cbind(x[, 1], 2)
  resolved <- gaussianprocesses:::.resolve_mean(
    polynomial_mean(1, center = c(0, 0)),
    constant_column
  )
  expect_identical(resolved$center, c(0, 0))
  expect_identical(resolved$scale, c(5, 1))
  expect_identical(
    gaussianprocesses:::.resolve_mean(linear_mean(), x),
    linear_mean()
  )

  expect_error(
    evaluate_mean(polynomial_mean(2, c(1, 2, 3)), 1:3),
    "no 'center' and 'scale'",
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    evaluate_mean(polynomial_mean(1, c(1, 2, 3), center = c(0, 0, 0),
                                  scale = 1), x),
    "one value or one per input column",
    class = "gaussianprocesses_validation_error"
  )
})


test_that("mean constructors validate their arguments", {
  invalid <- list(
    function() linear_mean("a"),
    function() linear_mean(c(1, NA)),
    function() linear_mean(numeric(0)),
    function() linear_mean(list(treatment = "estimated")),
    function() polynomial_mean(-1),
    function() polynomial_mean(1.5),
    function() polynomial_mean(c(1, 2)),
    function() polynomial_mean(2, center = "a"),
    function() polynomial_mean(2, center = Inf),
    function() polynomial_mean(2, scale = 0),
    function() polynomial_mean(2, scale = numeric(0)),
    function() basis_mean(1:3),
    function() constant_mean("a"),
    function() coefficient_prior(mean = NA),
    function() coefficient_prior(mean = matrix(1)),
    function() coefficient_prior(covariance = 0),
    function() coefficient_prior(covariance = c(1, -1)),
    function() coefficient_prior(covariance = "a"),
    function() coefficient_prior(covariance = matrix(c(1, 2, 0, 1), 2)),
    function() coefficient_prior(covariance = matrix(c(1, 2, 2, 1), 2)),
    function() coefficient_prior(covariance = matrix(1, 2, 3))
  )

  for (index in seq_along(invalid)) {
    expect_error(
      invalid[[index]](),
      class = "gaussianprocesses_validation_error",
      info = index
    )
  }

  expect_error(estimate_coefficients("bayes"))
  expect_s3_class(
    coefficient_prior(1, diag(2)),
    "gaussianprocesses_coefficients"
  )
})


test_that("coefficients and bases are checked against each other", {
  x <- cbind(c(0, 1, 2, 3), c(1, 0, 1, 0))
  y <- c(1, 2, 2, 4)

  expect_error(
    evaluate_mean(linear_mean(c(1, 2)), x),
    "3 basis column(s) (intercept, x1, x2) but 2 fixed coefficient(s)",
    fixed = TRUE,
    class = "gaussianprocesses_validation_error"
  )

  for (basis in list(
    function(x) "a",
    function(x) matrix(1, nrow(x) + 1, 1),
    function(x) matrix(0, nrow(x), 0),
    function(x) cbind(1, NA)[rep(1, nrow(x)), ]
  )) {
    expect_error(
      evaluate_mean(basis_mean(basis, 1), x),
      "finite numeric matrix",
      class = "gaussianprocesses_validation_error"
    )
  }

  for (prior in list(
    coefficient_prior(mean = c(1, 2)),
    coefficient_prior(covariance = c(1, 2)),
    coefficient_prior(covariance = diag(2))
  )) {
    expect_error(
      fit_gp(x, y, rbf_kernel(), noise_variance = 0.1,
             mean = linear_mean(prior)),
      "one per coefficient|must be 3 x 3",
      class = "gaussianprocesses_validation_error"
    )
  }

  expect_error(
    evaluate_mean(linear_mean(), x),
    "no value before a model is fitted",
    class = "gaussianprocesses_validation_error"
  )
})


test_that("every kind of mean is formatted and printed", {
  formats <- list(
    "ZeroMean()" = zero_mean(),
    "ConstantMean(value=2.5)" = constant_mean(2.5),
    "ConstantMean(estimated by ML)" = constant_mean(estimate_coefficients()),
    "LinearMean(estimated by REML)" =
      linear_mean(estimate_coefficients("reml")),
    "LinearMean(coefficients=c(1, 2))" = linear_mean(c(1, 2)),
    "PolynomialMean(degree=3, marginalized, vague prior)" =
      polynomial_mean(3, coefficient_prior()),
    "BasisMean(marginalized, Gaussian prior)" =
      basis_mean(identity, coefficient_prior(0, 1)),
    "BasisMean(coefficients=2)" = basis_mean(identity, 2)
  )

  for (expected in names(formats)) {
    expect_identical(
      gaussianprocesses:::.format_mean(formats[[expected]]),
      expected
    )
  }

  expect_output(print(linear_mean()), "LinearMean(estimated by ML)",
                fixed = TRUE)
  utils::capture.output(expect_invisible(print(zero_mean())))
})
