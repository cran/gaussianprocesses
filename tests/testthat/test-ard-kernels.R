test_that("scalar length scales remain isotropic in multiple dimensions", {
  x <- matrix(
    c(
      0, 0,
      1, 0,
      0, 1
    ),
    ncol = 2,
    byrow = TRUE
  )

  covariance <- kernel_rbf(
    x,
    variance = 2,
    length_scale = 2
  )

  expected_12 <- 2 * exp(-0.5 * (1 / 2)^2)
  expected_13 <- 2 * exp(-0.5 * (1 / 2)^2)

  expect_equal(covariance[1, 2], expected_12, tolerance = 1e-12)
  expect_equal(covariance[1, 3], expected_13, tolerance = 1e-12)
})


test_that("ARD RBF applies one length scale per input dimension", {
  x <- matrix(
    c(
      0, 0,
      1, 0,
      0, 1
    ),
    ncol = 2,
    byrow = TRUE
  )

  covariance <- kernel_rbf(
    x,
    variance = 1,
    length_scale = c(0.5, 2)
  )

  expect_equal(
    covariance[1, 2],
    exp(-0.5 * (1 / 0.5)^2),
    tolerance = 1e-12
  )
  expect_equal(
    covariance[1, 3],
    exp(-0.5 * (1 / 2)^2),
    tolerance = 1e-12
  )
  expect_lt(covariance[1, 2], covariance[1, 3])
})


test_that("ARD Matérn kernels use scaled Euclidean distance", {
  x <- matrix(
    c(
      0, 0,
      1, 2
    ),
    ncol = 2,
    byrow = TRUE
  )
  length_scale <- c(2, 4)
  scaled_distance <- sqrt((1 / 2)^2 + (2 / 4)^2)

  matern12 <- kernel_matern12(x, length_scale = length_scale)
  matern32 <- kernel_matern32(x, length_scale = length_scale)
  matern52 <- kernel_matern52(x, length_scale = length_scale)

  expect_equal(
    matern12[1, 2],
    exp(-scaled_distance),
    tolerance = 1e-12
  )

  r32 <- sqrt(3) * scaled_distance
  expect_equal(
    matern32[1, 2],
    (1 + r32) * exp(-r32),
    tolerance = 1e-12
  )

  r52 <- sqrt(5) * scaled_distance
  expect_equal(
    matern52[1, 2],
    (1 + r52 + 5 * scaled_distance^2 / 3) * exp(-r52),
    tolerance = 1e-12
  )
})


test_that("ARD rational quadratic uses dimension-scaled distance", {
  x <- matrix(c(0, 0, 1, 2), ncol = 2, byrow = TRUE)
  length_scale <- c(2, 4)
  alpha <- 1.5
  scaled_squared_distance <- (1 / 2)^2 + (2 / 4)^2

  covariance <- kernel_rational_quadratic(
    x,
    length_scale = length_scale,
    alpha = alpha
  )

  expect_equal(
    covariance[1, 2],
    (1 + scaled_squared_distance / (2 * alpha))^(-alpha),
    tolerance = 1e-12
  )
})


test_that("periodic ARD separates dimension-specific smoothness", {
  x <- matrix(
    c(
      0, 0,
      0.5, 0,
      0, 0.5
    ),
    ncol = 2,
    byrow = TRUE
  )

  covariance <- kernel_periodic(
    x,
    length_scale = c(0.5, 2),
    period = 2
  )

  first_dimension_distance <- sin(pi * 0.5 / 2)^2 / 0.5^2
  second_dimension_distance <- sin(pi * 0.5 / 2)^2 / 2^2

  expect_equal(
    covariance[1, 2],
    exp(-2 * first_dimension_distance),
    tolerance = 1e-12
  )
  expect_equal(
    covariance[1, 3],
    exp(-2 * second_dimension_distance),
    tolerance = 1e-12
  )
  expect_lt(covariance[1, 2], covariance[1, 3])
})


test_that("length-scale dimension mismatches fail explicitly", {
  x <- matrix(1:8, ncol = 2)

  expect_error(
    kernel_rbf(x, length_scale = c(1, 2, 3)),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    kernel_matern32(x, length_scale = c(1, 2, 3)),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    kernel_periodic(x, length_scale = c(1, 2, 3)),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("ARD length-scale values must be finite and positive", {
  expect_error(
    rbf_kernel(length_scale = c(1, 0)),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    matern52_kernel(length_scale = c(1, Inf)),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("ARD parameters flatten to stable indexed paths", {
  kernel <- sum_kernel(
    rbf_kernel(variance = 2, length_scale = c(0.5, 3)),
    white_noise_kernel(variance = 0.1)
  )

  parameters <- kernel_parameters(kernel, flatten = TRUE)

  expect_equal(
    names(parameters),
    c(
      "kernel1.variance",
      "kernel1.length_scale[1]",
      "kernel1.length_scale[2]",
      "kernel2.variance"
    )
  )
  expect_equal(
    unname(parameters),
    c(2, 0.5, 3, 0.1)
  )
})


test_that("individual ARD dimensions can be updated independently", {
  kernel <- rbf_kernel(
    variance = 1,
    length_scale = c(0.5, 2)
  )

  updated <- update_kernel_parameters(
    kernel,
    c("length_scale[2]" = 4)
  )

  expect_equal(
    kernel_parameters(updated)$length_scale,
    c(0.5, 4)
  )
  expect_equal(
    kernel_parameters(kernel)$length_scale,
    c(0.5, 2)
  )
})


test_that("ARD kernel specifications print vector length scales clearly", {
  kernel <- rbf_kernel(
    variance = 1.5,
    length_scale = c(0.5, 2)
  )

  output <- capture.output(print(kernel))

  expect_true(
    any(
      grepl(
        "length_scale=c(0.5, 2)",
        output,
        fixed = TRUE
      )
    )
  )

  # Elements are formatted independently, not padded to shared decimals.
  expect_identical(
    gaussianprocesses:::.format_kernel_parameter_value(c(0.5, 2, 1e-10)),
    "c(0.5, 2, 1e-10)"
  )
  expect_identical(
    gaussianprocesses:::.format_kernel_parameter_value(1.5),
    "1.5"
  )
})


test_that("GP fitting and prediction accept ARD kernels", {
  x <- matrix(
    c(
      0, 0,
      1, 0,
      0, 1,
      1, 1
    ),
    ncol = 2,
    byrow = TRUE
  )
  y <- c(0, 1, 0.2, 1.1)

  model <- fit_gp(
    x,
    y,
    kernel = rbf_kernel(
      variance = 1,
      length_scale = c(0.4, 2)
    ),
    noise_variance = 0.05
  )

  prediction <- predict_gp(
    model,
    matrix(c(0.5, 0.25), nrow = 1)
  )

  expect_length(prediction$mean, 1L)
  expect_true(is.finite(prediction$mean))
  expect_gte(prediction$latent_variance, 0)
})


test_that("hyperparameter optimization supports ARD flattened parameters", {
  x <- matrix(
    c(
      0, 0,
      0.5, 0,
      1, 0,
      0, 1,
      0.5, 1,
      1, 1
    ),
    ncol = 2,
    byrow = TRUE
  )
  y <- x[, 1] + 0.05 * x[, 2]

  model <- optimize_gp(
    x,
    y,
    kernel = rbf_kernel(
      variance = 1,
      length_scale = c(1, 1)
    ),
    noise_variance = 0.05,
    optimize_noise = FALSE,
    n_starts = 1,
    control = list(maxit = 5)
  )

  paths <- names(
    kernel_parameters(model$kernel, flatten = TRUE)
  )

  expect_equal(
    paths,
    c(
      "variance",
      "length_scale[1]",
      "length_scale[2]"
    )
  )
})
