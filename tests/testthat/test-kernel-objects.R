test_that("leaf kernel specifications reproduce direct kernels", {
  x <- c(0, 0.5, 1)

  specification <- rbf_kernel(variance = 1.7, length_scale = 0.4)

  expect_equal(
    evaluate_kernel(specification, x),
    kernel_rbf(x, variance = 1.7, length_scale = 0.4),
    tolerance = 1e-12
  )
})


test_that("sum kernels add child covariance matrices", {
  x <- c(0, 1, 2)
  first <- rbf_kernel(variance = 2, length_scale = 0.5)
  second <- linear_kernel(variance = 0.25)

  specification <- sum_kernel(first, second)

  expect_equal(
    evaluate_kernel(specification, x),
    evaluate_kernel(first, x) + evaluate_kernel(second, x),
    tolerance = 1e-12
  )
})


test_that("product kernels multiply child covariance matrices element-wise", {
  x <- c(0, 1, 2)
  first <- matern32_kernel(variance = 1.2, length_scale = 0.8)
  second <- periodic_kernel(
    variance = 0.7,
    length_scale = 1.5,
    period = 2
  )

  specification <- product_kernel(first, second)

  expect_equal(
    evaluate_kernel(specification, x),
    evaluate_kernel(first, x) * evaluate_kernel(second, x),
    tolerance = 1e-12
  )
})


test_that("scale kernels retain a tunable scalar multiplier", {
  x <- c(0, 1)
  child <- matern52_kernel(variance = 1.1, length_scale = 2)
  specification <- scale_kernel(child, scale = 3)

  expect_equal(
    evaluate_kernel(specification, x),
    3 * evaluate_kernel(child, x),
    tolerance = 1e-12
  )

  expect_equal(
    kernel_parameters(specification, flatten = TRUE),
    c(
      scale = 3,
      kernel.variance = 1.1,
      kernel.length_scale = 2
    )
  )
})


test_that("nested composite kernels evaluate recursively", {
  x <- c(0, 0.5, 1.5)

  trend <- linear_kernel(variance = 0.2)
  local <- product_kernel(
    rbf_kernel(variance = 1.5, length_scale = 0.6),
    periodic_kernel(variance = 0.8, length_scale = 1.2, period = 1)
  )
  specification <- sum_kernel(trend, scale_kernel(local, scale = 0.5))

  expected <- evaluate_kernel(trend, x) +
    0.5 * evaluate_kernel(local, x)

  expect_equal(
    evaluate_kernel(specification, x),
    expected,
    tolerance = 1e-12
  )
})


test_that("parameter extraction preserves nested structure", {
  specification <- sum_kernel(
    rbf_kernel(variance = 2, length_scale = 0.3),
    scale_kernel(
      periodic_kernel(variance = 1, length_scale = 0.8, period = 4),
      scale = 0.5
    )
  )

  nested <- kernel_parameters(specification)
  flat <- kernel_parameters(specification, flatten = TRUE)

  expect_equal(nested$kernel1$variance, 2)
  expect_equal(nested$kernel2$scale, 0.5)
  expect_equal(nested$kernel2$kernel$period, 4)

  expect_equal(
    names(flat),
    c(
      "kernel1.variance",
      "kernel1.length_scale",
      "kernel2.scale",
      "kernel2.kernel.variance",
      "kernel2.kernel.length_scale",
      "kernel2.kernel.period"
    )
  )
})


test_that("flat parameters support complete optimization-style round trips", {
  specification <- product_kernel(
    matern32_kernel(variance = 1.5, length_scale = 2),
    rational_quadratic_kernel(
      variance = 0.75,
      length_scale = 1.2,
      alpha = 0.9
    )
  )

  parameters <- kernel_parameters(specification, flatten = TRUE)
  optimized <- parameters * 1.1
  updated <- update_kernel_parameters(specification, optimized)

  expect_equal(
    kernel_parameters(updated, flatten = TRUE),
    optimized
  )
})


test_that("parameter updates can target a nested subset", {
  specification <- sum_kernel(
    rbf_kernel(variance = 1, length_scale = 2),
    white_noise_kernel(variance = 0.1)
  )

  updated <- update_kernel_parameters(
    specification,
    c(
      kernel1.length_scale = 3.5,
      kernel2.variance = 0.25
    )
  )

  expect_equal(
    kernel_parameters(updated, flatten = TRUE),
    c(
      kernel1.variance = 1,
      kernel1.length_scale = 3.5,
      kernel2.variance = 0.25
    )
  )

  expect_equal(
    kernel_parameters(specification, flatten = TRUE),
    c(
      kernel1.variance = 1,
      kernel1.length_scale = 2,
      kernel2.variance = 0.1
    )
  )
})


test_that("invalid composite kernels and updates fail clearly", {
  expect_error(
    sum_kernel(rbf_kernel()),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    product_kernel(rbf_kernel(), "not a kernel"),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    scale_kernel(rbf_kernel(), scale = 0),
    class = "gaussianprocesses_validation_error"
  )

  specification <- rbf_kernel()

  expect_error(
    update_kernel_parameters(
      specification,
      c(unknown = 2)
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    update_kernel_parameters(
      specification,
      c(length_scale = -1)
    ),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("kernel specifications have a readable recursive print method", {
  specification <- sum_kernel(
    rbf_kernel(variance = 2, length_scale = 0.5),
    scale_kernel(
      white_noise_kernel(variance = 0.1),
      scale = 0.25
    )
  )

  output <- capture.output(print(specification))

  expect_true(any(grepl("SumKernel", output, fixed = TRUE)))
  expect_true(any(grepl("RBF", output, fixed = TRUE)))
  expect_true(any(grepl("ScaleKernel", output, fixed = TRUE)))
  expect_true(any(grepl("WhiteNoise", output, fixed = TRUE)))
})
