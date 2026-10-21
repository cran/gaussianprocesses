test_that("analytical derivatives match central finite differences", {
  for (case in kernel_gradient_cases()) {
    analytical <- kernel_gradient(case$kernel, case$x, case$y)
    numerical <- finite_difference_kernel_gradient(
      case$kernel,
      case$x,
      case$y
    )

    expect_identical(names(analytical), names(numerical), info = case$name)

    for (path in names(analytical)) {
      expect_lte(
        scaled_gradient_discrepancy(analytical[[path]], numerical[[path]]),
        gradient_fd_tolerance,
        label = paste0(case$name, ": ", path)
      )
    }
  }
})


test_that("finite differences converge to the analytical value at second order", {
  # Central differences have O(h^2) error, so a tenfold smaller step should
  # shrink the discrepancy about a hundredfold; a wrong derivative would not.
  cases <- kernel_gradient_cases()
  names(cases) <- vapply(cases, `[[`, character(1), "name")

  for (name in c(
    "rbf 2D ARD",
    "matern12 1D",
    "rational quadratic small alpha",
    "periodic 1D cross",
    "product of three cross",
    "nested scale sum product"
  )) {
    case <- cases[[name]]
    analytical <- kernel_gradient(case$kernel, case$x, case$y)

    discrepancy <- function(step) {
      numerical <- finite_difference_kernel_gradient(
        case$kernel,
        case$x,
        case$y,
        step = step
      )
      max(mapply(scaled_gradient_discrepancy, analytical, numerical))
    }

    ratio <- discrepancy(1e-3) / discrepancy(1e-4)

    expect_gt(ratio, 50, label = name)
    expect_lt(ratio, 200, label = name)
  }
})


test_that("derivative names follow the flattened parameter paths", {
  kernel <- scale_kernel(
    sum_kernel(
      product_kernel(
        rbf_kernel(0.9, c(0.5, 1.5)),
        periodic_kernel(1.2, 0.7, period = 1.6)
      ),
      white_noise_kernel(0.05)
    ),
    scale = 0.7
  )
  x <- cbind(c(0, 1, 2), c(1, 0, 3))

  gradient <- kernel_gradient(kernel, x)

  expect_identical(
    names(gradient),
    c(
      "scale",
      "kernel.kernel1.kernel1.variance",
      "kernel.kernel1.kernel1.length_scale[1]",
      "kernel.kernel1.kernel1.length_scale[2]",
      "kernel.kernel1.kernel2.variance",
      "kernel.kernel1.kernel2.length_scale",
      "kernel.kernel1.kernel2.period",
      "kernel.kernel2.variance"
    )
  )
  expect_identical(
    names(gradient),
    names(kernel_parameters(kernel, flatten = TRUE))
  )

  for (case in kernel_gradient_cases()) {
    expect_identical(
      names(kernel_gradient(case$kernel, case$x, case$y)),
      names(kernel_parameters(case$kernel, flatten = TRUE)),
      info = case$name
    )
  }
})


test_that("derivatives have covariance dimensions and symmetry", {
  for (case in kernel_gradient_cases()) {
    gradient <- kernel_gradient(case$kernel, case$x, case$y)
    covariance <- evaluate_kernel(case$kernel, case$x, case$y)

    for (path in names(gradient)) {
      expect_identical(dim(gradient[[path]]), dim(covariance), info = path)
      expect_true(all(is.finite(gradient[[path]])), info = path)

      if (is.null(case$y)) {
        expect_equal(
          gradient[[path]],
          t(gradient[[path]]),
          tolerance = 1e-12,
          info = paste0(case$name, ": ", path)
        )
      }
    }
  }
})


test_that("omitting y matches passing the same inputs twice", {
  kernel <- sum_kernel(
    matern32_kernel(1.1, c(0.6, 1.4)),
    white_noise_kernel(0.2)
  )
  x <- cbind(c(0, 0.5, 0.5, 2), c(1, 0, 0, 3))

  expect_identical(kernel_gradient(kernel, x), kernel_gradient(kernel, x, x))
})


test_that("variance and scale derivatives equal the covariance they scale", {
  x <- cbind(c(-1, 0, 0.4, 2), c(0.5, 1, -0.3, 0))
  leaf <- periodic_kernel(1.7, c(0.6, 1.2), period = 1.5)

  expect_identical(
    kernel_gradient(leaf, x)[["variance"]],
    evaluate_kernel(leaf, x)
  )

  scaled <- scale_kernel(
    product_kernel(leaf, rbf_kernel(0.5, 2)),
    scale = 3
  )

  expect_identical(
    kernel_gradient(scaled, x)[["scale"]],
    evaluate_kernel(scaled, x)
  )
})


test_that("a shared length scale differentiates as the sum over dimensions", {
  x <- cbind(c(-1, 0, 0.4, 2), c(0.5, 1, -0.3, 0), c(2, 1, 0, -1))

  for (constructor in list(
    rbf_kernel,
    matern12_kernel,
    matern32_kernel,
    matern52_kernel
  )) {
    shared <- kernel_gradient(constructor(1.2, 0.9), x)
    ard <- kernel_gradient(constructor(1.2, c(0.9, 0.9, 0.9)), x)

    expect_equal(
      shared[["length_scale"]],
      ard[["length_scale[1]"]] + ard[["length_scale[2]"]] +
        ard[["length_scale[3]"]],
      tolerance = 1e-12
    )
  }
})


test_that("gradient evaluation reproduces evaluate_kernel exactly", {
  for (case in kernel_gradient_cases()) {
    inputs <- gaussianprocesses:::.prepare_kernel_inputs(case$x, case$y)
    result <- gaussianprocesses:::.kernel_value_and_gradient(
      case$kernel,
      inputs$x,
      inputs$y
    )

    expect_identical(
      result$value,
      evaluate_kernel(case$kernel, case$x, case$y),
      info = case$name
    )
  }
})


test_that("kernel_gradient validates its inputs", {
  expect_error(
    kernel_gradient(list(type = "rbf"), 1:3),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    kernel_gradient(rbf_kernel(length_scale = c(1, 2)), 1:3),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    kernel_gradient(rbf_kernel(), matrix(1:4, 2), 1:3),
    class = "gaussianprocesses_validation_error"
  )
})
