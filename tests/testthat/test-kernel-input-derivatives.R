# Covariance derivatives with respect to the inputs.

# Inputs whose paired rows are a distance r apart: y = x + r v for a unit
# vector v, so K(x, y) has diagonal entries at distance r and off-diagonal
# entries at ordinary distances.
paired_inputs <- function(dimension, distance) {
  x <- cbind(
    c(-1.1, -0.3, 0.4, 1.3),
    c(0.6, -0.8, 0.2, 1.1),
    c(1.5, 0.3, -0.6, -1.2)
  )[, seq_len(dimension), drop = FALSE]
  direction <- c(0.6, -0.48, 0.64)[seq_len(dimension)]
  direction <- direction / sqrt(sum(direction^2))

  list(x = x, y = sweep(x, 2, distance * direction, "+"))
}

input_derivative_distances <- c(0, 1e-8, 1e-4, 1)

smooth_kernel_cases <- function(dimension) {
  length_scale <- if (dimension == 1L) 0.8 else c(0.7, 1.3, 1.9)

  list(
    rbf = rbf_kernel(variance = 1.4, length_scale = length_scale),
    matern32 = matern32_kernel(variance = 0.9, length_scale = length_scale),
    matern52 = matern52_kernel(variance = 1.2, length_scale = length_scale),
    rational_quadratic = rational_quadratic_kernel(
      variance = 0.8,
      length_scale = length_scale,
      alpha = 0.7
    ),
    periodic = periodic_kernel(
      variance = 1.1,
      length_scale = length_scale,
      period = 1.7
    ),
    linear = linear_kernel(variance = 0.6),
    sum = sum_kernel(
      rbf_kernel(variance = 0.5, length_scale = length_scale),
      periodic_kernel(length_scale = 1.2, period = 2.3)
    ),
    product = product_kernel(
      matern52_kernel(length_scale = length_scale),
      linear_kernel(variance = 0.4)
    ),
    scale = scale_kernel(
      rational_quadratic_kernel(length_scale = 1.1, alpha = 2),
      2.5
    ),
    nested = product_kernel(
      sum_kernel(
        matern32_kernel(length_scale = length_scale),
        scale_kernel(linear_kernel(), 0.3)
      ),
      periodic_kernel(period = 3, length_scale = 0.9)
    )
  )
}


test_that("input derivatives match central differences at small and large r", {
  for (dimension in c(1L, 3L)) {
    kernels <- smooth_kernel_cases(dimension)

    for (distance in input_derivative_distances) {
      inputs <- paired_inputs(dimension, distance)

      for (name in names(kernels)) {
        expect_input_derivative_contract(
          kernels[[name]],
          inputs$x,
          inputs$y,
          sprintf("%s, D = %d, r = %g", name, dimension, distance)
        )
      }
    }
  }
})


test_that("hyperparameter derivatives of input derivatives are correct", {
  for (dimension in c(1L, 3L)) {
    kernels <- smooth_kernel_cases(dimension)
    kernels$additive <- sum_kernel(
      select_dimensions(rbf_kernel(length_scale = 0.7), 1),
      select_dimensions(periodic_kernel(period = 2), dimension)
    )

    for (distance in c(0, 1)) {
      inputs <- paired_inputs(dimension, distance)

      for (name in names(kernels)) {
        expect_input_parameter_gradient_contract(
          kernels[[name]],
          inputs$x,
          inputs$y,
          sprintf("%s, D = %d, r = %g", name, dimension, distance)
        )
      }
    }
  }
})


test_that("radial derivatives at coincident inputs take their limits", {
  length_scale <- c(0.7, 1.3, 1.9)
  x <- matrix(c(0.2, -0.5, 1), nrow = 1)
  # -kappa''(0) for kappa(r) the correlation as a function of r.
  curvature <- list(
    rbf = 1,
    matern32 = 3,
    matern52 = 5 / 3,
    rational_quadratic = 1
  )
  constructors <- list(
    rbf = rbf_kernel,
    matern32 = matern32_kernel,
    matern52 = matern52_kernel,
    rational_quadratic = rational_quadratic_kernel
  )

  for (name in names(constructors)) {
    kernel <- constructors[[name]](variance = 1.5, length_scale = length_scale)
    derivatives <- gaussianprocesses:::.kernel_input_derivatives(
      kernel,
      x,
      x,
      order = 2L
    )

    expect_identical(
      unlist(derivatives$gradient_x),
      c(0, 0, 0),
      info = name
    )

    for (d in 1:3) {
      for (e in 1:3) {
        expected <- 0

        if (d == e) {
          expected <- 1.5 * curvature[[name]] / length_scale[d]^2
        }

        expect_equal(
          drop(derivatives$cross[[d]][[e]]),
          expected,
          tolerance = 1e-14,
          info = sprintf("%s, d = %d, e = %d", name, d, e)
        )
      }
    }
  }
})


test_that("kernel_input_gradient() returns one matrix per input dimension", {
  x <- cbind(c(0, 0.5, 1.5), c(1, -1, 0.25))
  y <- x[1:2, ] + 0.1
  kernel <- rbf_kernel(length_scale = c(0.5, 2))

  gradient <- kernel_input_gradient(kernel, x, y)

  expect_length(gradient, 2L)
  expect_identical(dim(gradient[[1]]), c(3L, 2L))
  expect_identical(
    kernel_input_gradient(kernel, x),
    kernel_input_gradient(kernel, x, x)
  )
  expect_identical(
    kernel_input_gradient(rbf_kernel(), c(0, 1)),
    kernel_input_gradient(rbf_kernel(), matrix(c(0, 1)))
  )
  # A stationary kernel's derivative in x is minus its derivative in y.
  expect_equal(
    gradient[[2]],
    -gaussianprocesses:::.kernel_input_derivatives(kernel, x, y)$gradient_y[[2]]
  )
  expect_error(
    kernel_input_gradient(kernel, x, matrix(1:3, 1)),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    kernel_input_gradient(list(), x),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("values and gradients have a positive semidefinite covariance", {
  for (dimension in c(1L, 3L)) {
    x <- rbind(
      paired_inputs(dimension, 0.3)$x,
      paired_inputs(dimension, 0.3)$y,
      # A repeated input makes the joint covariance singular.
      paired_inputs(dimension, 0)$x[1, , drop = FALSE]
    )
    kernels <- smooth_kernel_cases(dimension)

    for (name in names(kernels)) {
      covariance <- joint_value_gradient_covariance(kernels[[name]], x)
      scale <- max(abs(diag(covariance)))
      label <- sprintf("%s, D = %d", name, dimension)

      expect_equal(covariance, t(covariance), tolerance = 1e-12, info = label)
      expect_gte(
        min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values),
        -1e-10 * scale,
        label = paste(label, "smallest eigenvalue")
      )
    }
  }
})


test_that("kernels that are not mean-square differentiable are rejected", {
  x <- cbind(c(0, 0.5, 1), c(1, 0, -1))
  rough <- list(
    matern12 = matern12_kernel(),
    white_noise = white_noise_kernel(),
    sum = sum_kernel(rbf_kernel(), white_noise_kernel(variance = 0.1)),
    product = product_kernel(matern12_kernel(), periodic_kernel()),
    scale = scale_kernel(matern12_kernel(), 2),
    nested = sum_kernel(
      product_kernel(rbf_kernel(), matern12_kernel()),
      white_noise_kernel()
    )
  )

  for (name in names(rough)) {
    expect_error(
      kernel_input_gradient(rough[[name]], x),
      class = "gaussianprocesses_smoothness_error",
      info = name
    )
  }

  expect_error(
    kernel_input_gradient(rough$nested, x),
    "Matern-1/2, WhiteNoise",
    fixed = TRUE
  )
  expect_error(
    kernel_input_gradient(rough$sum, x),
    class = "gaussianprocesses_validation_error"
  )

  # A once-differentiable component is enough.
  expect_identical(
    gaussianprocesses:::.kernel_smoothness(
      sum_kernel(matern32_kernel(), matern52_kernel())
    ),
    1
  )
  expect_length(
    kernel_input_gradient(sum_kernel(matern32_kernel(), rbf_kernel()), x),
    2L
  )
})
