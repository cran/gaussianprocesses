# Generates kernel-golden.rds, the reference values for test-kernel-golden.R.
#
# The committed fixture was generated from commit 84de3c0, before the kernel
# registry refactor, so the tests show that the refactor preserved every
# covariance, diagonal, gradient, parameter path, and printed form. Rerun this
# script only for an intentional change to kernel behaviour, and say so in
# NEWS.md. From the repository root:
#
#   Rscript -e 'pkgload::load_all(); source("tests/testthat/fixtures/make-kernel-golden.R")'

x_1d <- c(-1.2, -0.4, 0, 0, 0.35, 1.1, 2.4)
y_1d <- c(-0.4, 0.9, 3)

x_2d <- cbind(
  c(-1.2, -0.4, 0, 0, 0.35, 1.1),
  c(0.5, -0.3, 1.2, 1.2, 0, -0.8)
)
y_2d <- cbind(c(-0.4, 0.9, 3), c(-0.3, 0.1, 1))

x_3d <- cbind(x_2d, c(2, 0.1, -0.7, -0.7, 1.4, 0.3))
y_3d <- cbind(y_2d, c(0.1, -1, 0.6))

cases <- list(
  rbf_scalar = list(
    call = quote(rbf_kernel(variance = 1.3, length_scale = 0.7)),
    x = x_1d, y = y_1d
  ),
  rbf_ard = list(
    call = quote(rbf_kernel(variance = 0.8, length_scale = c(0.6, 1.4, 2.1))),
    x = x_3d, y = y_3d
  ),
  matern12_scalar = list(
    call = quote(matern12_kernel(variance = 1.1, length_scale = 0.9)),
    x = x_1d, y = y_1d
  ),
  matern12_ard = list(
    call = quote(matern12_kernel(variance = 0.7, length_scale = c(0.5, 1.3, 2))),
    x = x_3d, y = y_3d
  ),
  matern32_scalar = list(
    call = quote(matern32_kernel(variance = 1.4, length_scale = 0.8)),
    x = x_1d, y = y_1d
  ),
  matern32_ard = list(
    call = quote(matern32_kernel(variance = 0.9, length_scale = c(1.1, 0.4))),
    x = x_2d, y = y_2d
  ),
  matern52_scalar = list(
    call = quote(matern52_kernel(variance = 0.6, length_scale = 1.2)),
    x = x_1d, y = y_1d
  ),
  matern52_ard = list(
    call = quote(matern52_kernel(variance = 1.2, length_scale = c(0.7, 1.6, 0.9))),
    x = x_3d, y = y_3d
  ),
  rational_quadratic_scalar = list(
    call = quote(rational_quadratic_kernel(1.5, length_scale = 0.6, alpha = 0.8)),
    x = x_1d, y = y_1d
  ),
  rational_quadratic_ard = list(
    call = quote(rational_quadratic_kernel(0.9, length_scale = c(0.9, 1.7), alpha = 2.5)),
    x = x_2d, y = y_2d
  ),
  periodic_scalar = list(
    call = quote(periodic_kernel(variance = 1.1, length_scale = 0.8, period = 1.3)),
    x = x_1d, y = y_1d
  ),
  periodic_ard = list(
    call = quote(periodic_kernel(variance = 0.5, length_scale = c(0.6, 1.5), period = 2.2)),
    x = x_2d, y = y_2d
  ),
  linear = list(
    call = quote(linear_kernel(variance = 0.4)),
    x = x_2d, y = y_2d
  ),
  white_noise = list(
    call = quote(white_noise_kernel(variance = 0.05)),
    x = x_2d, y = y_2d
  ),
  sum = list(
    call = quote(sum_kernel(
      rbf_kernel(1, c(0.8, 1.9)),
      white_noise_kernel(0.1)
    )),
    x = x_2d, y = y_2d
  ),
  product = list(
    call = quote(product_kernel(
      periodic_kernel(0.9, 1.1, period = 1.7),
      matern32_kernel(1.2, c(0.5, 2.2)),
      linear_kernel(0.3)
    )),
    x = x_2d, y = y_2d
  ),
  scale = list(
    call = quote(scale_kernel(
      rational_quadratic_kernel(0.7, c(1.3, 0.5), alpha = 1.6),
      scale = 2.5
    )),
    x = x_2d, y = y_2d
  ),
  nested = list(
    call = quote(scale_kernel(
      sum_kernel(
        product_kernel(
          rbf_kernel(1.1, c(0.7, 1.8, 1.2)),
          periodic_kernel(0.8, c(1, 0.6, 2), period = 1.4)
        ),
        white_noise_kernel(0.02),
        linear_kernel(0.2)
      ),
      scale = 0.6
    )),
    x = x_3d, y = y_3d
  ),
  deeply_nested = list(
    call = quote(sum_kernel(
      scale_kernel(
        product_kernel(matern52_kernel(0.9, 0.8), matern12_kernel(1.3, 1.5)),
        scale = 1.7
      ),
      sum_kernel(
        rbf_kernel(0.5, 2),
        rational_quadratic_kernel(0.4, 0.9, alpha = 3)
      )
    )),
    x = x_1d, y = y_1d
  )
)

golden <- lapply(cases, function(case) {
  kernel <- eval(case$call)
  flat <- kernel_parameters(kernel, flatten = TRUE)

  # Change every parameter, then only the last one, to cover full and partial
  # updates.
  updated_all <- update_kernel_parameters(kernel, flat * 1.1)
  updated_last <- update_kernel_parameters(kernel, utils::tail(flat, 1) * 0.5)

  c(
    case,
    list(
      kernel = kernel,
      covariance = evaluate_kernel(kernel, case$x),
      cross_covariance = evaluate_kernel(kernel, case$x, case$y),
      diagonal = kernel_diagonal(kernel, case$x),
      gradient = kernel_gradient(kernel, case$x),
      cross_gradient = kernel_gradient(kernel, case$x, case$y),
      flat_parameters = flat,
      nested_parameters = kernel_parameters(kernel),
      updated_all = updated_all,
      updated_last = updated_last,
      printed = utils::capture.output(print(kernel))
    )
  )
})

saveRDS(
  golden,
  file.path("tests", "testthat", "fixtures", "kernel-golden.rds"),
  version = 2
)
