# kernel-golden.rds was generated before the kernel registry refactor (see
# fixtures/make-kernel-golden.R). These tests show that kernels built from the
# registry reproduce the earlier covariances, diagonals, gradients, parameter
# paths, and printed forms.

golden <- readRDS(test_path("fixtures", "kernel-golden.rds"))

# Largest absolute difference relative to the largest reference magnitude.
expect_matches_golden <- function(actual, expected, label) {
  scale <- max(abs(expected), .Machine$double.xmin)

  expect_lte(
    max(abs(actual - expected)) / scale,
    1e-14,
    label = label
  )
}


test_that("constructors build the same kernel objects", {
  for (name in names(golden)) {
    case <- golden[[name]]
    expect_identical(eval(case$call), case$kernel, info = name)
  }
})


test_that("covariances and diagonals match the golden values", {
  for (name in names(golden)) {
    case <- golden[[name]]
    kernel <- eval(case$call)

    expect_matches_golden(
      evaluate_kernel(kernel, case$x),
      case$covariance,
      paste(name, "covariance")
    )
    expect_matches_golden(
      evaluate_kernel(kernel, case$x, case$y),
      case$cross_covariance,
      paste(name, "cross-covariance")
    )
    expect_matches_golden(
      kernel_diagonal(kernel, case$x),
      case$diagonal,
      paste(name, "diagonal")
    )
  }
})


test_that("hyperparameter gradients match the golden values", {
  for (name in names(golden)) {
    case <- golden[[name]]
    kernel <- eval(case$call)

    for (which in c("gradient", "cross_gradient")) {
      gradient <- if (identical(which, "gradient")) {
        kernel_gradient(kernel, case$x)
      } else {
        kernel_gradient(kernel, case$x, case$y)
      }

      expect_identical(names(gradient), names(case[[which]]), info = name)

      for (path in names(gradient)) {
        expect_matches_golden(
          gradient[[path]],
          case[[which]][[path]],
          paste(name, which, path)
        )
      }
    }
  }
})


test_that("parameter paths, updates, and printing are unchanged", {
  for (name in names(golden)) {
    case <- golden[[name]]
    kernel <- eval(case$call)
    flat <- kernel_parameters(kernel, flatten = TRUE)

    expect_identical(flat, case$flat_parameters, info = name)
    expect_identical(
      kernel_parameters(kernel),
      case$nested_parameters,
      info = name
    )
    expect_identical(
      update_kernel_parameters(kernel, flat * 1.1),
      case$updated_all,
      info = name
    )
    expect_identical(
      update_kernel_parameters(kernel, utils::tail(flat, 1) * 0.5),
      case$updated_last,
      info = name
    )
    expect_identical(
      utils::capture.output(print(kernel)),
      case$printed,
      info = name
    )
  }
})
