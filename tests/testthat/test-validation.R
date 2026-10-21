test_that("numeric vector validation rejects invalid inputs", {
  expect_invisible(gaussianprocesses:::.assert_numeric_vector(c(1, 2, 3)))

  expect_error(
    gaussianprocesses:::.assert_numeric_vector(character()),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    gaussianprocesses:::.assert_numeric_vector(numeric()),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    gaussianprocesses:::.assert_numeric_vector(c(1, NA_real_)),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("numeric matrix validation is explicit", {
  x <- matrix(1:6, nrow = 2)

  expect_invisible(gaussianprocesses:::.assert_numeric_matrix(x))
  expect_error(
    gaussianprocesses:::.assert_numeric_matrix(data.frame(x = 1:3)),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    gaussianprocesses:::.assert_square_matrix(x),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("symmetry uses a relative tolerance", {
  x <- matrix(c(2, 1, 1, 2), nrow = 2)
  perturbed <- x
  perturbed[1, 2] <- perturbed[1, 2] + 1e-12

  expect_invisible(
    gaussianprocesses:::.assert_symmetric_matrix(
      perturbed,
      tolerance = 1e-10
    )
  )

  expect_error(
    gaussianprocesses:::.assert_symmetric_matrix(
      perturbed,
      tolerance = 1e-14
    ),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("scalar validation checks bounds and integer counts", {
  expect_invisible(
    gaussianprocesses:::.assert_numeric_scalar(0, "value", lower = 0)
  )

  expect_error(
    gaussianprocesses:::.assert_numeric_scalar(
      0,
      "value",
      lower = 0,
      lower_open = TRUE
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_equal(
    gaussianprocesses:::.assert_positive_integer(3, "count"),
    3L
  )

  expect_error(
    gaussianprocesses:::.assert_positive_integer(1.5, "count"),
    class = "gaussianprocesses_validation_error"
  )
})
