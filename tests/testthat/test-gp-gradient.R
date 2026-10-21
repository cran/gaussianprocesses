# Posterior gradient of the latent function.

gradient_models <- function() {
  x_1d <- seq(-2, 2, length.out = 15)
  x_3d <- cbind(
    seq(-1, 1, length.out = 20),
    sin(1:20),
    cos(seq(0, 3, length.out = 20))
  )
  y_3d <- sin(x_3d[, 1]) + x_3d[, 2] * x_3d[, 3]

  list(
    rbf_1d = fit_gp(
      x_1d,
      sin(2 * x_1d),
      rbf_kernel(variance = 1.3, length_scale = 0.7),
      noise_variance = 0.01
    ),
    matern52_1d = fit_gp(
      x_1d,
      abs(x_1d),
      matern52_kernel(length_scale = 0.9),
      noise_variance = 0.02,
      mean = constant_mean(0.5)
    ),
    ard_3d = fit_gp(
      x_3d,
      y_3d,
      rbf_kernel(length_scale = c(0.8, 1.2, 1.6)),
      noise_variance = 0.01
    ),
    composite_3d = fit_gp(
      x_3d,
      y_3d,
      sum_kernel(
        product_kernel(
          matern52_kernel(length_scale = c(1, 1.5, 2)),
          linear_kernel(variance = 0.5)
        ),
        periodic_kernel(variance = 0.3, period = 2.5, length_scale = 1.1)
      ),
      noise_variance = 0.02,
      mean = constant_mean(-0.2)
    )
  )
}

gradient_prediction_inputs <- function(model) {
  if (model$n_features == 1L) {
    return(c(-1.7, -0.45, 0.1, 0.83, 2.6))
  }

  cbind(c(-0.7, 0.05, 0.6, 1.4), c(0.3, -0.9, 0.5, 0), c(1, 0.2, -0.4, 0.7))
}

shift_prediction_inputs <- function(x, dimension, step) {
  x <- as.matrix(x)
  x[, dimension] <- x[, dimension] + step
  x
}


test_that("the posterior gradient mean differentiates the posterior mean", {
  models <- gradient_models()
  step <- 1e-5

  for (name in names(models)) {
    model <- models[[name]]
    x <- gradient_prediction_inputs(model)
    gradient <- predict_gradient_gp(model, x)

    numerical <- vapply(
      seq_len(model$n_features),
      function(d) {
        (predict_gp(model, shift_prediction_inputs(x, d, step))$mean -
          predict_gp(model, shift_prediction_inputs(x, d, -step))$mean) /
          (2 * step)
      },
      numeric(NROW(x))
    )

    expect_lte(
      max(abs(unname(gradient$mean) - matrix(numerical, NROW(x)))) /
        max(abs(gradient$mean)),
      1e-6,
      label = name
    )
  }
})


test_that("the posterior gradient variance is a difference-quotient limit", {
  models <- gradient_models()

  # Variance of (f(x + h e_d) - f(x - h e_d)) / (2 h) under the posterior.
  difference_quotient_variance <- function(model, x, dimension, step) {
    vapply(
      seq_len(nrow(x)),
      function(i) {
        point <- x[i, , drop = FALSE]
        covariance <- predict_gp(
          model,
          rbind(
            shift_prediction_inputs(point, dimension, step),
            shift_prediction_inputs(point, dimension, -step)
          ),
          include_covariance = TRUE
        )$latent_covariance

        (covariance[1, 1] + covariance[2, 2] - 2 * covariance[1, 2]) /
          (4 * step^2)
      },
      numeric(1)
    )
  }

  for (name in names(models)) {
    model <- models[[name]]
    x <- as.matrix(gradient_prediction_inputs(model))
    gradient <- predict_gradient_gp(model, x)

    for (d in seq_len(model$n_features)) {
      quotients <- lapply(
        c(0.02, 0.01, 0.005),
        function(step) difference_quotient_variance(model, x, d, step)
      )
      # Richardson extrapolation removes the O(step^2) error term and then the
      # O(step^3) term, which Matern-5/2 kernels have.
      first <- Map(
        function(coarse, fine) (4 * fine - coarse) / 3,
        quotients[-3],
        quotients[-1]
      )
      limit <- (8 * first[[2]] - first[[1]]) / 7

      expect_lte(
        max(abs(gradient$variance[, d] - limit)) / max(gradient$variance[, d]),
        1e-5,
        label = sprintf("%s, dimension %d", name, d)
      )
    }
  }
})


test_that("gradient predictions do not depend on the block size", {
  model <- gradient_models()$composite_3d
  x <- rbind(gradient_prediction_inputs(model), c(0, 0, 0), c(2, -1, 0.5))
  reference <- predict_gradient_gp(model, x)

  for (block_size in c(1L, 2L, 5L, 100L)) {
    blocked <- predict_gradient_gp(model, x, block_size = block_size)

    for (field in c("mean", "variance", "sd", "lower", "upper")) {
      expect_equal(
        blocked[[field]],
        reference[[field]],
        tolerance = 1e-12,
        info = sprintf("%s, block size %d", field, block_size)
      )
    }
  }
})


test_that("gradient predictions have named columns and credible intervals", {
  model <- gradient_models()$ard_3d
  x <- gradient_prediction_inputs(model)
  gradient <- predict_gradient_gp(model, x, interval_level = 0.9)

  expect_s3_class(gradient, "gaussianprocesses_gradient_prediction")
  expect_identical(gradient$x, x)
  expect_identical(gradient$interval_level, 0.9)

  for (field in c("mean", "variance", "sd", "lower", "upper")) {
    expect_identical(dim(gradient[[field]]), c(4L, 3L), info = field)
    expect_identical(colnames(gradient[[field]]), c("x1", "x2", "x3"))
  }

  expect_equal(gradient$sd, sqrt(gradient$variance))
  expect_equal(
    gradient$upper - gradient$mean,
    stats::qnorm(0.95) * gradient$sd
  )
  expect_equal(gradient$mean - gradient$lower, gradient$upper - gradient$mean)
  expect_true(all(gradient$variance > 0))

  named <- x
  colnames(named) <- c("a", "b", "c")
  expect_identical(
    colnames(predict_gradient_gp(model, named)$mean),
    c("a", "b", "c")
  )

  expect_output(print(gradient), "posterior gradient")
  expect_output(print(gradient), "input dimensions: 3")
  expect_output(print(gradient), "interval level: 90%")
  utils::capture.output(expect_invisible(print(gradient)))
})


test_that("the gradient of a smooth function is recovered from dense data", {
  x <- seq(0, 2 * pi, length.out = 40)
  model <- fit_gp(
    x,
    sin(x),
    rbf_kernel(length_scale = 1.5),
    noise_variance = 1e-6
  )
  points <- seq(1, 5, length.out = 9)

  gradient <- predict_gradient_gp(model, points)

  expect_equal(gradient$mean[, 1], cos(points), tolerance = 1e-3)
  expect_true(
    all(abs(gradient$mean[, 1] - cos(points)) < 3 * gradient$sd[, 1])
  )

  # Far from the data the gradient reverts to its prior, centred on the mean
  # function's zero gradient.
  far <- predict_gradient_gp(model, 50)
  expect_equal(unname(far$mean[1, 1]), 0, tolerance = 1e-12)
  expect_equal(unname(far$variance[1, 1]), 1 / 1.5^2, tolerance = 1e-12)
})


test_that("gradient prediction validates its inputs", {
  model <- gradient_models()$rbf_1d

  expect_error(
    predict_gradient_gp(list(), 1),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    predict_gradient_gp(model, cbind(1, 2)),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    predict_gradient_gp(model, 1, interval_level = 1),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    predict_gradient_gp(model, 1, variance_tolerance = -1),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    predict_gradient_gp(model, 1, block_size = 0),
    class = "gaussianprocesses_validation_error"
  )

  sparse <- fit_sparse_gp(
    seq(0, 1, length.out = 12),
    sin(1:12),
    rbf_kernel(),
    noise_variance = 0.01,
    n_inducing = 4
  )
  expect_error(
    predict_gradient_gp(sparse, 0.5),
    class = "gaussianprocesses_validation_error"
  )

  for (kernel in list(
    matern12_kernel(),
    sum_kernel(rbf_kernel(), white_noise_kernel(variance = 0.01))
  )) {
    rough <- fit_gp(1:6, sin(1:6), kernel, noise_variance = 0.01)
    expect_error(
      predict_gradient_gp(rough, 2.5),
      class = "gaussianprocesses_smoothness_error"
    )
  }
})

