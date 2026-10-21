test_that("kernel diagonal avoids full covariance construction semantics", {
  x <- matrix(
    c(
      1, 2,
      3, 4,
      5, 6
    ),
    ncol = 2,
    byrow = TRUE
  )

  kernels <- list(
    rbf_kernel(variance = 2, length_scale = c(1, 3)),
    linear_kernel(variance = 0.5),
    sum_kernel(
      matern32_kernel(variance = 1.2, length_scale = 2),
      white_noise_kernel(variance = 0.1)
    ),
    product_kernel(
      rbf_kernel(variance = 1.5, length_scale = 1),
      linear_kernel(variance = 0.2)
    )
  )

  for (kernel in kernels) {
    expect_equal(
      kernel_diagonal(kernel, x),
      diag(evaluate_kernel(kernel, x)),
      tolerance = 1e-12
    )
  }
})


test_that("farthest inducing-point selection is deterministic and unique", {
  x <- matrix(
    c(
      0, 0,
      1, 0,
      0, 1,
      1, 1,
      0.5, 0.5
    ),
    ncol = 2,
    byrow = TRUE
  )

  first <- select_inducing_points(
    x,
    n_inducing = 3,
    method = "farthest"
  )
  second <- select_inducing_points(
    x,
    n_inducing = 3,
    method = "farthest"
  )

  expect_equal(first, second)
  expect_equal(
    attr(first, "indices"),
    attr(second, "indices")
  )
  expect_equal(nrow(unique(first)), 3L)
})


test_that("quantile inducing selection spans one-dimensional inputs", {
  x <- seq(0, 10, by = 1)

  inducing <- select_inducing_points(
    x,
    n_inducing = 3,
    method = "quantile"
  )

  expect_equal(nrow(inducing), 3L)
  expect_equal(inducing[1, 1], 0)
  expect_equal(inducing[3, 1], 10)

  expect_error(
    select_inducing_points(
      cbind(x, x),
      n_inducing = 3,
      method = "quantile"
    ),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("FITC recovers exact GP when every training point is inducing", {
  x <- c(-1, -0.2, 0.5, 1.4)
  y <- c(0.4, -0.1, 0.8, 0.2)
  kernel <- rbf_kernel(
    variance = 1.3,
    length_scale = 0.7
  )
  noise <- 0.12

  exact <- fit_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = noise
  )
  sparse <- fit_sparse_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = noise,
    inducing_points = x
  )

  x_new <- c(-0.7, 0.1, 1.0)

  exact_prediction <- predict_gp(exact, x_new)
  sparse_prediction <- predict_sparse_gp(
    sparse,
    x_new
  )

  expect_equal(
    sparse_prediction$mean,
    exact_prediction$mean,
    tolerance = 1e-7
  )
  expect_equal(
    sparse_prediction$latent_variance,
    exact_prediction$latent_variance,
    tolerance = 1e-7
  )
})


test_that("sparse model retains low-rank rather than full training covariance", {
  x <- seq(-2, 2, length.out = 30)
  y <- sin(x)

  model <- fit_sparse_gp(
    x,
    y,
    kernel = matern32_kernel(length_scale = 0.8),
    noise_variance = 0.05,
    n_inducing = 6
  )

  expect_equal(dim(model$kuf), c(6L, 30L))
  expect_equal(dim(model$kuu), c(6L, 6L))
  expect_false("observation_covariance" %in% names(model))
  expect_equal(model$approximation, "FITC")
})


test_that("sparse prediction marginal and full covariance agree", {
  x <- seq(0, 1, length.out = 15)
  y <- sin(2 * pi * x)

  model <- fit_sparse_gp(
    x,
    y,
    kernel = rbf_kernel(length_scale = 0.2),
    noise_variance = 0.03,
    n_inducing = 5,
    selection = "quantile"
  )

  marginal <- predict_sparse_gp(
    model,
    c(0.2, 0.4, 0.8)
  )
  full <- predict_sparse_gp(
    model,
    c(0.2, 0.4, 0.8),
    include_covariance = TRUE
  )

  expect_equal(
    diag(full$latent_covariance),
    marginal$latent_variance,
    tolerance = 1e-12
  )
  expect_equal(
    diag(full$observation_covariance),
    marginal$observation_variance,
    tolerance = 1e-12
  )
})


test_that("FITC supports known heteroscedastic training noise", {
  x <- seq(0, 1, length.out = 12)
  y <- cos(2 * pi * x)
  noise <- 0.01 + 0.05 * x

  model <- fit_sparse_gp(
    x,
    y,
    kernel = rbf_kernel(length_scale = 0.25),
    noise_variance = noise,
    n_inducing = 4
  )

  expect_equal(
    model$training_noise_variance,
    noise
  )

  expect_error(
    predict_sparse_gp(model, c(0.2, 0.8)),
    class = "gaussianprocesses_validation_error"
  )

  prediction <- predict_sparse_gp(
    model,
    c(0.2, 0.8),
    observation_noise_variance = c(0.02, 0.05)
  )

  expect_equal(
    prediction$observation_variance -
      prediction$latent_variance,
    c(0.02, 0.05),
    tolerance = 1e-12
  )
})


test_that("benchmark reports runtime memory and predictive differences", {
  x <- seq(-1, 1, length.out = 20)
  y <- sin(2 * x)
  kernel <- matern32_kernel(
    variance = 1,
    length_scale = 0.6
  )

  benchmark <- benchmark_sparse_gp(
    x,
    y,
    kernel = kernel,
    n_inducing = 5,
    prediction_x = seq(-1, 1, length.out = 10),
    noise_variance = 0.05
  )

  comparison <- benchmark$comparison

  expect_equal(comparison$model, c("exact", "FITC"))
  expect_true(all(comparison$retained_model_bytes > 0))
  expect_gte(
    comparison$posterior_mean_rmse_vs_exact[2],
    0
  )
  expect_gte(
    comparison$latent_variance_mae_vs_exact[2],
    0
  )
})


test_that("sparse model validates inducing and floor settings", {
  x <- matrix(1:12, ncol = 2)
  y <- seq_len(nrow(x))

  expect_error(
    fit_sparse_gp(
      x,
      y,
      kernel = rbf_kernel(),
      n_inducing = nrow(x) + 1
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    fit_sparse_gp(
      x,
      y,
      kernel = rbf_kernel(),
      inducing_points = matrix(1:3, ncol = 1)
    ),
    class = "gaussianprocesses_validation_error"
  )

  expect_error(
    fit_sparse_gp(
      x,
      y,
      kernel = rbf_kernel(),
      n_inducing = 3,
      fitc_diagonal_floor = 0
    ),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("automatic inducing selection never returns duplicated inputs", {
  x <- rep(c(0, 1, 2, 3), each = 3)

  farthest <- select_inducing_points(
    x,
    n_inducing = 4,
    method = "farthest"
  )
  expect_equal(sort(farthest[, 1]), c(0, 1, 2, 3))
  expect_equal(x[attr(farthest, "indices")], farthest[, 1])

  # Quantile targets 0, 0, 2 previously selected two copies of 0.
  quantile <- select_inducing_points(
    c(0, 0, 0, 0, 1, 2),
    n_inducing = 3,
    method = "quantile"
  )
  expect_identical(anyDuplicated(quantile[, 1]), 0L)
  expect_equal(sort(quantile[, 1]), c(0, 1, 2))

  expect_error(
    select_inducing_points(x, n_inducing = 5),
    class = "gaussianprocesses_validation_error"
  )

  repeated_rows <- rbind(c(0, 0), c(0, 0), c(1, 1), c(1, 1))
  expect_equal(
    nrow(unique(select_inducing_points(repeated_rows, 2))),
    2L
  )
  expect_error(
    select_inducing_points(repeated_rows, 3),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("default inducing count respects distinct training inputs", {
  x <- rep(seq(0, 5), each = 4)

  model <- fit_sparse_gp(
    x,
    sin(x),
    kernel = rbf_kernel(),
    noise_variance = 0.01
  )

  expect_equal(model$n_inducing, 6L)
  expect_equal(nrow(unique(model$inducing_points)), 6L)
  expect_equal(model$kuu_jitter, 0)
})


test_that("FITC log marginal likelihood matches the dense FITC density", {
  x <- seq(-2, 2, length.out = 25)
  y <- sin(1.5 * x) + 0.1 * cos(5 * x)
  kernel <- matern52_kernel(variance = 1.1, length_scale = 0.9)
  noise <- 0.05

  model <- fit_sparse_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = noise,
    n_inducing = 7
  )
  expect_equal(model$kuu_jitter, 0)

  kuu <- model$kuu
  kuf <- model$kuf
  low_rank <- crossprod(kuf, solve(kuu, kuf))
  lambda <- diag(evaluate_kernel(kernel, x)) - diag(low_rank) + noise
  covariance <- low_rank + diag(lambda)

  expected <- -0.5 * drop(crossprod(y, solve(covariance, y))) -
    0.5 * as.numeric(determinant(covariance)$modulus) -
    0.5 * length(y) * log(2 * pi)

  expect_equal(model$log_marginal_likelihood, expected, tolerance = 1e-8)
  expect_equal(
    model$training_alpha,
    drop(solve(covariance, y)),
    tolerance = 1e-8
  )
  expect_equal(
    crossprod(model$a_cholesky),
    kuu + kuf %*% (t(kuf) / lambda),
    tolerance = 1e-8
  )
})


test_that("FITC is invariant to an aliased inducing point", {
  # Under a periodic kernel, u and u + period have the same covariance with
  # every input. The alias adds no information but makes K_uu singular, so
  # K_uu needs jitter. That jitter must be used consistently: previously A
  # needed its own jitter and the approximate log marginal likelihood moved.
  x <- seq(0.02, 0.98, length.out = 25)
  y <- cos(2 * pi * x) + 0.3 * sin(6 * pi * x)
  kernel <- periodic_kernel(
    variance = 1,
    length_scale = 0.6,
    period = 1
  )
  x_new <- c(0.33, 0.81)

  base <- fit_sparse_gp(
    x,
    y,
    kernel = kernel,
    noise_variance = 1e-6,
    inducing_points = c(0.1, 0.4, 0.7)
  )
  # 1.1 aliases 0.1 under a kernel of period 1, which the collapse check
  # reports.
  expect_warning(
    aliased <- fit_sparse_gp(
      x,
      y,
      kernel = kernel,
      noise_variance = 1e-6,
      inducing_points = c(0.1, 0.4, 0.7, 1.1)
    ),
    class = "gaussianprocesses_inducing_collapse_warning"
  )
  expect_identical(c(aliased$inducing_collapse$first,
                     aliased$inducing_collapse$second), c(1L, 4L))

  expect_equal(base$kuu_jitter, 0)
  expect_gt(aliased$kuu_jitter, 0)
  expect_equal(aliased$a_jitter, 0)

  expect_equal(
    aliased$log_marginal_likelihood,
    base$log_marginal_likelihood,
    tolerance = 1e-6
  )

  base_prediction <- predict_sparse_gp(base, x_new)
  aliased_prediction <- predict_sparse_gp(aliased, x_new)

  expect_equal(
    aliased_prediction$mean,
    base_prediction$mean,
    tolerance = 1e-6
  )
  expect_equal(
    aliased_prediction$latent_variance,
    base_prediction$latent_variance,
    tolerance = 1e-6
  )
})


test_that("sparse objects have concise print methods", {
  x <- seq(0, 1, length.out = 12)
  y <- sin(2 * pi * x)

  model <- fit_sparse_gp(
    x,
    y,
    kernel = rbf_kernel(length_scale = 0.2),
    noise_variance = 0.05,
    n_inducing = 4
  )
  prediction <- predict_sparse_gp(
    model,
    c(0.25, 0.75)
  )

  model_output <- capture.output(print(model))
  prediction_output <- capture.output(print(prediction))

  expect_true(
    any(grepl("FITC", model_output, fixed = TRUE))
  )
  expect_true(
    any(grepl("inducing points: 4", model_output, fixed = TRUE))
  )
  expect_true(
    any(grepl("Sparse Gaussian-process prediction", prediction_output, fixed = TRUE))
  )
})
