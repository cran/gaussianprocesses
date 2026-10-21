# VFE, sparse gradients, sparse optimization, and sparse sampling.

sparse_data <- function(n = 60, seed = 35) {
  set.seed(seed)
  x <- cbind(stats::runif(n, 0, 6), stats::runif(n, -1, 1))

  list(x = x, y = sin(x[, 1]) + 0.5 * x[, 2] + stats::rnorm(n, sd = 0.2))
}


test_that("the VFE bound lies below the exact log marginal likelihood", {
  set.seed(1)
  x <- sort(stats::runif(40, 0, 6))
  y <- sin(x) + stats::rnorm(40, sd = 0.2)
  gaps <- vapply(seq_len(50), function(case) {
    kernel <- rbf_kernel(exp(stats::runif(1, -2, 2)),
                         exp(stats::runif(1, -1.5, 1)))
    noise <- if (case %% 5 == 0) {
      exp(stats::runif(40, -5, 0))
    } else {
      exp(stats::runif(1, -5, 0))
    }
    model <- fit_sparse_gp(x, y, kernel, noise_variance = noise,
                           n_inducing = sample(3:15, 1), method = "vfe")
    model$log_marginal_likelihood -
      log_marginal_likelihood(fit_gp(x, y, kernel, noise_variance = noise))
  }, numeric(1))

  expect_lt(max(gaps), 1e-10)
})


test_that("VFE equals the exact model when every input is inducing", {
  x <- seq(0, 6, length.out = 25)
  y <- sin(x) + 0.1 * cos(5 * x)
  kernel <- matern52_kernel(1.2, 0.8)

  vfe <- fit_sparse_gp(x, y, kernel, noise_variance = 0.04,
                       inducing_points = x, method = "vfe")
  exact <- fit_gp(x, y, kernel, noise_variance = 0.04)
  expect_identical(vfe$kuu_jitter, 0)
  expect_lt(abs(vfe$log_marginal_likelihood - log_marginal_likelihood(exact)),
            1e-8)
  expect_lt(vfe$trace_term, 1e-8)

  x_new <- c(-0.5, 2.2, 7)
  expect_equal(predict_sparse_gp(vfe, x_new)$mean,
               predict_gp(exact, x_new)$mean, tolerance = 1e-8)

  # With jitter in K_uu, the bound uses Q = K_fu (K_uu + delta I)^{-1} K_uf
  # throughout: it equals the dense formula with that Q and stays below the
  # exact likelihood.
  dense_x <- seq(0, 1, length.out = 30)
  dense_y <- sin(3 * dense_x)
  smooth <- rbf_kernel(1, 2)
  jittered <- fit_sparse_gp(dense_x, dense_y, smooth, noise_variance = 0.01,
                            inducing_points = dense_x, method = "vfe")
  expect_gt(jittered$kuu_jitter, 0)

  kuu <- evaluate_kernel(smooth, dense_x) + diag(jittered$kuu_jitter, 30)
  kuf <- evaluate_kernel(smooth, dense_x)
  low_rank <- crossprod(kuf, solve(kuu, kuf))
  covariance <- low_rank + diag(0.01, 30)
  expected <- -0.5 * drop(crossprod(dense_y, solve(covariance, dense_y))) -
    0.5 * as.numeric(determinant(covariance)$modulus) -
    15 * log(2 * pi) -
    0.5 * sum(diag(evaluate_kernel(smooth, dense_x)) - diag(low_rank)) / 0.01
  expect_equal(jittered$log_marginal_likelihood, expected, tolerance = 1e-8)
  expect_lte(
    jittered$log_marginal_likelihood,
    log_marginal_likelihood(fit_gp(dense_x, dense_y, smooth,
                                   noise_variance = 0.01))
  )
})


test_that("VFE matches its dense bound and predictive equations", {
  data <- sparse_data()
  kernel <- rbf_kernel(1.2, c(0.8, 1.5))
  model <- fit_sparse_gp(data$x, data$y, kernel, noise_variance = 0.04,
                         n_inducing = 9, method = "vfe")
  expect_identical(model$method, "vfe")
  expect_identical(model$approximation, "VFE")

  kuu <- model$kuu
  kuf <- model$kuf
  low_rank <- crossprod(kuf, solve(kuu, kuf))
  residual <- diag(evaluate_kernel(kernel, data$x)) - diag(low_rank)
  covariance <- low_rank + diag(0.04, 60)
  expected <- -0.5 * drop(crossprod(data$y, solve(covariance, data$y))) -
    0.5 * as.numeric(determinant(covariance)$modulus) -
    30 * log(2 * pi) - 0.5 * sum(residual) / 0.04
  expect_equal(model$log_marginal_likelihood, expected, tolerance = 1e-8)
  expect_equal(model$trace_term, sum(residual) / 0.04, tolerance = 1e-8)

  # Mean sigma^{-2} K_*u Sigma K_uf y and covariance
  # K_** - Q_** + K_*u Sigma K_u* with Sigma = (K_uu + sigma^{-2} K_uf K_fu)^{-1}.
  x_new <- cbind(c(0.5, 3, 5.5), c(-0.5, 0, 0.8))
  ksu <- evaluate_kernel(kernel, x_new, model$inducing_points)
  sigma <- solve(kuu + tcrossprod(kuf) / 0.04)
  expected_mean <- drop(ksu %*% sigma %*% kuf %*% data$y) / 0.04
  expected_covariance <- evaluate_kernel(kernel, x_new) -
    ksu %*% solve(kuu, t(ksu)) + ksu %*% sigma %*% t(ksu)
  prediction <- predict_sparse_gp(model, x_new, include_covariance = TRUE)
  expect_equal(prediction$mean, expected_mean, tolerance = 1e-8)
  expect_equal(prediction$latent_covariance, expected_covariance,
               tolerance = 1e-8)
  expect_identical(prediction$approximation, "VFE")

  # Marginal predictions in any block size agree with the full covariance.
  for (block_size in list(NULL, 1L, 2L)) {
    marginal <- predict_sparse_gp(model, x_new, block_size = block_size)
    expect_equal(marginal$mean, prediction$mean, tolerance = 1e-12)
    expect_equal(marginal$latent_variance, prediction$latent_variance,
                 tolerance = 1e-12)
  }

  expect_error(
    fit_sparse_gp(data$x, data$y, kernel, noise_variance = 0, method = "vfe"),
    class = "gaussianprocesses_validation_error"
  )
  expect_true(is.na(fit_sparse_gp(data$x, data$y, kernel, noise_variance = 0,
                                  n_inducing = 9)$trace_term))
})


test_that("FITC and VFE gradients pass finite-difference checks", {
  data <- sparse_data()
  inducing <- data$x[seq(1, 60, by = 6), ]
  cases <- list(
    list(kernel = rbf_kernel(1.2, c(0.8, 1.5)), noise = 0.04),
    list(kernel = sum_kernel(matern52_kernel(1, 0.9), linear_kernel(0.3)),
         noise = 0.1),
    list(kernel = rbf_kernel(0.7, 1.1), noise = seq(0.02, 0.2, length.out = 60))
  )

  for (method in c("fitc", "vfe")) {
    for (case in cases) {
      pattern <- length(case$noise) > 1L
      n_kernel <- length(kernel_parameters(case$kernel, flatten = TRUE))
      fit_at <- function(theta) {
        fit_sparse_gp(
          data$x,
          data$y,
          update_kernel_parameters(case$kernel, exp(theta[seq_len(n_kernel)])),
          noise_variance = exp(theta[[n_kernel + 1L]]) *
            if (pattern) case$noise else 1,
          inducing_points = inducing,
          method = method
        )
      }
      theta <- c(log(kernel_parameters(case$kernel, flatten = TRUE)),
                 if (pattern) 0 else log(case$noise))
      analytical <- .sparse_gradient(
        fit_at(theta),
        include_noise = TRUE,
        noise_parameter = if (pattern) "noise_scale" else "noise_variance"
      )
      numerical <- vapply(seq_along(theta), function(j) {
        shifted <- function(step) {
          up <- theta
          down <- theta
          up[j] <- up[j] + step
          down[j] <- down[j] - step
          (fit_at(up)$log_marginal_likelihood -
              fit_at(down)$log_marginal_likelihood) / (2 * step)
        }
        (4 * shifted(5e-6) - shifted(1e-5)) / 3
      }, numeric(1))

      expect_lt(
        max(abs(analytical - numerical) / pmax(abs(numerical), 1e-3)),
        1e-6,
        label = paste(method, utils::capture.output(print(case$kernel))[1])
      )
    }
  }
})


test_that("sparse hyperparameters are optimized with the shared optimizer", {
  set.seed(3)
  x <- sort(stats::runif(30, 0, 6))
  y <- sin(x) + stats::rnorm(30, sd = 0.2)

  # With every input inducing, VFE and FITC reach the exact optimum.
  exact <- optimize_gp(x, y, rbf_kernel(), noise_variance = 0.1, n_starts = 1)

  for (method in c("vfe", "fitc")) {
    sparse <- optimize_sparse_gp(x, y, rbf_kernel(), noise_variance = 0.1,
                                 inducing_points = x, method = method,
                                 n_starts = 1)
    expect_equal(
      unname(sparse$optimization$optimized_parameters),
      unname(exact$optimization$optimized_parameters),
      tolerance = 1e-4,
      label = method
    )
    expect_identical(sparse$optimization$approximation, toupper(method))
  }

  model <- optimize_sparse_gp(x, y, rbf_kernel(), noise_variance = 0.1,
                              n_inducing = 8, n_starts = 2)
  record <- model$optimization
  expect_identical(model$method, "vfe")
  expect_true(record$converged)
  expect_named(record$optimized_parameters,
               c("variance", "length_scale", "noise_variance"))
  expect_identical(record$log_marginal_likelihood,
                   model$log_marginal_likelihood)
  expect_identical(attr(logLik(model), "df"), 3L)
  expect_gt(record$model_fits, 0L)
  expect_output(print(model), "optimization: converged")

  numerical <- optimize_sparse_gp(x, y, rbf_kernel(), noise_variance = 0.1,
                                  n_inducing = 8, n_starts = 2,
                                  gradient = "numerical")
  expect_equal(numerical$log_marginal_likelihood,
               model$log_marginal_likelihood, tolerance = 1e-6)

  # Fixed noise, and the scale of a noise pattern.
  fixed <- optimize_sparse_gp(x, y, rbf_kernel(), noise_variance = 0.05,
                              n_inducing = 8, optimize_noise = FALSE,
                              n_starts = 1)
  expect_identical(fixed$noise_variance, 0.05)
  expect_named(fixed$optimization$optimized_parameters,
               c("variance", "length_scale"))
  pattern <- optimize_sparse_gp(x, y, rbf_kernel(),
                                noise_variance = seq(0.5, 1.5, length.out = 30),
                                n_inducing = 8, n_starts = 1)
  expect_identical(pattern$optimization$noise_parameterization, "scale")

  for (call in list(
    function() optimize_sparse_gp(x, y, rbf_kernel(), noise_variance = 0,
                                  optimize_noise = FALSE, n_inducing = 8),
    function() optimize_sparse_gp(x, y, rbf_kernel(), n_inducing = 8,
                                  optimize_noise = NA),
    function() optimize_sparse_gp(x, y, rbf_kernel(), n_inducing = 8,
                                  mean = linear_mean())
  )) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }
})


test_that("sparse posterior draws pass moment tests", {
  data <- sparse_data()
  x_new <- cbind(c(0.5, 0.7, 3, 5.5), c(-0.5, -0.4, 0, 0.8))

  for (method in c("fitc", "vfe")) {
    model <- fit_sparse_gp(data$x, data$y, rbf_kernel(1.2, c(0.8, 1.5)),
                           noise_variance = 0.04, n_inducing = 9,
                           method = method)
    prediction <- predict_sparse_gp(model, x_new, include_covariance = TRUE)
    draws <- sample_gp_posterior(model, x_new, n_draws = n_statistical_draws,
                                 type = "observation", seed = 9)

    expect_identical(draws$model, model$approximation)
    expect_equal(draws$latent_covariance, prediction$latent_covariance)
    expect_gaussian_moments(draws$latent, prediction$mean,
                            prediction$latent_covariance,
                            paste(method, "latent"))
    expect_gaussian_moments(draws$observation, prediction$mean,
                            prediction$observation_covariance,
                            paste(method, "observation"))
  }

  heteroscedastic <- fit_sparse_gp(data$x, data$y, rbf_kernel(),
                                   noise_variance = rep(0.05, 60),
                                   n_inducing = 9)
  expect_null(sample_gp_posterior(heteroscedastic, x_new, seed = 1)$observation)
  expect_error(
    sample_gp_posterior(heteroscedastic, x_new, type = "observation"),
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    sample_gp_posterior(heteroscedastic, x_new,
                        observation_noise_variance = 0.1),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("the trace term measures inducing-point coverage", {
  set.seed(4)
  x <- sort(stats::runif(200, 0, 10))
  y <- sin(x) + stats::rnorm(200, sd = 0.1)
  trace <- vapply(c(5, 10, 20, 40), function(m) {
    fit_sparse_gp(x, y, rbf_kernel(1, 1), noise_variance = 0.01,
                  n_inducing = m, method = "vfe")$trace_term
  }, numeric(1))

  expect_true(all(diff(trace) < 0))
  expect_lt(trace[4], 1e-3)
  expect_gt(trace[1], 100)

  fitc <- fit_sparse_gp(x, y, rbf_kernel(1, 1), noise_variance = 0.01,
                        n_inducing = 10)
  expect_equal(fitc$trace_term, sum(fitc$residual_diagonal) / 0.01)
  expect_output(print(fitc), "trace term")
  expect_match(
    utils::capture.output(print(summary(fitc)))[1],
    "Sparse Gaussian-process regression \\(FITC\\) model"
  )
})


test_that("schema-3 sparse models upgrade to FITC with their trace term", {
  data <- sparse_data()
  model <- fit_sparse_gp(data$x, data$y, rbf_kernel(), noise_variance = 0.04,
                         n_inducing = 9)
  old <- model
  old$method <- NULL
  old$trace_term <- NULL
  old$numerical_policy <- NULL
  old$schema_version <- 3L

  upgraded <- .prepare_sparse_gp_model(old)
  expect_identical(upgraded$method, "fitc")
  expect_equal(upgraded$trace_term, model$trace_term)
  expect_identical(upgraded$schema_version, .model_schema_version)
  expect_identical(predict_sparse_gp(old, data$x[1:3, ]),
                   predict_sparse_gp(model, data$x[1:3, ]))
})
