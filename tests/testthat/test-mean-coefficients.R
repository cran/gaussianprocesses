# Parametric means with estimated or marginalized coefficients.

coefficient_data <- function(n = 30, seed = 3) {
  set.seed(seed)
  x <- cbind(runif(n), runif(n))
  y <- 1 + 2 * x[, 1] - x[, 2] + sin(5 * x[, 1]) + rnorm(n, sd = 0.05)

  list(
    x = x,
    y = y,
    kernel = rbf_kernel(variance = 1.5, length_scale = c(0.3, 0.45)),
    noise_variance = 0.02,
    x_new = cbind(c(0.2, 0.9, 1.5, 0.5), c(0.5, 0.1, -0.5, 0.5))
  )
}

coefficient_fit <- function(data, mean, ...) {
  fit_gp(
    data$x,
    data$y,
    data$kernel,
    noise_variance = data$noise_variance,
    mean = mean,
    ...
  )
}

# Dense reference quantities, with solve() on small matrices.
dense_reference <- function(data) {
  covariance <- evaluate_kernel(data$kernel, data$x) +
    diag(data$noise_variance, nrow(data$x))
  basis <- cbind(1, data$x)
  precision_basis <- solve(covariance, basis)
  gls <- solve(
    crossprod(basis, precision_basis),
    crossprod(precision_basis, data$y)
  )

  list(
    covariance = covariance,
    basis = basis,
    coefficients = drop(gls),
    coefficient_covariance = solve(crossprod(basis, precision_basis))
  )
}

gaussian_log_density <- function(y, mean, covariance) {
  residual <- y - mean
  -0.5 * drop(crossprod(residual, solve(covariance, residual))) -
    0.5 * as.numeric(determinant(covariance)$modulus) -
    0.5 * length(y) * log(2 * pi)
}


test_that("estimated coefficients are the generalized least-squares estimate", {
  data <- coefficient_data()
  reference <- dense_reference(data)

  for (method in c("ml", "reml")) {
    model <- coefficient_fit(data, linear_mean(estimate_coefficients(method)))

    expect_equal(
      unname(coef(model)),
      reference$coefficients,
      tolerance = 1e-10
    )
    expect_identical(names(coef(model)), c("intercept", "x1", "x2"))
    expect_equal(
      unname(vcov(model)),
      reference$coefficient_covariance,
      tolerance = 1e-10
    )
    expect_equal(
      model$alpha,
      drop(solve(
        reference$covariance,
        data$y - reference$basis %*% reference$coefficients
      )),
      tolerance = 1e-10
    )
  }

  profile <- coefficient_fit(data, linear_mean())
  expect_equal(
    log_marginal_likelihood(profile),
    gaussian_log_density(
      data$y,
      drop(reference$basis %*% reference$coefficients),
      reference$covariance
    ),
    tolerance = 1e-10
  )

  # The restricted likelihood is the likelihood of the error contrasts A'y,
  # for an orthonormal basis A of the complement of the column space of H,
  # up to the constant (1/2) log |H'H|.
  restricted <- coefficient_fit(
    data,
    linear_mean(estimate_coefficients("reml"))
  )
  basis <- reference$basis
  contrasts <- qr.Q(qr(basis), complete = TRUE)[, -seq_len(ncol(basis))]
  contrast_likelihood <- gaussian_log_density(
    drop(crossprod(contrasts, data$y)),
    0,
    crossprod(contrasts, reference$covariance %*% contrasts)
  )
  expect_equal(
    log_marginal_likelihood(restricted) +
      0.5 * as.numeric(determinant(crossprod(basis))$modulus),
    contrast_likelihood,
    tolerance = 1e-10
  )
})


test_that("a Gaussian coefficient prior adds h(x)' B h(x') to the kernel", {
  data <- coefficient_data()
  basis <- cbind(1, data$x)
  new_basis <- cbind(1, data$x_new)

  priors <- list(
    diagonal = list(mean = c(1, 0, 0), covariance = c(4, 2, 0.5)),
    matrix = list(
      mean = 0.5,
      covariance = matrix(c(2, 0.3, 0, 0.3, 1, -0.2, 0, -0.2, 1.5), 3)
    )
  )

  for (name in names(priors)) {
    prior <- priors[[name]]
    model <- coefficient_fit(
      data,
      linear_mean(coefficient_prior(prior$mean, prior$covariance))
    )
    b <- rep_len(prior$mean, 3)
    B <- if (is.matrix(prior$covariance)) {
      prior$covariance
    } else {
      diag(prior$covariance)
    }

    covariance <- evaluate_kernel(data$kernel, data$x) +
      diag(data$noise_variance, nrow(data$x)) +
      basis %*% B %*% t(basis)
    cross <- evaluate_kernel(data$kernel, data$x_new, data$x) +
      new_basis %*% B %*% t(basis)
    prior_new <- evaluate_kernel(data$kernel, data$x_new) +
      new_basis %*% B %*% t(new_basis)
    residual <- data$y - drop(basis %*% b)

    prediction <- predict_gp(model, data$x_new, include_covariance = TRUE)
    expect_equal(
      prediction$mean,
      drop(new_basis %*% b + cross %*% solve(covariance, residual)),
      tolerance = 1e-10,
      info = name
    )
    expect_equal(
      prediction$latent_covariance,
      prior_new - cross %*% solve(covariance, t(cross)),
      tolerance = 1e-10,
      info = name
    )
    expect_equal(
      predict_gp(model, data$x_new, block_size = 1)$latent_variance,
      prediction$latent_variance,
      tolerance = 1e-12,
      info = name
    )
    expect_equal(
      log_marginal_likelihood(model),
      gaussian_log_density(data$y, drop(basis %*% b), covariance),
      tolerance = 1e-10,
      info = name
    )
  }
})


test_that("in the vague limit the mean is the plug-in mean of the estimate", {
  data <- coefficient_data()
  estimated <- predict_gp(coefficient_fit(data, linear_mean()), data$x_new)
  vague_model <- coefficient_fit(data, linear_mean(coefficient_prior()))
  vague <- predict_gp(vague_model, data$x_new)
  wide <- predict_gp(
    coefficient_fit(
      data,
      linear_mean(coefficient_prior(covariance = 1e8))
    ),
    data$x_new
  )

  expect_equal(vague$mean, estimated$mean, tolerance = 1e-10)
  expect_equal(
    coef(vague_model),
    coef(coefficient_fit(data, linear_mean())),
    tolerance = 1e-10
  )
  expect_equal(wide$mean, vague$mean, tolerance = 1e-6)
  expect_equal(wide$latent_variance, vague$latent_variance, tolerance = 1e-5)

  # Plugging in the estimate understates the uncertainty, most of all away
  # from the data.
  expect_true(all(vague$latent_variance >= estimated$latent_variance))
  expect_gt(vague$latent_variance[3], 2 * estimated$latent_variance[3])

  # The vague-prior likelihood is the restricted likelihood.
  expect_equal(
    log_marginal_likelihood(vague_model),
    log_marginal_likelihood(
      coefficient_fit(data, linear_mean(estimate_coefficients("reml")))
    ),
    tolerance = 1e-12
  )
})


test_that("likelihood gradients match finite differences for every treatment", {
  data <- coefficient_data(n = 20)
  means <- list(
    profile = linear_mean(),
    restricted = polynomial_mean(2, estimate_coefficients("reml")),
    vague = constant_mean(coefficient_prior()),
    gaussian = linear_mean(coefficient_prior(0, c(3, 1, 1)))
  )
  step <- 1e-5

  for (name in names(means)) {
    model <- coefficient_fit(data, means[[name]])
    analytical <- log_marginal_likelihood_gradient(model)
    parameters <- c(
      kernel_parameters(data$kernel, flatten = TRUE),
      noise_variance = data$noise_variance
    )

    numerical <- vapply(
      names(parameters),
      function(path) {
        likelihood_at <- function(direction) {
          shifted <- parameters
          shifted[[path]] <- shifted[[path]] * exp(direction * step)
          kernel <- update_kernel_parameters(
            data$kernel,
            shifted[names(shifted) != "noise_variance"]
          )
          log_marginal_likelihood(
            fit_gp(
              data$x,
              data$y,
              kernel,
              noise_variance = shifted[["noise_variance"]],
              mean = means[[name]]
            )
          )
        }

        (likelihood_at(1) - likelihood_at(-1)) / (2 * step)
      },
      numeric(1)
    )

    expect_identical(names(analytical), names(numerical), info = name)
    expect_lte(
      max(abs(analytical - numerical)) / max(1, max(abs(analytical))),
      1e-6,
      label = name
    )
  }
})


test_that("optimize_gp() maximizes profile and restricted likelihoods", {
  data <- coefficient_data(n = 25)
  fits <- list()

  for (method in c("ml", "reml")) {
    mean <- linear_mean(estimate_coefficients(method))
    fit <- function() {
      optimize_gp(
        data$x,
        data$y,
        rbf_kernel(length_scale = c(0.5, 0.5)),
        noise_variance = 0.05,
        mean = mean,
        n_starts = 3
      )
    }
    model <- fit()

    expect_identical(fit()$optimization, model$optimization, info = method)
    expect_true(model$optimization$converged, info = method)
    expect_identical(
      model$optimization$likelihood,
      if (method == "ml") "profile" else "restricted"
    )
    expect_identical(model$optimization$mean_coefficients, coef(model))
    expect_lt(
      max(abs(log_marginal_likelihood_gradient(model))),
      1e-3,
      label = method
    )
    expect_gt(
      log_marginal_likelihood(model),
      log_marginal_likelihood(
        fit_gp(
          data$x,
          data$y,
          rbf_kernel(length_scale = c(0.5, 0.5)),
          noise_variance = 0.05,
          mean = mean
        )
      )
    )
    fits[[method]] <- model
  }

  # A vague prior has the restricted likelihood, so the same optimum.
  vague <- optimize_gp(
    data$x,
    data$y,
    rbf_kernel(length_scale = c(0.5, 0.5)),
    noise_variance = 0.05,
    mean = linear_mean(coefficient_prior()),
    n_starts = 3
  )
  expect_equal(
    vague$optimization$optimized_parameters,
    fits$reml$optimization$optimized_parameters,
    tolerance = 1e-4
  )
  expect_identical(vague$optimization$likelihood, "restricted")
})


test_that("rank-deficient bases raise classed errors", {
  data <- coefficient_data(n = 12)
  collinear <- function(x) cbind(1, x[, 1], 2 * x[, 1] - 1)

  cases <- list(
    collinear = list(
      x = data$x,
      mean = basis_mean(collinear)
    ),
    constant_column = list(
      x = cbind(data$x[, 1], 1),
      mean = linear_mean()
    ),
    too_few_values = list(
      x = rep(c(0, 1), 6),
      mean = polynomial_mean(2)
    ),
    more_coefficients_than_observations = list(
      x = data$x[1:2, ],
      mean = linear_mean()
    ),
    vague_prior = list(
      x = data$x,
      mean = basis_mean(collinear, coefficient_prior())
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    expect_error(
      fit_gp(
        case$x,
        data$y[seq_len(NROW(case$x))],
        rbf_kernel(),
        noise_variance = 0.01,
        mean = case$mean
      ),
      class = "gaussianprocesses_rank_deficiency_error",
      info = name
    )
  }

  # A proper prior identifies the coefficients.
  model <- fit_gp(
    data$x,
    data$y,
    rbf_kernel(),
    noise_variance = 0.01,
    mean = basis_mean(collinear, coefficient_prior(0, 1))
  )
  expect_length(coef(model), 3L)
  expect_true(all(is.finite(predict_gp(model, data$x_new)$latent_variance)))
})


test_that("posterior draws include marginalized coefficient uncertainty", {
  data <- coefficient_data()
  model <- coefficient_fit(data, linear_mean(coefficient_prior()))
  points <- data$x_new[c(1, 3), ]
  prediction <- predict_gp(model, points)
  plug_in <- predict_gp(coefficient_fit(data, linear_mean()), points)
  n_draws <- 4000

  draws <- sample_gp_posterior(model, points, n_draws = n_draws, seed = 26)
  sample_mean <- rowMeans(draws$latent)
  sample_variance <- apply(draws$latent, 1, stats::var)

  mean_error <- sqrt(prediction$latent_variance / n_draws)
  variance_error <- prediction$latent_variance * sqrt(2 / (n_draws - 1))

  expect_true(all(abs(sample_mean - prediction$mean) < 5 * mean_error))
  expect_true(
    all(abs(sample_variance - prediction$latent_variance) < 5 * variance_error)
  )
  # The plug-in variance at the distant point is far outside that range.
  expect_gt(
    abs(sample_variance[2] - plug_in$latent_variance[2]),
    5 * variance_error[2]
  )
})


test_that("prior draws include a Gaussian coefficient prior", {
  x <- c(-1, 0, 2)
  kernel <- rbf_kernel(variance = 0.5, length_scale = 0.7)
  mean <- linear_mean(coefficient_prior(c(1, -0.5), c(2, 0.5)))
  n_draws <- 4000

  draws <- sample_gp_prior(x, kernel, mean = mean, n_draws = n_draws, seed = 7)
  expected_mean <- 1 - 0.5 * x
  expected_variance <- 0.5 + 2 + 0.5 * x^2

  expect_equal(draws$mean, expected_mean)
  expect_true(
    all(abs(rowMeans(draws$latent) - expected_mean) <
          5 * sqrt(expected_variance / n_draws))
  )
  expect_true(
    all(abs(apply(draws$latent, 1, stats::var) - expected_variance) <
          5 * expected_variance * sqrt(2 / (n_draws - 1)))
  )

  simulated <- simulate_gp_data(x, kernel, mean = mean, seed = 7)
  expect_length(simulated$latent, 3L)

  for (unknown in list(linear_mean(), linear_mean(coefficient_prior()))) {
    expect_error(
      sample_gp_prior(x, kernel, mean = unknown),
      "no value before a model is fitted",
      class = "gaussianprocesses_validation_error"
    )
  }
})


test_that("leave-one-out predictions match refitting without each point", {
  data <- coefficient_data(n = 15)
  means <- list(
    vague = linear_mean(coefficient_prior()),
    gaussian = linear_mean(coefficient_prior(c(1, 0, 0), c(2, 1, 1)))
  )

  for (name in names(means)) {
    loo <- loo_gp(coefficient_fit(data, means[[name]]))

    for (index in c(1, 7, 15)) {
      refit <- fit_gp(
        data$x[-index, ],
        data$y[-index],
        data$kernel,
        noise_variance = data$noise_variance,
        mean = means[[name]]
      )
      prediction <- predict_gp(refit, data$x[index, , drop = FALSE])

      expect_equal(loo$mean[index], prediction$mean, tolerance = 1e-8,
                   info = name)
      expect_equal(loo$variance[index], prediction$observation_variance,
                   tolerance = 1e-8, info = name)
    }
  }

  # Estimated coefficients stay at their full-data estimate.
  estimated <- coefficient_fit(data, linear_mean())
  loo <- loo_gp(estimated)
  refit <- fit_gp(
    data$x[-4, ],
    data$y[-4],
    data$kernel,
    noise_variance = data$noise_variance,
    mean = linear_mean(coef(estimated))
  )
  expect_equal(
    loo$mean[4],
    predict_gp(refit, data$x[4, , drop = FALSE])$mean,
    tolerance = 1e-8
  )
})


test_that("gradient predictions include parametric means", {
  data <- coefficient_data()
  means <- list(
    linear = linear_mean(),
    polynomial = polynomial_mean(2, coefficient_prior()),
    constant = constant_mean(coefficient_prior(0, 2))
  )
  points <- data$x_new[c(1, 2, 4), ]
  step <- 1e-5

  for (name in names(means)) {
    model <- coefficient_fit(data, means[[name]])
    gradient <- predict_gradient_gp(model, points)

    for (d in 1:2) {
      shift <- matrix(0, nrow(points), 2)
      shift[, d] <- step
      numerical <- (predict_gp(model, points + shift)$mean -
        predict_gp(model, points - shift)$mean) / (2 * step)

      expect_equal(unname(gradient$mean[, d]), numerical, tolerance = 1e-6,
                   info = name)

      # Variance of the central difference quotient, extrapolated to h = 0.
      quotient_variance <- function(h) {
        shift[, d] <- h
        vapply(
          seq_len(nrow(points)),
          function(i) {
            covariance <- predict_gp(
              model,
              rbind(points[i, ] + shift[i, ], points[i, ] - shift[i, ]),
              include_covariance = TRUE
            )$latent_covariance
            (covariance[1, 1] + covariance[2, 2] - 2 * covariance[1, 2]) /
              (4 * h^2)
          },
          numeric(1)
        )
      }
      limit <- (4 * quotient_variance(0.005) - quotient_variance(0.01)) / 3

      expect_lte(
        max(abs(gradient$variance[, d] - limit)) / max(gradient$variance[, d]),
        1e-5,
        label = paste(name, d)
      )
    }
  }

  expect_error(
    predict_gradient_gp(
      coefficient_fit(data, basis_mean(function(x) cbind(1, x[, 1]))),
      points
    ),
    "no derivatives",
    class = "gaussianprocesses_validation_error"
  )
})


test_that("sparse and heteroscedastic models need fixed coefficients", {
  data <- coefficient_data()

  for (mean in list(linear_mean(), linear_mean(coefficient_prior(0, 1)))) {
    expect_error(
      fit_sparse_gp(data$x, data$y, data$kernel, noise_variance = 0.02,
                    mean = mean, n_inducing = 5),
      "fixed coefficients only",
      class = "gaussianprocesses_validation_error"
    )
    expect_error(
      fit_heteroscedastic_gp(data$x, data$y, data$kernel,
                             noise_kernel = rbf_kernel(), mean = mean),
      "fixed coefficients only",
      class = "gaussianprocesses_validation_error"
    )
  }

  fixed <- polynomial_mean(1, coefficients = c(1, 2, -1))
  sparse <- fit_sparse_gp(data$x, data$y, data$kernel, noise_variance = 0.02,
                          mean = fixed, n_inducing = 8)
  expect_equal(sparse$mean$center, c(0.5, 0.5), tolerance = 0.1)
  expect_length(predict_sparse_gp(sparse, data$x_new)$mean, 4L)
})


test_that("coef(), vcov(), and logLik() report the coefficients", {
  data <- coefficient_data()

  zero <- coefficient_fit(data, zero_mean())
  expect_identical(coef(zero), stats::setNames(numeric(0), character(0)))
  expect_identical(dim(vcov(zero)), c(0L, 0L))

  constant <- coefficient_fit(data, constant_mean(2))
  expect_identical(coef(constant), c(intercept = 2))
  expect_identical(unname(vcov(constant)), matrix(0, 1, 1))

  fixed <- coefficient_fit(data, linear_mean(c(1, 2, -1)))
  expect_identical(coef(fixed), c(intercept = 1, x1 = 2, x2 = -1))
  expect_null(fixed$mean_fit)
  expect_equal(
    log_marginal_likelihood(fixed),
    gaussian_log_density(
      data$y,
      drop(cbind(1, data$x) %*% c(1, 2, -1)),
      dense_reference(data)$covariance
    ),
    tolerance = 1e-10
  )

  estimated <- coefficient_fit(data, linear_mean())
  expect_identical(attr(logLik(estimated), "df"), 3L)
  expect_identical(attr(logLik(estimated), "nobs"), 30L)

  restricted <- coefficient_fit(
    data,
    linear_mean(estimate_coefficients("reml"))
  )
  expect_identical(attr(logLik(restricted), "nobs"), 27L)

  vague <- coefficient_fit(data, linear_mean(coefficient_prior()))
  expect_identical(attr(logLik(vague), "df"), 0L)
  expect_identical(attr(logLik(vague), "nobs"), 27L)

  output <- utils::capture.output(print(estimated))
  expect_true(any(grepl("mean: LinearMean(estimated by ML)", output,
                        fixed = TRUE)))
  expect_true(any(grepl("mean coefficients: intercept=", output,
                        fixed = TRUE)))
})


test_that("models from schema 1 are read as the current schema", {
  data <- coefficient_data()
  model <- coefficient_fit(data, constant_mean(1))
  old <- model
  old$schema_version <- 1L
  old$mean_fit <- NULL

  upgraded <- gaussianprocesses:::.prepare_gp_model(old)

  expect_identical(
    upgraded$schema_version,
    gaussianprocesses:::.model_schema_version
  )
  expect_identical(upgraded$package_version, model$package_version)
  expect_equal(predict_gp(old, data$x_new), predict_gp(model, data$x_new))
})
