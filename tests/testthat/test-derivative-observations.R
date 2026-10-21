# Conditioning on observed derivatives.

derivative_data <- function(dimension = 1L) {
  if (dimension == 1L) {
    x <- c(-2, -1.1, 0, 0.9, 2.2, -1.6, 0.4, 1.5)
    types <- c(0, 0, 0, 0, 0, 1, 1, 1)
    f <- function(x) sin(1.5 * x[, 1])
    gradient <- function(x) cbind(1.5 * cos(1.5 * x[, 1]))
    x_new <- c(-0.7, 0.3, 1.2, 2.6)
  } else {
    x <- cbind(
      c(-1, 0.2, 1.1, -0.4, 0.7, 0, -0.8, 0.9, 0.3),
      c(0.5, -0.9, 0.3, 1, -0.2, 0, -0.5, 0.8, 1.2)
    )
    types <- c(0, 0, 0, 0, 0, 1, 2, 1, 2)
    f <- function(x) sin(x[, 1]) * cos(x[, 2])
    gradient <- function(x) {
      cbind(cos(x[, 1]) * cos(x[, 2]), -sin(x[, 1]) * sin(x[, 2]))
    }
    x_new <- cbind(c(-0.5, 0.4, 1.3), c(0.2, -0.6, 0.9))
  }

  x <- as.matrix(x)
  y <- f(x)
  derivative_rows <- types > 0
  y[derivative_rows] <- gradient(x)[cbind(
    which(derivative_rows),
    types[derivative_rows]
  )]

  list(
    x = x,
    y = y,
    types = types,
    x_new = as.matrix(x_new),
    noise = ifelse(types == 0, 1e-3, 4e-3)
  )
}


# Observations (f(x + h e_d) - f(x - h e_d)) / (2 h) in place of the
# derivatives, conditioned on with dense matrices.
difference_posterior <- function(kernel, data, h) {
  x <- data$x
  types <- data$types

  # Each observation is a combination of function values at points.
  points <- list()
  weights <- list()

  for (i in seq_len(nrow(x))) {
    if (types[i] == 0) {
      points[[i]] <- x[i, , drop = FALSE]
      weights[[i]] <- 1
    } else {
      shift <- matrix(0, 1, ncol(x))
      shift[types[i]] <- h
      points[[i]] <- rbind(x[i, ] + shift, x[i, ] - shift)
      weights[[i]] <- c(1, -1) / (2 * h)
    }
  }

  covariance <- function(a, b) {
    drop(crossprod(weights[[a]], evaluate_kernel(kernel, points[[a]],
                                                  points[[b]]) %*%
                     weights[[b]]))
  }
  n <- nrow(x)
  observation_covariance <- outer(
    seq_len(n),
    seq_len(n),
    Vectorize(covariance)
  ) + diag(data$noise)
  cross <- vapply(
    seq_len(n),
    function(a) {
      drop(evaluate_kernel(kernel, data$x_new, points[[a]]) %*% weights[[a]])
    },
    numeric(nrow(data$x_new))
  )
  cross <- matrix(cross, nrow(data$x_new))

  solved <- solve(observation_covariance, data$y)
  list(
    mean = drop(cross %*% solved),
    covariance = evaluate_kernel(kernel, data$x_new) -
      cross %*% solve(observation_covariance, t(cross)),
    log_likelihood = -0.5 * sum(data$y * solved) -
      0.5 * as.numeric(determinant(observation_covariance)$modulus) -
      0.5 * n * log(2 * pi)
  )
}


test_that("derivative observations are the limit of difference observations", {
  cases <- list(
    list(dimension = 1L, kernel = rbf_kernel(1.2, 0.8)),
    list(dimension = 1L, kernel = matern52_kernel(1, 1.1)),
    list(
      dimension = 2L,
      kernel = product_kernel(
        rbf_kernel(length_scale = c(0.9, 1.4)),
        periodic_kernel(variance = 0.8, period = 3, length_scale = 1.2)
      )
    )
  )

  for (case in cases) {
    data <- derivative_data(case$dimension)
    model <- fit_gp(
      data$x,
      data$y,
      case$kernel,
      noise_variance = data$noise,
      derivative = data$types
    )
    prediction <- predict_gp(model, data$x_new, include_covariance = TRUE,
                             observation_noise_variance = 0)

    posteriors <- lapply(
      c(0.02, 0.01, 0.005),
      function(h) difference_posterior(case$kernel, data, h)
    )
    fine <- posteriors[[3]]
    # The difference quotient has error O(h^2), plus O(h^3) for the
    # Matern 5/2 kernel; Richardson extrapolation removes both.
    extrapolate <- function(coarse, fine, power) {
      Map(function(a, b) (2^power * b - a) / (2^power - 1), coarse, fine)
    }
    limit <- extrapolate(
      extrapolate(posteriors[[1]], posteriors[[2]], 2),
      extrapolate(posteriors[[2]], posteriors[[3]], 2),
      3
    )
    label <- paste("dimension", case$dimension)

    expect_lte(
      max(abs(prediction$mean - limit$mean)) / max(abs(limit$mean)),
      1e-6,
      label = paste(label, "mean")
    )
    expect_lte(
      max(abs(prediction$latent_covariance - limit$covariance)) /
        max(abs(limit$covariance)),
      1e-6,
      label = paste(label, "covariance")
    )
    expect_lte(
      abs(log_marginal_likelihood(model) - limit$log_likelihood) /
        abs(limit$log_likelihood),
      1e-6,
      label = paste(label, "log likelihood")
    )
    # Without extrapolation the error is visibly O(h^2).
    expect_gt(
      max(abs(prediction$mean - fine$mean)),
      10 * max(abs(prediction$mean - limit$mean))
    )
  }
})


test_that("the likelihood gradient includes the derivative blocks", {
  data <- derivative_data(2L)
  kernels <- list(
    ard = rbf_kernel(1.3, c(0.8, 1.2)),
    matern = matern52_kernel(0.9, 1.1),
    periodic = periodic_kernel(1, c(0.9, 1.3), period = 2.5),
    product = product_kernel(
      rbf_kernel(length_scale = 1.4),
      linear_kernel(0.6)
    ),
    sum = sum_kernel(
      rational_quadratic_kernel(0.7, 1, alpha = 1.5),
      scale_kernel(matern32_kernel(length_scale = c(1.2, 0.9)), 0.4)
    ),
    additive = sum_kernel(
      select_dimensions(rbf_kernel(length_scale = 0.8), 1),
      select_dimensions(matern52_kernel(length_scale = 1.3), 2)
    )
  )
  step <- 1e-5

  for (name in names(kernels)) {
    kernel <- kernels[[name]]
    model <- fit_gp(data$x, data$y, kernel, noise_variance = 0.01,
                    derivative = data$types)
    analytical <- log_marginal_likelihood_gradient(model)
    parameters <- c(
      kernel_parameters(kernel, flatten = TRUE),
      noise_variance = 0.01
    )

    numerical <- vapply(
      names(parameters),
      function(path) {
        likelihood_at <- function(direction) {
          shifted <- parameters
          shifted[[path]] <- shifted[[path]] * exp(direction * step)
          log_marginal_likelihood(
            fit_gp(
              data$x,
              data$y,
              update_kernel_parameters(
                kernel,
                shifted[names(shifted) != "noise_variance"]
              ),
              noise_variance = shifted[["noise_variance"]],
              derivative = data$types
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


test_that("gradient observations reduce the posterior variance", {
  x <- c(-2, -1, 0, 1, 2)
  kernel <- rbf_kernel(length_scale = 0.8)
  values_only <- fit_gp(x, sin(x), kernel, noise_variance = 1e-4)
  with_gradients <- fit_gp(
    c(x, x),
    c(sin(x), cos(x)),
    kernel,
    noise_variance = 1e-4,
    derivative = rep(c(0, 1), each = 5)
  )
  grid <- seq(-2.5, 2.5, length.out = 41)
  before <- predict_gp(values_only, grid)
  after <- predict_gp(with_gradients, grid)

  expect_true(all(after$latent_variance <= before$latent_variance + 1e-12))
  # Between the observations the standard deviation falls by more than half.
  between <- abs(grid - round(grid)) > 0.4
  expect_true(
    all(after$latent_sd[between] < 0.5 * before$latent_sd[between])
  )
  expect_lt(
    max(abs(after$mean - sin(grid))[abs(grid) <= 2]),
    max(abs(before$mean - sin(grid))[abs(grid) <= 2]) / 2
  )

  covariance <- gaussianprocesses:::.model_observation_covariance(
    with_gradients
  )
  expect_equal(covariance, t(covariance), tolerance = 1e-14)
  expect_gte(
    min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values),
    -1e-10 * max(diag(covariance))
  )
})


test_that("derivative predictions and draws use derivative observations", {
  data <- derivative_data(2L)
  model <- fit_gp(data$x, data$y, rbf_kernel(1, c(0.9, 1.1)),
                  noise_variance = data$noise, derivative = data$types)
  gradient <- predict_gradient_gp(model, data$x_new)
  step <- 1e-5

  for (d in 1:2) {
    shift <- matrix(0, nrow(data$x_new), 2)
    shift[, d] <- step
    mean_at <- function(x) {
      predict_gp(model, x, observation_noise_variance = 0)$mean
    }
    numerical <- (mean_at(data$x_new + shift) -
      mean_at(data$x_new - shift)) / (2 * step)
    expect_equal(unname(gradient$mean[, d]), numerical, tolerance = 1e-6)

    quotient_variance <- function(h) {
      shift[, d] <- h
      vapply(
        seq_len(nrow(data$x_new)),
        function(i) {
          covariance <- predict_gp(
            model,
            rbind(data$x_new[i, ] + shift[i, ], data$x_new[i, ] - shift[i, ]),
            include_covariance = TRUE,
            observation_noise_variance = 0
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
      1e-5
    )
  }

  # At an observed derivative the posterior derivative is close to it.
  at_data <- predict_gradient_gp(model, data$x[data$types > 0, ])
  observed <- data$y[data$types > 0]
  types <- data$types[data$types > 0]
  expect_lt(
    max(abs(at_data$mean[cbind(seq_along(types), types)] - observed)),
    0.05
  )

  expect_equal(
    predict_gp(model, data$x_new, block_size = 1,
               observation_noise_variance = 0)$latent_variance,
    predict_gp(model, data$x_new,
               observation_noise_variance = 0)$latent_variance,
    tolerance = 1e-12
  )
  draws <- sample_gp_posterior(model, data$x_new, n_draws = 3, seed = 4)
  expect_identical(dim(draws$latent), c(3L, 3L))
})


test_that("means contribute their input derivatives", {
  data <- derivative_data(1L)
  kernel <- rbf_kernel(length_scale = 0.9)
  shifted <- data$y + ifelse(data$types == 0, 3 + 2 * data$x[, 1], 2)

  fixed <- fit_gp(data$x, shifted, kernel, noise_variance = data$noise,
                  mean = linear_mean(c(3, 2)), derivative = data$types)
  plain <- fit_gp(data$x, data$y, kernel, noise_variance = data$noise,
                  derivative = data$types)
  expect_equal(fixed$alpha, plain$alpha, tolerance = 1e-12)
  expect_equal(
    predict_gp(fixed, data$x_new, observation_noise_variance = 0)$mean,
    predict_gp(plain, data$x_new, observation_noise_variance = 0)$mean +
      3 + 2 * data$x_new[, 1],
    tolerance = 1e-12
  )
  expect_equal(
    unname(predict_gradient_gp(fixed, data$x_new)$mean[, 1]),
    unname(predict_gradient_gp(plain, data$x_new)$mean[, 1]) + 2,
    tolerance = 1e-12
  )

  # A constant mean has zero derivative.
  constant <- fit_gp(data$x, data$y + ifelse(data$types == 0, 5, 0), kernel,
                     noise_variance = data$noise, mean = constant_mean(5),
                     derivative = data$types)
  expect_equal(constant$alpha, plain$alpha, tolerance = 1e-12)

  # Estimated coefficients use derivative rows dh(x)' / dx of the basis.
  estimated <- fit_gp(data$x, shifted, kernel, noise_variance = data$noise,
                      mean = linear_mean(), derivative = data$types)
  basis <- cbind(1, ifelse(data$types == 0, data$x[, 1], 1))
  basis[data$types > 0, 1] <- 0
  covariance <- gaussianprocesses:::.model_observation_covariance(estimated)
  expect_equal(
    unname(coef(estimated)),
    drop(solve(
      crossprod(basis, solve(covariance, basis)),
      crossprod(basis, solve(covariance, shifted))
    )),
    tolerance = 1e-10
  )

  expect_error(
    fit_gp(data$x, data$y, kernel, mean = basis_mean(function(x) x, 1),
           derivative = data$types),
    "derivative observations",
    class = "gaussianprocesses_validation_error"
  )
})


test_that("leave-one-out and diagnostics cover every observation type", {
  data <- derivative_data(2L)
  kernel <- rbf_kernel(1, c(0.9, 1.1))
  model <- fit_gp(data$x, data$y, kernel, noise_variance = data$noise,
                  derivative = data$types)
  loo <- loo_gp(model)

  for (index in c(2, 6, 7)) {
    refit <- fit_gp(data$x[-index, ], data$y[-index], kernel,
                    noise_variance = data$noise[-index],
                    derivative = data$types[-index])
    type <- data$types[index]
    point <- data$x[index, , drop = FALSE]
    prediction <- if (type == 0) {
      p <- predict_gp(refit, point, observation_noise_variance = 0)
      c(p$mean, p$latent_variance)
    } else {
      p <- predict_gradient_gp(refit, point)
      unname(c(p$mean[1, type], p$variance[1, type]))
    }

    expect_equal(loo$mean[index], prediction[1], tolerance = 1e-8)
    expect_equal(loo$variance[index], prediction[2] + data$noise[index],
                 tolerance = 1e-8)
  }

  conditioning <- gp_numerical_diagnostics(model)$observation_types
  expect_identical(conditioning$type,
                   c("value", "derivative x1", "derivative x2"))
  expect_identical(conditioning$observations, c(5L, 2L, 2L))
  # Derivative variances are variance / length_scale^2.
  expect_equal(conditioning$covariance_scale[2:3],
               c(1 / 0.9^2, 1 / 1.1^2) + 4e-3, tolerance = 1e-12)
  expect_true(all(conditioning$reciprocal_condition_number > 0))

  plain <- gp_numerical_diagnostics(
    fit_gp(data$x, data$y, kernel, noise_variance = 0.01)
  )$observation_types
  expect_identical(plain$type, "value")

  output <- utils::capture.output(print(model))
  expect_true(any(grepl("observations: 9 (5 values, 4 derivatives)",
                        output, fixed = TRUE)))
})


test_that("optimize_gp() uses derivative observations", {
  data <- derivative_data(1L)
  fit <- function() {
    optimize_gp(data$x, data$y, rbf_kernel(), noise_variance = data$noise,
                derivative = data$types, n_starts = 2)
  }
  model <- fit()

  expect_identical(fit()$optimization, model$optimization)
  expect_true(model$optimization$converged)
  expect_identical(model$derivative, as.integer(data$types))
  expect_lt(max(abs(log_marginal_likelihood_gradient(model))), 1e-3)
  expect_gt(
    log_marginal_likelihood(model),
    log_marginal_likelihood(
      fit_gp(data$x, data$y, rbf_kernel(), noise_variance = data$noise,
             derivative = data$types)
    )
  )
})


test_that("derivative observations are validated", {
  data <- derivative_data(2L)

  for (derivative in list(
    data$types[-1],
    c(data$types[-1], 0.5),
    c(data$types[-1], -1),
    c(data$types[-1], 3),
    c(data$types[-1], NA),
    as.character(data$types)
  )) {
    expect_error(
      fit_gp(data$x, data$y, rbf_kernel(), derivative = derivative),
      "'derivative' must give one integer per observation",
      class = "gaussianprocesses_validation_error"
    )
  }

  for (kernel in list(
    matern12_kernel(),
    sum_kernel(rbf_kernel(), white_noise_kernel(0.1))
  )) {
    expect_error(
      fit_gp(data$x, data$y, kernel, derivative = data$types),
      class = "gaussianprocesses_smoothness_error"
    )
    expect_error(
      optimize_gp(data$x, data$y, kernel, derivative = data$types),
      class = "gaussianprocesses_smoothness_error"
    )
  }

  # Only function values: the same model as without the argument.
  values <- fit_gp(data$x, data$y, rbf_kernel(), noise_variance = 0.01,
                   derivative = rep(0, nrow(data$x)))
  expect_null(values$derivative)
  expect_identical(
    values$cholesky,
    fit_gp(data$x, data$y, rbf_kernel(), noise_variance = 0.01)$cholesky
  )
})


test_that("specifications and older schemas keep the observation types", {
  data <- derivative_data(1L)
  model <- fit_gp(data$x, data$y, rbf_kernel(), noise_variance = data$noise,
                  derivative = data$types)

  specification <- gaussianprocesses:::.gp_model_specification(model)
  expect_identical(specification$derivative, model$derivative)
  refit <- gaussianprocesses:::.fit_gp_from_specification(
    specification,
    model$x,
    model$y
  )
  expect_identical(refit$cholesky, model$cholesky)

  plain <- fit_gp(data$x, data$y, rbf_kernel(), noise_variance = 0.01)
  old <- plain
  old$schema_version <- 2L
  old$derivative <- NULL
  upgraded <- gaussianprocesses:::.prepare_gp_model(old)
  expect_identical(upgraded$schema_version, .model_schema_version)
  expect_identical(predict_gp(upgraded, data$x_new),
                   predict_gp(plain, data$x_new))
})
