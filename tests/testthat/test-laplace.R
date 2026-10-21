# The Laplace approximation for latent Gaussian-process models.

laplace_data <- function() {
  set.seed(32)
  x <- seq(0, 5, length.out = 25)
  f <- sin(x)

  list(
    x = x,
    continuous = f + stats::rnorm(25, sd = 0.2),
    binary = as.numeric(f + stats::rnorm(25, sd = 0.3) > 0),
    counts = stats::rpois(25, 3 * exp(f)),
    exposure = rep(c(1, 2.5), length.out = 25)
  )
}


# Psi(f) = log p(y | f) - (f - m)' K^{-1} (f - m) / 2 and its gradient,
# written directly with solve(), independently of the Laplace code.
direct_objective <- function(likelihood, y, prior_mean, covariance, exposure) {
  precision <- solve(covariance)

  list(
    value = function(f) {
      sum(evaluate_likelihood(likelihood, y, f, exposure)$log_density) -
        0.5 * drop(crossprod(f - prior_mean, precision %*% (f - prior_mean)))
    },
    gradient = function(f) {
      evaluate_likelihood(likelihood, y, f, exposure)$first -
        drop(precision %*% (f - prior_mean))
    },
    precision = precision
  )
}


# log p(y | X) by adaptive tensor-product Gauss-Hermite quadrature: the nodes
# are placed at f_hat + sqrt(2) L t with L L' the Laplace covariance, so
#   int g(f) df = 2^(n/2) |L| sum_j w_j exp(|t_j|^2) g(f(t_j)).
quadrature_log_marginal <- function(model, order = 40) {
  rule <- .gauss_hermite(order)
  n <- model$n_observations
  covariance <- .latent_prior_covariance(model$kernel, model$mean, model$x)
  laplace <- solve(solve(covariance) + diag(model$curvature, n))
  factor <- t(chol(laplace))
  nodes <- as.matrix(expand.grid(rep(list(rule$nodes), n)))
  log_weights <- rowSums(as.matrix(expand.grid(
    rep(list(log(rule$weights) + rule$nodes^2), n)
  )))
  latent <- sweep(sqrt(2) * nodes %*% t(factor), 2, model$mode, "+")
  exposure <- if (is.null(model$exposure)) rep(1, n) else model$exposure
  entry <- .likelihood_entry(model$likelihood)
  log_likelihood <- rowSums(vapply(
    seq_len(n),
    function(i) {
      entry$log_density(
        model$y[i],
        latent[, i],
        model$likelihood$parameters,
        exposure[i]
      )
    },
    numeric(nrow(latent))
  ))
  prior_factor <- chol(covariance)
  whitened <- forwardsolve(
    t(prior_factor),
    t(sweep(latent, 2, model$prior_mean, "-"))
  )
  log_prior <- -0.5 * colSums(whitened^2) - sum(log(diag(prior_factor))) -
    0.5 * n * log(2 * pi)
  terms <- log_weights + log_likelihood + log_prior
  largest <- max(terms)

  largest + log(sum(exp(terms - largest))) + 0.5 * n * log(2) +
    sum(log(diag(factor)))
}


# Richardson-extrapolated central differences of a function of a vector.
richardson_gradient <- function(f, theta, step = 1e-4) {
  vapply(
    seq_along(theta),
    function(j) {
      shifted <- function(h) {
        up <- theta
        down <- theta
        up[j] <- up[j] + h
        down[j] <- down[j] - h
        (f(up) - f(down)) / (2 * h)
      }

      (4 * shifted(step / 2) - shifted(step)) / 3
    },
    numeric(1)
  )
}


test_that("with a Gaussian likelihood the Laplace model is the exact GP", {
  data <- laplace_data()
  kernel <- rbf_kernel(variance = 1.5, length_scale = 0.8)
  x_new <- c(-1, 0.3, 2.2, 7)
  means <- list(
    zero_mean(),
    constant_mean(0.3),
    linear_mean(coefficient_prior(c(0, 0.2), c(1, 0.5)))
  )

  for (mean in means) {
    exact <- fit_gp(data$x, data$continuous, kernel, noise_variance = 0.04,
                    mean = mean)
    latent <- fit_latent_gp(data$x, data$continuous, kernel,
                            gaussian_likelihood(0.04), mean = mean)
    label <- .format_mean(mean)

    # Psi is quadratic, so one Newton step reaches the mode.
    expect_identical(latent$laplace$iterations, 1L)
    expect_true(latent$laplace$converged)
    expect_lt(max(abs(latent$mode - fitted(exact))), 1e-10, label = label)
    expect_lt(
      abs(log_marginal_likelihood(latent) - log_marginal_likelihood(exact)),
      1e-10,
      label = label
    )
    expect_equal(
      unname(log_marginal_likelihood_gradient(latent)),
      unname(log_marginal_likelihood_gradient(exact)),
      tolerance = 1e-8,
      label = label
    )

    exact_prediction <- predict_gp(exact, x_new, include_covariance = TRUE)

    for (include_covariance in c(FALSE, TRUE)) {
      prediction <- predict_latent_gp(latent, x_new,
                                      include_covariance = include_covariance)
      expect_lt(max(abs(prediction$mean - exact_prediction$mean)), 1e-10)
      expect_lt(
        max(abs(prediction$latent_variance -
                  exact_prediction$latent_variance)),
        1e-10
      )
      expect_lt(
        max(abs(prediction$response_variance -
                  exact_prediction$observation_variance)),
        1e-10
      )
    }

    expect_lt(
      max(abs(prediction$latent_covariance -
                exact_prediction$latent_covariance)),
      1e-10
    )
  }
})


test_that("the mode and curvature agree with direct maximization of Psi", {
  set.seed(5)
  x <- c(0, 0.4, 1.1, 1.5, 2.3, 3, 3.2, 4.4)
  kernel <- rbf_kernel(variance = 2, length_scale = 0.9)
  cases <- list(
    list(likelihood = bernoulli_likelihood(), y = c(1, 1, 0, 1, 0, 0, 1, 0)),
    list(likelihood = bernoulli_likelihood("probit"),
         y = c(0, 1, 1, 1, 0, 0, 0, 1)),
    list(likelihood = poisson_likelihood(), y = c(0, 3, 5, 2, 9, 4, 1, 0),
         exposure = c(1, 1, 2, 2, 1, 3, 1, 1), mean = constant_mean(0.5))
  )

  for (case in cases) {
    mean <- if (is.null(case$mean)) zero_mean() else case$mean
    model <- fit_latent_gp(x, case$y, kernel, case$likelihood, mean = mean,
                           exposure = case$exposure)
    covariance <- evaluate_kernel(kernel, x)
    prior_mean <- evaluate_mean(mean, x)
    direct <- direct_objective(case$likelihood, case$y, prior_mean,
                               covariance, case$exposure)

    optimum <- stats::optim(
      prior_mean,
      function(f) -direct$value(f),
      function(f) -direct$gradient(f),
      method = "BFGS",
      control = list(reltol = 1e-16, maxit = 2000)
    )
    hessian <- stats::optimHess(
      optimum$par,
      function(f) -direct$value(f),
      function(f) -direct$gradient(f)
    )
    label <- utils::capture.output(print(case$likelihood))

    expect_lt(max(abs(model$mode - optimum$par)), 1e-6, label = label)
    expect_equal(
      hessian,
      direct$precision + diag(model$curvature),
      tolerance = 1e-6,
      label = label
    )
    expect_lt(model$laplace$gradient_norm, 1e-6)
  }
})


test_that("the analytical gradient passes finite differences and needs its implicit term", {
  data <- laplace_data()
  cases <- list(
    list(
      y = data$binary,
      likelihood = bernoulli_likelihood(),
      kernel = sum_kernel(rbf_kernel(1.5, 0.8), linear_kernel(0.1))
    ),
    list(
      y = data$binary,
      likelihood = bernoulli_likelihood("probit"),
      kernel = matern52_kernel(2, 1.2)
    ),
    list(
      y = data$counts,
      likelihood = poisson_likelihood(),
      kernel = rbf_kernel(0.8, 1.1),
      exposure = data$exposure,
      mean = constant_mean(0.5)
    ),
    list(
      y = data$continuous,
      likelihood = gaussian_likelihood(0.05),
      kernel = rbf_kernel(1.5, 0.8)
    )
  )

  for (case in cases) {
    mean <- if (is.null(case$mean)) zero_mean() else case$mean
    model <- fit_latent_gp(data$x, case$y, case$kernel, case$likelihood,
                           mean = mean, exposure = case$exposure)
    analytical <- log_marginal_likelihood_gradient(model)
    explicit <- .laplace_gradient(model, implicit = FALSE)
    n_kernel <- length(kernel_parameters(case$kernel, flatten = TRUE))

    log_marginal_at <- function(theta) {
      likelihood <- case$likelihood

      if (length(theta) > n_kernel) {
        likelihood$parameters$variance <- exp(theta[[n_kernel + 1L]])
      }

      log_marginal_likelihood(fit_latent_gp(
        data$x,
        case$y,
        update_kernel_parameters(case$kernel, exp(theta[seq_len(n_kernel)])),
        likelihood,
        mean = mean,
        exposure = case$exposure
      ))
    }

    numerical <- richardson_gradient(log_marginal_at, log(c(
      kernel_parameters(case$kernel, flatten = TRUE),
      unlist(case$likelihood$parameters)
    )))
    label <- utils::capture.output(print(case$likelihood))

    expect_identical(
      names(analytical),
      c(
        names(kernel_parameters(case$kernel, flatten = TRUE)),
        if (length(case$likelihood$parameters) > 0L) "likelihood.variance"
      )
    )
    expect_lt(
      max(abs(analytical - numerical) / pmax(abs(numerical), 1e-3)),
      1e-5,
      label = label
    )

    # Leaving out the implicit term breaks the check for every non-Gaussian
    # likelihood; with a Gaussian likelihood the third derivative is zero.
    explicit_error <- max(
      abs(explicit - numerical) / pmax(abs(numerical), 1e-3)
    )

    if (identical(case$likelihood$family, "gaussian")) {
      expect_lt(explicit_error, 1e-5)
    } else {
      expect_gt(explicit_error, 1e-2, label = label)
    }
  }
})


test_that("Laplace errors for n <= 3 match the documented values", {
  x <- c(0, 0.5, 2)
  kernel <- rbf_kernel(variance = 2, length_scale = 1)

  # The quadrature reference is exact for a Gaussian likelihood.
  gaussian <- fit_latent_gp(x, c(0.3, -0.2, 1), kernel,
                            gaussian_likelihood(0.5))
  expect_lt(
    abs(quadrature_log_marginal(gaussian, order = 5) -
          log_marginal_likelihood(gaussian)),
    1e-10
  )

  # The values in ?fit_latent_gp, for n = 1, 2, 3.
  documented <- list(
    list(likelihood = bernoulli_likelihood(), y = c(1, 1, 0),
         error = c(-0.0171, -0.0192, -0.0387)),
    list(likelihood = bernoulli_likelihood("probit"), y = c(1, 1, 0),
         error = c(-0.0196, -0.0298, -0.0450)),
    list(likelihood = poisson_likelihood(), y = c(0, 2, 7),
         error = c(0.0125, -0.0101, -0.0197)),
    list(likelihood = poisson_likelihood(), y = c(120, 150, 200),
         mean = constant_mean(5), error = c(-0.000681, -0.00118, -0.00158))
  )

  for (case in documented) {
    for (n in 1:3) {
      model <- fit_latent_gp(
        x[seq_len(n)],
        case$y[seq_len(n)],
        kernel,
        case$likelihood,
        mean = if (is.null(case$mean)) zero_mean() else case$mean
      )
      exact <- quadrature_log_marginal(model)

      # The quadrature has converged: 30 and 40 nodes agree.
      expect_lt(abs(quadrature_log_marginal(model, order = 30) - exact), 1e-7)
      expect_identical(
        signif(log_marginal_likelihood(model) - exact, 3),
        case$error[[n]],
        label = paste(utils::capture.output(print(case$likelihood)), n)
      )
    }
  }

  # One binary observation under a zero-mean prior has probability 1/2.
  single <- fit_latent_gp(0, 1, kernel, bernoulli_likelihood())
  expect_equal(quadrature_log_marginal(single), log(0.5), tolerance = 1e-12)
})


test_that("Newton iterations converge on linearly separable data", {
  x <- seq(1, 20)
  y <- as.numeric(x > 10)

  for (likelihood in list(bernoulli_likelihood(),
                          bernoulli_likelihood("probit"))) {
    expect_no_warning(
      model <- fit_latent_gp(x, y, rbf_kernel(variance = 25, length_scale = 3),
                             likelihood)
    )

    expect_true(model$laplace$converged)
    expect_identical(model$laplace$message, "converged")
    expect_lte(abs(model$laplace$decrement), 1e-14)
    expect_lt(model$laplace$gradient_norm, 1e-6)
    expect_true(all(sign(model$mode) == ifelse(y == 1, 1, -1)))
    expect_true(is.finite(log_marginal_likelihood(model)))
  }
})


test_that("the line search and warm starts reach the same mode", {
  data <- laplace_data()
  kernel <- rbf_kernel(variance = 3, length_scale = 0.5)
  counts <- stats::rpois(25, exp(5 + sin(data$x)))

  # Large counts make plain Newton steps overshoot from the prior mean.
  model <- fit_latent_gp(data$x, counts, kernel, poisson_likelihood())
  expect_true(model$laplace$converged)
  expect_gt(model$laplace$step_halvings, 0L)
  expect_lt(
    max(abs(model$mode - log(pmax(counts, 1)))),
    0.2
  )

  entry <- .likelihood_entry(poisson_likelihood())
  covariance <- evaluate_kernel(kernel, data$x)
  policy <- model$numerical_policy
  mode_from <- function(initial_alpha) {
    .laplace_mode(covariance, rep(0, 25), counts, entry, list(), rep(1, 25),
                  policy, 100L, 1e-14, initial_alpha = initial_alpha)
  }

  warm <- mode_from(model$alpha)
  expect_lte(warm$diagnostics$iterations, 1L)
  expect_lt(max(abs(warm$mode - model$mode)), 1e-8)
  expect_equal(warm$log_marginal_likelihood, model$log_marginal_likelihood,
               tolerance = 1e-12)

  # A warm start worse than the prior mean is ignored.
  poor <- mode_from(rep(50, 25))
  expect_identical(poor$diagnostics$iterations,
                   model$laplace$iterations)
  expect_equal(poor$mode, model$mode, tolerance = 1e-12)
})


test_that("Newton iterations that do not converge are reported", {
  data <- laplace_data()

  expect_warning(
    model <- fit_latent_gp(data$x, data$binary, rbf_kernel(4, 0.5),
                           bernoulli_likelihood(), max_iterations = 1),
    class = "gaussianprocesses_convergence_warning"
  )
  expect_false(model$laplace$converged)
  expect_identical(model$laplace$message, "iteration limit reached")
  expect_identical(model$laplace$iterations, 1L)
  expect_gt(model$laplace$decrement, 1e-14)
  expect_output(print(model), "Newton iterations: 1 \\(iteration limit")

  # A prior mean whose rate overflows leaves no finite starting point.
  expect_error(
    fit_latent_gp(data$x, data$counts, rbf_kernel(), poisson_likelihood(),
                  mean = constant_mean(800)),
    class = "gaussianprocesses_numerical_error"
  )
})


test_that("the factorization of B follows the jitter policy", {
  data <- laplace_data()
  model <- fit_latent_gp(data$x, data$binary, rbf_kernel(2, 1),
                         bernoulli_likelihood(), initial_jitter = 1e-6)
  covariance <- evaluate_kernel(rbf_kernel(2, 1), data$x)
  b_matrix <- diag(25) + outer(sqrt(model$curvature), sqrt(model$curvature)) *
    covariance

  expect_equal(model$jitter, 1e-6 * mean(diag(b_matrix)))
  expect_equal(crossprod(model$b_cholesky),
               b_matrix + diag(model$jitter, 25))
  expect_identical(fit_latent_gp(data$x, data$binary, rbf_kernel(2, 1),
                                 bernoulli_likelihood())$jitter, 0)
})


test_that("latent models validate their inputs", {
  data <- laplace_data()
  kernel <- rbf_kernel()

  for (call in list(
    function() fit_latent_gp(data$x, data$binary, kernel, "logit"),
    function() fit_latent_gp(data$x, data$binary[-1], kernel,
                             bernoulli_likelihood()),
    function() fit_latent_gp(data$x, data$counts, kernel,
                             bernoulli_likelihood()),
    function() fit_latent_gp(data$x, data$binary, kernel,
                             bernoulli_likelihood(), exposure = 2),
    function() fit_latent_gp(data$x, data$counts, kernel, poisson_likelihood(),
                             exposure = c(1, 2)),
    function() fit_latent_gp(data$x, data$binary, kernel,
                             bernoulli_likelihood(),
                             mean = linear_mean(estimate_coefficients("reml"))),
    function() fit_latent_gp(data$x, data$binary, kernel,
                             bernoulli_likelihood(),
                             mean = constant_mean(coefficient_prior())),
    function() fit_latent_gp(data$x, data$binary, kernel,
                             bernoulli_likelihood(), method = "ep"),
    function() fit_latent_gp(data$x, data$binary, kernel,
                             bernoulli_likelihood(), max_iterations = 0),
    function() fit_latent_gp(data$x, data$binary, kernel,
                             bernoulli_likelihood(), convergence_tolerance = 0),
    function() fit_latent_gp(cbind(data$x, data$x), data$binary,
                             select_dimensions(kernel, 3),
                             bernoulli_likelihood())
  )) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }

  # Factor responses keep their levels; the second level is coded 1.
  labels <- factor(ifelse(data$binary == 1, "yes", "no"),
                   levels = c("no", "yes"))
  model <- fit_latent_gp(data$x, labels, kernel, bernoulli_likelihood())
  expect_identical(model$response_levels, c("no", "yes"))
  expect_identical(model$y, data$binary)
  expect_null(model$exposure)
})


test_that("likelihoods that are not log-concave are refused", {
  original <- .likelihood_entry
  local_mocked_bindings(.likelihood_entry = function(likelihood) {
    entry <- original(likelihood)
    entry$log_concave <- FALSE
    entry
  })

  expect_error(
    fit_latent_gp(1:3, c(0, 1, 1), rbf_kernel(), bernoulli_likelihood()),
    class = "gaussianprocesses_log_concavity_error"
  )
})
