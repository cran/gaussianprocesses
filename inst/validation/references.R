# Independent references for the inference methods of gaussianprocesses
#
# Each reference compares a quantity computed by the package with an
# independent computation of the same quantity:
#   implementation  another R package (DiceKriging, kernlab, nlme, gplite,
#                   or LAPACK through base R);
#   brute force     dense linear algebra, refitting, numerical integration,
#                   numerical differentiation, or a general-purpose
#                   optimizer applied to the defining formula;
#   closed form     a formula that the package does not use to compute the
#                   quantity, such as the exact GP that an approximation
#                   reduces to.
# A reference's run() returns its discrepancy, a non-negative number; the
# reference passes when the discrepancy is at most its tolerance. Unless a
# reference says otherwise, the discrepancy is the largest absolute
# difference divided by the largest absolute reference value.
#
# tests/testthat/test-validation-references.R runs every reference and
# requires it to pass, skipping references whose package is not installed.
# The validation article (vignettes/articles/validation.Rmd) runs them again
# and tabulates the discrepancies. `recorded_version` is the version of the
# reference package with which the tolerance was set.
#
# Tolerances depend on what limits the agreement: 1e-10 where both sides are
# exact up to rounding, which leaves room for other platforms and BLAS
# libraries; 1e-6 for first and 1e-4 for second finite differences; 1e-7
# for numerical integration against the package's quadrature; 1e-6 in the
# log likelihood for two optimizers' maxima; the convergence of iterative
# references (gplite, the heteroscedastic fixed point); and five Monte Carlo
# standard errors for sampling. Except for sampling, each tolerance is at
# least 30 times the discrepancy observed with the recorded versions.
#
# The file defines data and functions only. Source it with gaussianprocesses
# attached; it uses the package's exported functions and, where a
# reference needs a quantity that is not exported, gaussianprocesses:::.


# Methods and the exported functions they cover -------------------------------

# Every exported function that performs inference belongs to one method, and
# every method has at least one reference. The test file checks both.
validation_methods <- list(
  "Exact regression" = c(
    "fit_gp", "predict_gp", "fit_time_series_gp", "forecast_gp"
  ),
  "Marginal likelihood and gradient" = c(
    "log_marginal_likelihood", "log_marginal_likelihood_gradient"
  ),
  "Hyperparameter optimization" = c("optimize_gp", "optimize_time_series_gp"),
  "Mean coefficients" = c(
    "zero_mean", "constant_mean", "linear_mean", "polynomial_mean",
    "basis_mean", "estimate_coefficients", "coefficient_prior",
    "evaluate_mean"
  ),
  "Leave-one-out prediction" = c(
    "loo_gp", "gp_loo_scores", "gp_diagnostics", "gp_calibration_diagnostics"
  ),
  "Derivative observations and gradient prediction" = c(
    "predict_gradient_gp", "kernel_input_gradient"
  ),
  "Heteroscedastic regression" = c(
    "fit_heteroscedastic_gp", "predict_heteroscedastic_gp"
  ),
  "Sparse approximations (FITC and VFE)" = c(
    "fit_sparse_gp", "predict_sparse_gp", "optimize_sparse_gp",
    "select_inducing_points"
  ),
  "State-space inference" = character(0),
  "Laplace approximation" = c(
    "fit_latent_gp", "optimize_latent_gp", "predict_latent_gp",
    "gaussian_likelihood", "bernoulli_likelihood", "poisson_likelihood",
    "evaluate_likelihood", "likelihood_predictive"
  ),
  "Posterior and prior sampling" = c(
    "sample_gp_prior", "sample_gp_posterior", "simulate_gp_data"
  ),
  "Hyperparameter uncertainty" = c(
    "gp_hyperparameter_uncertainty", "gp_profile_likelihood"
  ),
  "Predictive scores" = c(
    "gp_scores", "gp_holdout_scores", "gp_classification_scores",
    "gp_count_scores", "forecast_metrics", "rolling_origin_gp"
  ),
  "Kernels" = c(
    "evaluate_kernel", "kernel_diagonal", "kernel_gradient",
    "kernel_parameters", "update_kernel_parameters", "rbf_kernel",
    "matern12_kernel", "matern32_kernel", "matern52_kernel",
    "rational_quadratic_kernel", "periodic_kernel", "linear_kernel",
    "white_noise_kernel", "spectral_mixture_kernel", "changepoint_kernel",
    "sum_kernel", "product_kernel", "scale_kernel", "select_dimensions",
    "time_series_kernel", "kernel_rbf", "kernel_matern12", "kernel_matern32",
    "kernel_matern52", "kernel_rational_quadratic", "kernel_periodic",
    "kernel_linear", "kernel_white_noise", "initialize_spectral_mixture"
  )
)

# Exported functions that perform no inference of their own.
validation_exclusions <- c(
  gp_time_index = "converts time values to a numeric index",
  gp_numerical_diagnostics = "reports conditioning and jitter",
  gp_simulation_scenario = "generates simulation data",
  benchmark_exact_gp_scaling = "times other functions",
  benchmark_sparse_gp = "times and compares other functions",
  run_gp_benchmark_suite = "times and compares other functions",
  compare_gp_predictions = "computes differences between predictions"
)


# Helpers ---------------------------------------------------------------------

validation_scaled_error <- function(actual, expected) {
  max(abs(actual - expected)) / max(abs(expected))
}


# Evaluates code() with a fixed seed and restores the caller's random state.
validation_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  set.seed(seed)
  code()
}


# Central differences with one Richardson extrapolation, error O(h^4).
validation_numerical_gradient <- function(f, theta, step = 1e-3) {
  vapply(seq_along(theta), function(j) {
    difference <- function(h) {
      shift <- replace(numeric(length(theta)), j, h)
      (f(theta + shift) - f(theta - shift)) / (2 * h)
    }
    (4 * difference(step / 2) - difference(step)) / 3
  }, numeric(1))
}


validation_data <- function(n, n_features = 1L, seed = 1L, noise_sd = 0.2) {
  validation_with_seed(seed, function() {
    x <- matrix(stats::runif(n * n_features, 0, 10), n, n_features)
    signal <- sin(x[, 1]) + 0.3 * cos(3 * x[, 1])

    if (n_features > 1L) {
      signal <- signal + 0.5 * sin(x[, 2] / 2)
    }

    list(
      x = if (n_features == 1L) sort(drop(x)) else x,
      y = if (n_features == 1L) {
        (signal + stats::rnorm(n, sd = noise_sd))[order(drop(x))]
      } else {
        signal + stats::rnorm(n, sd = noise_sd)
      },
      x_new = if (n_features == 1L) {
        seq(-1, 11, length.out = 25)
      } else {
        matrix(stats::runif(25 * n_features, -1, 11), 25, n_features)
      }
    )
  })
}


# Dense Gaussian conditioning: the posterior of f(x_new) given y with
# covariance K + diag(noise), from solve().
validation_dense_posterior <- function(kernel, x, y, noise, x_new,
                                       prior_mean = 0, new_mean = 0) {
  covariance <- evaluate_kernel(kernel, x) + diag(noise, length(y))
  cross <- evaluate_kernel(kernel, x, x_new)
  list(
    mean = new_mean + drop(crossprod(cross, solve(covariance, y - prior_mean))),
    covariance = evaluate_kernel(kernel, x_new) -
      crossprod(cross, solve(covariance, cross)),
    log_marginal_likelihood = -0.5 * sum((y - prior_mean) *
                                           solve(covariance, y - prior_mean)) -
      0.5 * as.numeric(determinant(covariance)$modulus) -
      0.5 * length(y) * log(2 * pi)
  )
}


# Kernel parameters on the log scale (all parameters used below are
# positive), and the kernel at given log-scale values.
validation_log_parameters <- function(kernel) {
  log(kernel_parameters(kernel, flatten = TRUE))
}

validation_kernel_at <- function(kernel, log_values) {
  update_kernel_parameters(kernel, exp(log_values))
}


validation_reference <- function(id, method, quantity, reference, type,
                                 tolerance, run, package = NA_character_,
                                 recorded_version = NA_character_) {
  list(
    id = id, method = method, quantity = quantity, reference = reference,
    type = type, tolerance = tolerance, run = run, package = package,
    recorded_version = recorded_version
  )
}


# References ------------------------------------------------------------------

validation_references <- function() {
  list(
    # Exact regression ------------------------------------------------------
    validation_reference(
      "exact-dicekriging",
      "Exact regression",
      "posterior mean and covariance at 25 inputs; ARD RBF kernel, known per-observation noise, fixed constant mean",
      "simple kriging with DiceKriging::km() and predict(type = \"SK\")",
      "implementation",
      tolerance = 1e-10,
      package = "DiceKriging",
      recorded_version = "1.6.1",
      run = function() {
        data <- validation_data(40, 2L, seed = 2L)
        noise <- seq(0.01, 0.05, length.out = 40)
        length_scale <- c(1.5, 2.5)
        model <- fit_gp(data$x, data$y, rbf_kernel(1.4, length_scale),
                        noise_variance = noise, mean = constant_mean(0.2))
        ours <- predict_gp(model, data$x_new, include_covariance = TRUE,
                           observation_noise_variance = 0)
        # DiceKriging's Gaussian covariance exp(-h^2 / (2 theta^2)) has
        # theta equal to the length scale.
        reference_model <- DiceKriging::km(
          ~1,
          design = data.frame(x1 = data$x[, 1], x2 = data$x[, 2]),
          response = data$y,
          covtype = "gauss",
          coef.trend = 0.2,
          coef.cov = length_scale,
          coef.var = 1.4,
          noise.var = noise,
          control = list(trace = FALSE)
        )
        reference <- stats::predict(
          reference_model,
          newdata = data.frame(x1 = data$x_new[, 1], x2 = data$x_new[, 2]),
          type = "SK",
          se.compute = TRUE,
          cov.compute = TRUE,
          light.return = FALSE
        )
        max(validation_scaled_error(ours$mean, reference$mean),
            validation_scaled_error(ours$latent_covariance, reference$cov))
      }
    ),
    validation_reference(
      "exact-kernlab",
      "Exact regression",
      "posterior mean and variance at 25 inputs; RBF kernel, known noise, zero mean",
      "kernlab::gausspr() with the rbfdot kernel",
      "implementation",
      tolerance = 1e-10,
      package = "kernlab",
      recorded_version = "0.9-33",
      run = function() {
        data <- validation_data(50, seed = 3L)
        model <- fit_gp(data$x, data$y, rbf_kernel(1, 1.2),
                        noise_variance = 0.04)
        ours <- predict_gp(model, data$x_new)
        # rbfdot is exp(-sigma |x - y|^2) with unit variance.
        reference_model <- kernlab::gausspr(
          matrix(data$x), data$y,
          type = "regression", scaled = FALSE, kernel = "rbfdot",
          kpar = list(sigma = 1 / (2 * 1.2^2)), var = 0.04,
          variance.model = TRUE
        )
        reference_mean <- as.numeric(kernlab::predict(
          reference_model, matrix(data$x_new), type = "response"
        ))
        reference_variance <- as.numeric(kernlab::predict(
          reference_model, matrix(data$x_new), type = "variance"
        ))
        max(validation_scaled_error(ours$mean, reference_mean),
            validation_scaled_error(ours$latent_variance, reference_variance))
      }
    ),
    validation_reference(
      "exact-dense",
      "Exact regression",
      "posterior mean and covariance at 25 inputs; composite ARD kernel on two inputs, per-observation noise, fixed linear mean",
      "Gaussian conditioning with dense matrices and solve()",
      "closed form",
      tolerance = 1e-10,
      run = function() {
        data <- validation_data(60, 2L, seed = 4L)
        kernel <- sum_kernel(
          matern52_kernel(1, c(2, 4)),
          product_kernel(
            select_dimensions(periodic_kernel(0.5, 1, 3), 1),
            select_dimensions(rbf_kernel(1, 5), 2)
          )
        )
        mean <- linear_mean(c(0.3, 0.05, -0.02))
        noise <- seq(0.02, 0.08, length.out = 60)
        model <- fit_gp(data$x, data$y, kernel, noise_variance = noise,
                        mean = mean)
        ours <- predict_gp(model, data$x_new, include_covariance = TRUE,
                           observation_noise_variance = 0)
        reference <- validation_dense_posterior(
          kernel, data$x, data$y, noise, data$x_new,
          evaluate_mean(mean, data$x), evaluate_mean(mean, data$x_new)
        )
        max(validation_scaled_error(ours$mean, reference$mean),
            validation_scaled_error(ours$latent_covariance,
                                    reference$covariance))
      }
    ),

    # Marginal likelihood and gradient --------------------------------------
    validation_reference(
      "likelihood-dense",
      "Marginal likelihood and gradient",
      "log marginal likelihood; composite ARD kernel, per-observation noise, fixed linear mean (relative error)",
      "the Gaussian log density with determinant() and solve()",
      "closed form",
      tolerance = 1e-10,
      run = function() {
        data <- validation_data(60, 2L, seed = 4L)
        kernel <- sum_kernel(matern32_kernel(1, c(2, 4)),
                             rational_quadratic_kernel(0.3, 1.5, 0.7))
        mean <- linear_mean(c(0.3, 0.05, -0.02))
        noise <- seq(0.02, 0.08, length.out = 60)
        model <- fit_gp(data$x, data$y, kernel, noise_variance = noise,
                        mean = mean)
        reference <- validation_dense_posterior(
          kernel, data$x, data$y, noise, data$x_new,
          evaluate_mean(mean, data$x)
        )$log_marginal_likelihood
        abs(log_marginal_likelihood(model) - reference) / abs(reference)
      }
    ),
    validation_reference(
      "likelihood-gradient",
      "Marginal likelihood and gradient",
      "gradient in the log kernel parameters and log noise variance; sum of RBF, periodic x RBF, and rational quadratic kernels",
      "Richardson-extrapolated central differences of log_marginal_likelihood()",
      "brute force",
      tolerance = 1e-6,
      run = function() {
        data <- validation_data(40, seed = 5L)
        kernel <- sum_kernel(
          rbf_kernel(1, 3),
          product_kernel(periodic_kernel(0.5, 1, 2), rbf_kernel(1, 8)),
          rational_quadratic_kernel(0.2, 0.7, 1.5)
        )
        theta <- c(validation_log_parameters(kernel),
                   noise_variance = log(0.05))
        likelihood <- function(values) {
          kernel_values <- values[-length(values)]
          names(kernel_values) <- names(theta)[-length(theta)]
          log_marginal_likelihood(fit_gp(
            data$x, data$y, validation_kernel_at(kernel, kernel_values),
            noise_variance = exp(values[[length(values)]])
          ))
        }
        gradient <- log_marginal_likelihood_gradient(fit_gp(
          data$x, data$y, kernel, noise_variance = 0.05
        ))
        validation_scaled_error(
          gradient[names(theta)],
          validation_numerical_gradient(likelihood, theta)
        )
      }
    ),

    # Hyperparameter optimization -------------------------------------------
    validation_reference(
      "optimization-dicekriging",
      "Hyperparameter optimization",
      "maximized profile log likelihood (absolute difference); RBF kernel with estimated variance, length scale, noise, and constant mean",
      "maximum likelihood with DiceKriging::km(nugget.estim = TRUE)",
      "implementation",
      tolerance = 1e-6,
      package = "DiceKriging",
      recorded_version = "1.6.1",
      run = function() {
        data <- validation_data(50, seed = 6L)
        model <- optimize_gp(
          data$x, data$y, rbf_kernel(1, 1), noise_variance = 0.1,
          mean = constant_mean(estimate_coefficients("ml")), n_starts = 3
        )
        # One start: several would need foreach, which DiceKriging only
        # suggests. Its random start reaches the same maximum for every seed
        # tried.
        reference <- validation_with_seed(1L, function() {
          DiceKriging::km(
            ~1, design = data.frame(x = data$x), response = data$y,
            covtype = "gauss", nugget.estim = TRUE,
            control = list(trace = FALSE)
          )
        })
        abs(log_marginal_likelihood(model) - reference@logLik)
      }
    ),
    validation_reference(
      "optimization-optim",
      "Hyperparameter optimization",
      "maximized log marginal likelihood (absolute difference); Matern-5/2 kernel with estimated variance, length scale, and noise",
      "stats::optim() (BFGS, numerical gradients, three starts) on the dense log density",
      "brute force",
      tolerance = 1e-6,
      run = function() {
        data <- validation_data(50, seed = 7L)
        model <- optimize_gp(data$x, data$y, matern52_kernel(1, 1),
                             noise_variance = 0.1, n_starts = 3)
        objective <- function(log_values) {
          -validation_dense_posterior(
            matern52_kernel(exp(log_values[[1]]), exp(log_values[[2]])),
            data$x, data$y, exp(log_values[[3]]), 0
          )$log_marginal_likelihood
        }
        optima <- lapply(list(c(0, 0, log(0.1)), c(1, 1, log(0.3)),
                              c(-1, -0.5, log(0.03))), function(start) {
          stats::optim(start, objective, method = "BFGS",
                       control = list(reltol = 1e-14, maxit = 1000))$value
        })
        abs(log_marginal_likelihood(model) + min(unlist(optima)))
      }
    ),

    # Mean coefficients -----------------------------------------------------
    validation_reference(
      "mean-nlme",
      "Mean coefficients",
      "maximum-likelihood and REML log likelihoods and coefficients of a linear mean; Gaussian correlation with a nugget",
      "generalized least squares with nlme::gls() and corGaus(fixed = TRUE)",
      "implementation",
      tolerance = 1e-10,
      package = "nlme",
      recorded_version = "3.1-168",
      run = function() {
        data <- validation_data(40, 2L, seed = 8L)
        frame <- data.frame(y = data$y, x1 = data$x[, 1], x2 = data$x[, 2])
        length_scale <- 2.5
        nugget <- 0.05
        errors <- vapply(c("ML", "REML"), function(method) {
          # corGaus is (1 - nugget) exp(-(d / range)^2), scaled by the
          # estimated sigma^2: range = sqrt(2) times the length scale.
          reference <- nlme::gls(
            y ~ x1 + x2, data = frame, method = method,
            correlation = nlme::corGaus(
              value = c(sqrt(2) * length_scale, nugget),
              form = ~ x1 + x2, nugget = TRUE, fixed = TRUE
            )
          )
          variance <- stats::sigma(reference)^2
          model <- fit_gp(
            data$x, data$y,
            rbf_kernel((1 - nugget) * variance, length_scale),
            noise_variance = nugget * variance,
            mean = linear_mean(estimate_coefficients(tolower(method)))
          )
          max(
            abs(log_marginal_likelihood(model) -
                  as.numeric(stats::logLik(reference))) /
              abs(as.numeric(stats::logLik(reference))),
            validation_scaled_error(unname(coef(model)),
                                    unname(stats::coef(reference)))
          )
        }, numeric(1))
        max(errors)
      }
    ),
    validation_reference(
      "mean-dicekriging",
      "Mean coefficients",
      "generalized least-squares coefficients, concentrated log likelihood, and universal-kriging mean and covariance (vague coefficient prior)",
      "DiceKriging::km() with a linear trend, logLikFun(), and predict(type = \"UK\")",
      "implementation",
      tolerance = 1e-10,
      package = "DiceKriging",
      recorded_version = "1.6.1",
      run = function() {
        data <- validation_data(40, 2L, seed = 9L)
        length_scale <- c(2, 3)
        design <- data.frame(x1 = data$x[, 1], x2 = data$x[, 2])
        reference <- DiceKriging::km(
          ~ x1 + x2, design = design, response = data$y, covtype = "gauss",
          coef.cov = length_scale, coef.var = 1.5,
          noise.var = rep(0.04, 40), control = list(trace = FALSE)
        )
        kernel <- rbf_kernel(1.5, length_scale)
        estimated <- fit_gp(data$x, data$y, kernel, noise_variance = 0.04,
                            mean = linear_mean())
        # km() leaves the likelihood empty when every parameter is given;
        # logLikFun() evaluates it, concentrated on the trend.
        reference@method <- "MLE"
        reference@case <- "LLconcentration_beta"
        reference_likelihood <- as.numeric(
          DiceKriging::logLikFun(c(length_scale, 1.5), reference)
        )
        universal <- stats::predict(
          reference,
          newdata = data.frame(x1 = data$x_new[, 1], x2 = data$x_new[, 2]),
          type = "UK", se.compute = TRUE, cov.compute = TRUE,
          light.return = FALSE
        )
        vague <- predict_gp(
          fit_gp(data$x, data$y, kernel, noise_variance = 0.04,
                 mean = linear_mean(coefficient_prior())),
          data$x_new, include_covariance = TRUE
        )
        max(
          validation_scaled_error(unname(coef(estimated)),
                                  unname(reference@trend.coef)),
          abs(log_marginal_likelihood(estimated) - reference_likelihood) /
            abs(reference_likelihood),
          validation_scaled_error(vague$mean, universal$mean),
          validation_scaled_error(vague$latent_covariance, universal$cov)
        )
      }
    ),
    validation_reference(
      "mean-gaussian-prior",
      "Mean coefficients",
      "posterior mean, covariance, and log marginal likelihood with a Gaussian prior on the coefficients of a quadratic mean",
      "the exact GP with covariance k(x, x') + h(x)' B h(x') and mean h(x)' b, by dense conditioning",
      "closed form",
      tolerance = 1e-10,
      run = function() {
        data <- validation_data(40, seed = 10L)
        prior_mean <- c(0.5, -0.1, 0.02)
        prior_covariance <- diag(c(2, 0.5, 0.1))
        kernel <- matern32_kernel(0.8, 1.5)
        mean <- basis_mean(
          function(x) cbind(1, x[, 1], x[, 1]^2),
          coefficient_prior(prior_mean, prior_covariance)
        )
        model <- fit_gp(data$x, data$y, kernel, noise_variance = 0.04,
                        mean = mean)
        ours <- predict_gp(model, data$x_new, include_covariance = TRUE)
        basis <- cbind(1, data$x, data$x^2)
        new_basis <- cbind(1, data$x_new, data$x_new^2)
        covariance <- evaluate_kernel(kernel, data$x) +
          basis %*% prior_covariance %*% t(basis) + diag(0.04, 40)
        cross <- evaluate_kernel(kernel, data$x, data$x_new) +
          basis %*% prior_covariance %*% t(new_basis)
        residual <- data$y - drop(basis %*% prior_mean)
        reference_mean <- drop(new_basis %*% prior_mean) +
          drop(crossprod(cross, solve(covariance, residual)))
        reference_covariance <- evaluate_kernel(kernel, data$x_new) +
          new_basis %*% prior_covariance %*% t(new_basis) -
          crossprod(cross, solve(covariance, cross))
        reference_likelihood <- -0.5 * sum(residual *
                                             solve(covariance, residual)) -
          0.5 * as.numeric(determinant(covariance)$modulus) -
          20 * log(2 * pi)
        max(
          validation_scaled_error(ours$mean, reference_mean),
          validation_scaled_error(ours$latent_covariance,
                                  reference_covariance),
          abs(log_marginal_likelihood(model) - reference_likelihood) /
            abs(reference_likelihood)
        )
      }
    ),

    # Leave-one-out prediction ----------------------------------------------
    validation_reference(
      "loo-refit",
      "Leave-one-out prediction",
      "leave-one-out predictive means and variances of all 30 observations, with marginalized linear-mean coefficients (vague prior), which each refit re-estimates",
      "30 refits, each without one observation, predicting it with fit_gp() and predict_gp()",
      "brute force",
      tolerance = 1e-10,
      run = function() {
        data <- validation_data(30, seed = 11L)
        kernel <- matern52_kernel(1, 1.5)
        mean <- linear_mean(coefficient_prior())
        model <- fit_gp(data$x, data$y, kernel, noise_variance = 0.05,
                        mean = mean)
        loo <- loo_gp(model)
        refits <- vapply(seq_along(data$y), function(i) {
          refit <- fit_gp(data$x[-i], data$y[-i], kernel,
                          noise_variance = 0.05, mean = mean)
          prediction <- predict_gp(refit, data$x[i])
          c(prediction$mean, prediction$observation_variance)
        }, numeric(2))
        max(validation_scaled_error(data$y - loo$residual, refits[1, ]),
            validation_scaled_error(loo$variance, refits[2, ]))
      }
    ),

    # Derivative observations and gradient prediction ------------------------
    validation_reference(
      "derivative-observations",
      "Derivative observations and gradient prediction",
      "posterior mean and covariance of function values given 12 values and 6 derivative observations; RBF + Matern-5/2 kernel",
      "dense conditioning on a joint covariance built from central differences of evaluate_kernel() (step 1e-4)",
      "brute force",
      tolerance = 1e-6,
      run = function() {
        validation_with_seed(12L, function() {
          kernel <- sum_kernel(rbf_kernel(1, 1.2), matern52_kernel(0.3, 2))
          x_values <- sort(stats::runif(12, 0, 6))
          x_slopes <- sort(stats::runif(6, 0, 6))
          y <- c(sin(x_values), cos(x_slopes)) +
            stats::rnorm(18, sd = 0.05)
          x_new <- seq(-0.5, 6.5, length.out = 15)
          model <- fit_gp(c(x_values, x_slopes), y, kernel,
                          noise_variance = 0.0025,
                          derivative = rep(c(0L, 1L), c(12, 6)))
          ours <- predict_gp(model, x_new, include_covariance = TRUE)

          h <- 1e-4
          k <- function(a, b) evaluate_kernel(kernel, a, b)
          value_slope <- function(a, b) (k(a, b + h) - k(a, b - h)) / (2 * h)
          slope_slope <- (k(x_slopes + h, x_slopes + h) -
                            k(x_slopes + h, x_slopes - h) -
                            k(x_slopes - h, x_slopes + h) +
                            k(x_slopes - h, x_slopes - h)) / (4 * h^2)
          covariance <- rbind(
            cbind(k(x_values, x_values), value_slope(x_values, x_slopes)),
            cbind(t(value_slope(x_values, x_slopes)), slope_slope)
          ) + diag(0.0025, 18)
          cross <- rbind(k(x_values, x_new), t(value_slope(x_new, x_slopes)))
          reference_mean <- drop(crossprod(cross, solve(covariance, y)))
          reference_covariance <- k(x_new, x_new) -
            crossprod(cross, solve(covariance, cross))
          max(validation_scaled_error(ours$mean, reference_mean),
              validation_scaled_error(ours$latent_covariance,
                                      reference_covariance))
        })
      }
    ),
    validation_reference(
      "gradient-prediction",
      "Derivative observations and gradient prediction",
      "posterior mean of the gradient at 25 inputs on two dimensions; ARD Matern-5/2 kernel",
      "Richardson-extrapolated central differences of the posterior mean from predict_gp()",
      "brute force",
      tolerance = 1e-6,
      run = function() {
        data <- validation_data(50, 2L, seed = 13L)
        model <- fit_gp(data$x, data$y, matern52_kernel(1, c(2, 3)),
                        noise_variance = 0.04)
        ours <- predict_gradient_gp(model, data$x_new)$mean
        reference <- t(vapply(seq_len(nrow(data$x_new)), function(i) {
          validation_numerical_gradient(function(point) {
            predict_gp(model, matrix(point, 1))$mean
          }, data$x_new[i, ], step = 1e-2)
        }, numeric(2)))
        validation_scaled_error(ours, reference)
      }
    ),

    # Heteroscedastic regression --------------------------------------------
    validation_reference(
      "heteroscedastic-fixed-point",
      "Heteroscedastic regression",
      "estimated noise variances (largest difference of log variances): the documented update applied to the returned estimate reproduces it",
      "one step of the update recomputed with leave-one-out refits and dense conditioning of the log-noise GP",
      "brute force",
      tolerance = 1e-8,
      run = function() {
        simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 60,
                                             seed = 14L)
        x <- drop(simulation$x)
        y <- simulation$observed
        kernel <- matern52_kernel(1, 0.2)
        noise_kernel <- rbf_kernel(1, 0.4)
        model <- fit_heteroscedastic_gp(
          x, y, kernel, noise_kernel, max_iterations = 500,
          convergence_tolerance = 1e-10
        )
        noise <- model$training_noise_variance
        loo <- vapply(seq_along(y), function(i) {
          posterior <- validation_dense_posterior(kernel, x[-i], y[-i],
                                                  noise[-i], x[i])
          c(y[[i]] - posterior$mean, posterior$covariance + noise[[i]])
        }, numeric(2))
        # log(e^2) - E[log chi^2_1] - log((v + s^2) / s^2), smoothed by a GP
        # with noise variance Var[log chi^2_1] = trigamma(1/2) and a constant
        # mean equal to the targets' average.
        target <- log(pmax(loo[1, ]^2, 1e-6)) - digamma(0.5) - log(2) -
          log(loo[2, ] / noise)
        smoothed <- validation_dense_posterior(
          noise_kernel, x, target, trigamma(0.5), x,
          mean(target), mean(target)
        )$mean
        max(abs(log(pmax(exp(smoothed), 1e-6)) - log(noise)))
      }
    ),

    # Sparse approximations -------------------------------------------------
    validation_reference(
      "fitc-dense",
      "Sparse approximations (FITC and VFE)",
      "FITC log marginal likelihood (relative error), predictive mean, and covariance; 60 inputs, 9 inducing points",
      "the FITC density N(y; 0, Q + diag(K - Q) + sigma^2 I) and predictive equations with dense matrices",
      "closed form",
      tolerance = 1e-10,
      run = function() {
        data <- validation_data(60, 2L, seed = 15L)
        kernel <- rbf_kernel(1.2, c(2, 3))
        model <- fit_sparse_gp(data$x, data$y, kernel, noise_variance = 0.04,
                               n_inducing = 9, method = "fitc")
        z <- model$inducing_points
        kuu <- evaluate_kernel(kernel, z)
        kuf <- evaluate_kernel(kernel, z, data$x)
        low_rank <- crossprod(kuf, solve(kuu, kuf))
        lambda <- diag(evaluate_kernel(kernel, data$x)) - diag(low_rank) + 0.04
        covariance <- low_rank + diag(lambda)
        likelihood <- -0.5 * sum(data$y * solve(covariance, data$y)) -
          0.5 * as.numeric(determinant(covariance)$modulus) -
          30 * log(2 * pi)
        ksu <- evaluate_kernel(kernel, data$x_new, z)
        sigma <- solve(kuu + kuf %*% (t(kuf) / lambda))
        reference_mean <- drop(ksu %*% sigma %*% kuf %*% (data$y / lambda))
        reference_covariance <- evaluate_kernel(kernel, data$x_new) -
          ksu %*% solve(kuu, t(ksu)) + ksu %*% sigma %*% t(ksu)
        ours <- predict_sparse_gp(model, data$x_new, include_covariance = TRUE)
        max(
          abs(model$log_marginal_likelihood - likelihood) / abs(likelihood),
          validation_scaled_error(ours$mean, reference_mean),
          validation_scaled_error(ours$latent_covariance, reference_covariance)
        )
      }
    ),
    validation_reference(
      "vfe-dense",
      "Sparse approximations (FITC and VFE)",
      "VFE bound (relative error), predictive mean, and covariance; 60 inputs, 9 inducing points",
      "Titsias's bound log N(y; 0, Q + sigma^2 I) - tr(K - Q) / (2 sigma^2) and the predictive equations with dense matrices",
      "closed form",
      tolerance = 1e-10,
      run = function() {
        data <- validation_data(60, 2L, seed = 16L)
        kernel <- matern52_kernel(1.2, c(2, 3))
        model <- fit_sparse_gp(data$x, data$y, kernel, noise_variance = 0.04,
                               n_inducing = 9, method = "vfe")
        z <- model$inducing_points
        kuu <- evaluate_kernel(kernel, z)
        kuf <- evaluate_kernel(kernel, z, data$x)
        low_rank <- crossprod(kuf, solve(kuu, kuf))
        covariance <- low_rank + diag(0.04, 60)
        bound <- -0.5 * sum(data$y * solve(covariance, data$y)) -
          0.5 * as.numeric(determinant(covariance)$modulus) -
          30 * log(2 * pi) -
          0.5 * sum(diag(evaluate_kernel(kernel, data$x)) - diag(low_rank)) /
          0.04
        ksu <- evaluate_kernel(kernel, data$x_new, z)
        sigma <- solve(kuu + tcrossprod(kuf) / 0.04)
        reference_mean <- drop(ksu %*% sigma %*% kuf %*% data$y) / 0.04
        reference_covariance <- evaluate_kernel(kernel, data$x_new) -
          ksu %*% solve(kuu, t(ksu)) + ksu %*% sigma %*% t(ksu)
        ours <- predict_sparse_gp(model, data$x_new, include_covariance = TRUE)
        max(
          abs(model$log_marginal_likelihood - bound) / abs(bound),
          validation_scaled_error(ours$mean, reference_mean),
          validation_scaled_error(ours$latent_covariance, reference_covariance)
        )
      }
    ),
    validation_reference(
      "sparse-exact-limit",
      "Sparse approximations (FITC and VFE)",
      "FITC and VFE log marginal likelihoods (relative error) and predictions when every input is an inducing point",
      "the exact GP, to which both approximations reduce when Z = X",
      "closed form",
      tolerance = 1e-10,
      run = function() {
        data <- validation_data(30, seed = 17L)
        kernel <- matern32_kernel(1, 1.5)
        exact <- fit_gp(data$x, data$y, kernel, noise_variance = 0.04)
        expected <- predict_gp(exact, data$x_new)
        errors <- vapply(c("fitc", "vfe"), function(method) {
          model <- fit_sparse_gp(data$x, data$y, kernel, noise_variance = 0.04,
                                 inducing_points = data$x, method = method)
          prediction <- predict_sparse_gp(model, data$x_new)
          max(
            abs(model$log_marginal_likelihood -
                  log_marginal_likelihood(exact)) /
              abs(log_marginal_likelihood(exact)),
            validation_scaled_error(prediction$mean, expected$mean),
            validation_scaled_error(prediction$latent_variance,
                                    expected$latent_variance)
          )
        }, numeric(1))
        max(errors)
      }
    ),
    validation_reference(
      "sparse-gradient",
      "Sparse approximations (FITC and VFE)",
      "gradient of the VFE bound and the FITC objective in the log hyperparameters and the 6 x 2 inducing inputs",
      "Richardson-extrapolated central differences of the objective, refitting with fit_sparse_gp()",
      "brute force",
      tolerance = 1e-6,
      run = function() {
        data <- validation_data(40, 2L, seed = 18L)
        kernel <- rbf_kernel(1.2, c(2, 3))
        z <- data$x[c(3, 9, 15, 22, 30, 37), ]
        errors <- vapply(c("vfe", "fitc"), function(method) {
          objective <- function(values) {
            fit_sparse_gp(
              data$x, data$y,
              rbf_kernel(exp(values[[1]]), exp(values[2:3])),
              noise_variance = exp(values[[4]]),
              inducing_points = matrix(values[-(1:4)], 6),
              method = method
            )$log_marginal_likelihood
          }
          theta <- c(log(c(1.2, 2, 3, 0.05)), as.vector(z))
          model <- fit_sparse_gp(data$x, data$y, kernel, noise_variance = 0.05,
                                 inducing_points = z, method = method)
          gradient <- gaussianprocesses:::.sparse_gradient(
            model, include_noise = TRUE, include_inducing = TRUE
          )
          validation_scaled_error(
            unname(gradient),
            validation_numerical_gradient(objective, theta)
          )
        }, numeric(1))
        max(errors)
      }
    ),
    validation_reference(
      "inducing-greedy",
      "Sparse approximations (FITC and VFE)",
      "the first 10 points chosen by greedy variance selection (number of differences)",
      "the pivot order of LAPACK's pivoted Cholesky decomposition, chol(pivot = TRUE)",
      "implementation",
      tolerance = 0,
      run = function() {
        data <- validation_data(60, 2L, seed = 19L)
        kernel <- matern32_kernel(1, c(1.5, 2))
        selected <- select_inducing_points(data$x, 10, method = "variance",
                                           kernel = kernel)
        pivot <- attr(suppressWarnings(
          chol(evaluate_kernel(kernel, data$x), pivot = TRUE)
        ), "pivot")
        sum(attr(selected, "indices") != pivot[1:10])
      }
    ),

    # State-space inference -------------------------------------------------
    validation_reference(
      "state-space-exact",
      "State-space inference",
      "log marginal likelihood (relative error), posterior means, and variances at 25 inputs; Matern-5/2 + scaled Matern-1/2 kernel on 200 irregular inputs with ties",
      "the exact GP with dense matrices (fit_gp() and predict_gp())",
      "closed form",
      tolerance = 1e-10,
      run = function() {
        data <- validation_data(200, seed = 20L)
        x <- data$x
        x[51:53] <- x[[50]]
        kernel <- sum_kernel(matern52_kernel(1, 1.5),
                             scale_kernel(matern12_kernel(1, 0.5), 0.2))
        exact <- fit_gp(x, data$y, kernel, noise_variance = 0.04)
        kalman <- gaussianprocesses:::.fit_state_space_model(
          x, data$y, kernel, 0.04
        )
        expected <- predict_gp(exact, data$x_new)
        actual <- predict_gp(kalman, data$x_new)
        max(
          abs(kalman$log_marginal_likelihood -
                log_marginal_likelihood(exact)) /
            abs(log_marginal_likelihood(exact)),
          validation_scaled_error(actual$mean, expected$mean),
          validation_scaled_error(actual$latent_variance,
                                  expected$latent_variance)
        )
      }
    ),

    # Laplace approximation -------------------------------------------------
    validation_reference(
      "laplace-gplite-bernoulli",
      "Laplace approximation",
      "posterior mode, approximate log marginal likelihood (relative error), and latent predictive means and variances (absolute differences); logit and probit links",
      "gplite with approx_laplace() and the same kernel and jitter",
      "implementation",
      tolerance = 1e-5,
      package = "gplite",
      recorded_version = "0.13.0",
      run = function() {
        validation_with_seed(21L, function() {
          x <- sort(stats::runif(40, 0, 6))
          y <- stats::rbinom(40, 1, stats::plogis(2 * sin(x)))
          x_new <- c(-0.5, 1.3, 3.7, 8)
          # gplite adds `jitter` to the diagonal of K, which a white-noise
          # kernel reproduces; its magnitude is the standard deviation.
          jitter <- 1e-9
          errors <- vapply(c("logit", "probit"), function(link) {
            reference <- gplite::gp_init(
              cfs = gplite::cf_sexp(lscale = 0.9, magn = sqrt(2.5)),
              lik = gplite::lik_bernoulli(link),
              approx = gplite::approx_laplace(maxiter = 200, tol = 1e-12)
            )
            reference <- gplite::gp_fit(reference, matrix(x), y,
                                        jitter = jitter)
            model <- fit_latent_gp(
              x, y, sum_kernel(rbf_kernel(2.5, 0.9), white_noise_kernel(jitter)),
              bernoulli_likelihood(link)
            )
            reference_prediction <- gplite::gp_pred(reference, matrix(x_new),
                                                    var = TRUE, jitter = jitter)
            prediction <- predict_latent_gp(model, x_new)
            max(
              abs(model$mode - reference$fit$fmean),
              abs(model$log_marginal_likelihood +
                    gplite::gp_energy(reference, include_prior = FALSE)) /
                abs(model$log_marginal_likelihood),
              abs(prediction$mean - reference_prediction$mean),
              abs(prediction$latent_variance - reference_prediction$var)
            )
          }, numeric(1))
          max(errors)
        })
      }
    ),
    validation_reference(
      "laplace-gplite-poisson",
      "Laplace approximation",
      "posterior mode, approximate log marginal likelihood (relative error), and latent predictive means and variances (absolute differences); Poisson counts with exposure",
      "gplite with approx_laplace() and an offset",
      "implementation",
      tolerance = 1e-5,
      package = "gplite",
      recorded_version = "0.13.0",
      run = function() {
        validation_with_seed(22L, function() {
          x <- sort(stats::runif(40, 0, 10))
          exposure <- stats::runif(40, 0.5, 2)
          y <- stats::rpois(40, exposure * exp(1 + sin(x)))
          x_new <- c(-0.5, 2.2, 6.1, 11)
          jitter <- 1e-9
          reference <- gplite::gp_init(
            cfs = gplite::cf_sexp(lscale = 0.9, magn = sqrt(1.5)),
            lik = gplite::lik_poisson(),
            approx = gplite::approx_laplace(maxiter = 200, tol = 1e-12)
          )
          reference <- gplite::gp_fit(reference, matrix(x), y,
                                      offset = log(exposure), jitter = jitter)
          model <- fit_latent_gp(
            x, y, sum_kernel(rbf_kernel(1.5, 0.9), white_noise_kernel(jitter)),
            poisson_likelihood(), exposure = exposure
          )
          reference_prediction <- gplite::gp_pred(reference, matrix(x_new),
                                                  var = TRUE, jitter = jitter)
          prediction <- predict_latent_gp(model, x_new)
          max(
            abs(model$mode - reference$fit$fmean),
            abs(model$log_marginal_likelihood +
                  gplite::gp_energy(reference, include_prior = FALSE)) /
              abs(model$log_marginal_likelihood),
            abs(prediction$mean - reference_prediction$mean),
            abs(prediction$latent_variance - reference_prediction$var)
          )
        })
      }
    ),
    validation_reference(
      "laplace-gaussian-exact",
      "Laplace approximation",
      "log marginal likelihood (relative error), latent predictive means, and variances with a Gaussian likelihood",
      "the exact GP, which the Laplace approximation equals for Gaussian observations",
      "closed form",
      tolerance = 1e-10,
      run = function() {
        data <- validation_data(40, seed = 23L)
        kernel <- matern52_kernel(1, 1.5)
        exact <- fit_gp(data$x, data$y, kernel, noise_variance = 0.04)
        latent <- fit_latent_gp(data$x, data$y, kernel,
                                gaussian_likelihood(0.04))
        expected <- predict_gp(exact, data$x_new)
        actual <- predict_latent_gp(latent, data$x_new)
        max(
          abs(latent$log_marginal_likelihood -
                log_marginal_likelihood(exact)) /
            abs(log_marginal_likelihood(exact)),
          validation_scaled_error(actual$mean, expected$mean),
          validation_scaled_error(actual$latent_variance,
                                  expected$latent_variance)
        )
      }
    ),
    validation_reference(
      "laplace-gradient",
      "Laplace approximation",
      "gradient of the approximate log marginal likelihood in the log kernel parameters; Bernoulli (logit) and Poisson likelihoods",
      "Richardson-extrapolated central differences of the approximate log marginal likelihood, refitting with fit_latent_gp()",
      "brute force",
      tolerance = 1e-6,
      run = function() {
        validation_with_seed(24L, function() {
          x <- sort(stats::runif(30, 0, 6))
          cases <- list(
            list(y = stats::rbinom(30, 1, stats::plogis(2 * sin(x))),
                 likelihood = bernoulli_likelihood("logit")),
            list(y = stats::rpois(30, exp(1 + sin(x))),
                 likelihood = poisson_likelihood())
          )
          errors <- vapply(cases, function(case) {
            kernel <- matern52_kernel(1.5, 1.2)
            theta <- validation_log_parameters(kernel)
            objective <- function(values) {
              names(values) <- names(theta)
              fit_latent_gp(x, case$y, validation_kernel_at(kernel, values),
                            case$likelihood)$log_marginal_likelihood
            }
            gradient <- log_marginal_likelihood_gradient(
              fit_latent_gp(x, case$y, kernel, case$likelihood)
            )
            validation_scaled_error(
              gradient[names(theta)],
              validation_numerical_gradient(objective, theta)
            )
          }, numeric(1))
          max(errors)
        })
      }
    ),
    validation_reference(
      "latent-predictive-integration",
      "Laplace approximation",
      "class probabilities (logit and probit links) and Poisson-lognormal probabilities, distribution function, and predictive mean (largest absolute difference)",
      "stats::integrate() of the likelihood against the Gaussian latent predictive density",
      "brute force",
      tolerance = 1e-7,
      run = function() {
        mean <- c(-2, 0.3, 1.5)
        variance <- c(0.5, 1, 4)
        # Twenty standard deviations either side; exp(z) overflows at
        # infinite limits.
        integral <- function(f, m, v) {
          stats::integrate(function(z) f(z) * stats::dnorm(z, m, sqrt(v)),
                           m - 20 * sqrt(v), m + 20 * sqrt(v),
                           rel.tol = 1e-12)$value
        }
        probability_errors <- vapply(c("logit", "probit"), function(link) {
          link_function <- if (link == "logit") stats::plogis else stats::pnorm
          ours <- likelihood_predictive(bernoulli_likelihood(link), mean,
                                        variance)$probability
          reference <- mapply(function(m, v) integral(link_function, m, v),
                              mean, variance)
          max(abs(ours - reference))
        }, numeric(1))
        counts <- 0:12
        ours_density <- exp(likelihood_predictive(
          poisson_likelihood(), rep(log(3), 13), rep(0.4, 13),
          observed = counts, exposure = rep(1.5, 13)
        )$log_density)
        reference_density <- vapply(counts, function(k) {
          integral(function(z) stats::dpois(k, 1.5 * exp(z)), log(3), 0.4)
        }, numeric(1))
        ours_mean <- likelihood_predictive(poisson_likelihood(), log(3), 0.4,
                                           exposure = 1.5)$mean
        reference_mean <- integral(function(z) 1.5 * exp(z), log(3), 0.4)
        distribution <- gaussianprocesses:::.count_distribution(
          counts, rep(log(3), 13), rep(0.4, 13), rep(1.5, 13), 40L
        )
        max(probability_errors,
            abs(ours_density - reference_density),
            abs(distribution - cumsum(reference_density)),
            abs(ours_mean - reference_mean) / reference_mean)
      }
    ),

    # Posterior and prior sampling ------------------------------------------
    validation_reference(
      "sampling-moments",
      "Posterior and prior sampling",
      "sample means, variances, and correlations of 20 000 posterior and prior draws at 8 inputs (largest deviation in Monte Carlo standard errors)",
      "the Gaussian posterior and prior moments from predict_gp() and evaluate_kernel()",
      "closed form",
      tolerance = 5,
      run = function() {
        data <- validation_data(30, seed = 25L)
        kernel <- matern52_kernel(1, 1.5)
        model <- fit_gp(data$x, data$y, kernel, noise_variance = 0.04)
        x_new <- c(-1, 1.5, 2, 4.4, 6, 8.5, 9, 11)
        posterior <- predict_gp(model, x_new, include_covariance = TRUE)
        n_draws <- 20000
        deviation <- function(draws, mean, covariance) {
          sample_covariance <- stats::cov(t(draws))
          max(
            abs(rowMeans(draws) - mean) / sqrt(diag(covariance) / n_draws),
            abs(diag(sample_covariance) / diag(covariance) - 1) /
              sqrt(2 / n_draws),
            abs(stats::cov2cor(sample_covariance) -
                  stats::cov2cor(covariance)) * sqrt(n_draws)
          )
        }
        max(
          deviation(sample_gp_posterior(model, x_new, n_draws, seed = 1)$latent,
                    posterior$mean, posterior$latent_covariance),
          deviation(sample_gp_prior(x_new, kernel, n_draws = n_draws,
                                    seed = 2)$latent,
                    rep(0, 8), evaluate_kernel(kernel, x_new))
        )
      }
    ),

    # Hyperparameter uncertainty --------------------------------------------
    validation_reference(
      "uncertainty-information",
      "Hyperparameter uncertainty",
      "observed information matrix of the log hyperparameters at the maximum",
      "second-order central differences of log_marginal_likelihood() (step 1e-3)",
      "brute force",
      tolerance = 1e-4,
      run = function() {
        data <- validation_data(50, seed = 26L)
        model <- optimize_gp(data$x, data$y, matern52_kernel(1, 1),
                             noise_variance = 0.1, n_starts = 2)
        uncertainty <- gp_hyperparameter_uncertainty(model)
        estimate <- log(model$optimization$optimized_parameters)
        likelihood <- function(values) {
          log_marginal_likelihood(fit_gp(
            data$x, data$y, matern52_kernel(exp(values[[1]]), exp(values[[2]])),
            noise_variance = exp(values[[3]])
          ))
        }
        h <- 1e-3
        hessian <- matrix(0, 3, 3)
        for (i in 1:3) {
          for (j in 1:3) {
            shift <- function(a, b) {
              likelihood(estimate + h * (a * (seq_len(3) == i) +
                                           b * (seq_len(3) == j)))
            }
            hessian[i, j] <- (shift(1, 1) - shift(1, -1) - shift(-1, 1) +
                                shift(-1, -1)) / (4 * h^2)
          }
        }
        validation_scaled_error(unname(uncertainty$information), -hessian)
      }
    ),
    validation_reference(
      "uncertainty-profile",
      "Hyperparameter uncertainty",
      "profile log likelihood of the length scale at 5 values (absolute difference)",
      "stats::optim() (BFGS, numerical gradients) over the other log parameters of the dense log density",
      "brute force",
      tolerance = 1e-6,
      run = function() {
        data <- validation_data(50, seed = 27L)
        model <- optimize_gp(data$x, data$y, matern52_kernel(1, 1),
                             noise_variance = 0.1, n_starts = 2)
        estimate <- model$optimization$optimized_parameters[["length_scale"]]
        values <- estimate * c(0.6, 0.8, 1, 1.25, 1.6)
        profile <- gp_profile_likelihood(model, "length_scale", values = values)
        reference <- vapply(values, function(length_scale) {
          objective <- function(log_values) {
            -validation_dense_posterior(
              matern52_kernel(exp(log_values[[1]]), length_scale),
              data$x, data$y, exp(log_values[[2]]), 0
            )$log_marginal_likelihood
          }
          start <- log(model$optimization$optimized_parameters[
            c("variance", "noise_variance")
          ])
          -stats::optim(start, objective, method = "BFGS",
                        control = list(reltol = 1e-14, maxit = 1000))$value
        }, numeric(1))
        max(abs(profile$log_marginal_likelihood - reference))
      }
    ),

    # Predictive scores -----------------------------------------------------
    validation_reference(
      "scores-gaussian",
      "Predictive scores",
      "pointwise CRPS and log scores of Gaussian predictions",
      "stats::integrate() of (F(z) - 1{z >= y})^2 over z, and stats::dnorm()",
      "brute force",
      tolerance = 1e-10,
      run = function() {
        observed <- c(-1.2, 0.1, 0.4, 2.5, 3)
        mean <- c(-1, 0.5, 0.4, 1.5, -0.5)
        sd <- c(0.3, 1, 0.2, 2, 1.5)
        scores <- gp_scores(observed, mean, sd)$pointwise
        crps <- mapply(function(y, m, s) {
          stats::integrate(function(z) {
            (stats::pnorm(z, m, s) - (z >= y))^2
          }, -Inf, y, rel.tol = 1e-12)$value +
            stats::integrate(function(z) {
              (stats::pnorm(z, m, s) - 1)^2
            }, y, Inf, rel.tol = 1e-12)$value
        }, observed, mean, sd)
        max(abs(scores$crps - crps),
            abs(scores$log_score - stats::dnorm(observed, mean, sd, log = TRUE)))
      }
    ),
    validation_reference(
      "scores-classification",
      "Predictive scores",
      "log loss, Brier score, and ROC AUC of 200 probability forecasts",
      "direct formulas and the AUC as the fraction of correctly ordered positive-negative pairs",
      "brute force",
      tolerance = 1e-10,
      run = function() {
        validation_with_seed(28L, function() {
          probability <- stats::runif(200)
          observed <- stats::rbinom(200, 1, probability)
          probability[1:10] <- round(probability[1:10], 1)
          scores <- gp_classification_scores(observed, probability)$summary
          pairs <- outer(probability[observed == 1],
                         probability[observed == 0], "-")
          auc <- mean((pairs > 0) + 0.5 * (pairs == 0))
          max(
            abs(scores[["log_loss"]] -
                  -mean(stats::dbinom(observed, 1, probability, log = TRUE))),
            abs(scores[["brier_score"]] - mean((probability - observed)^2)),
            abs(scores[["auc"]] - auc)
          )
        })
      }
    ),
    validation_reference(
      "scores-counts",
      "Predictive scores",
      "pointwise log scores and predictive means and variances of Poisson-lognormal count forecasts (largest absolute difference)",
      "stats::integrate() of the Poisson probabilities and moments against the latent predictive density",
      "brute force",
      tolerance = 1e-10,
      run = function() {
        observed <- c(0, 2, 5, 9, 30)
        latent_mean <- c(-0.5, 0.7, 1.4, 2, 2.5)
        latent_variance <- c(0.3, 0.5, 0.2, 0.8, 0.4)
        scores <- gp_count_scores(observed, latent_mean, latent_variance,
                                  seed = 1)$pointwise
        integral <- function(f, m, v) {
          stats::integrate(function(z) f(z) * stats::dnorm(z, m, sqrt(v)),
                           m - 20 * sqrt(v), m + 20 * sqrt(v),
                           rel.tol = 1e-12)$value
        }
        log_score <- mapply(function(k, m, v) {
          log(integral(function(z) stats::dpois(k, exp(z)), m, v))
        }, observed, latent_mean, latent_variance)
        # E[y] = E[exp(f)], Var[y] = E[exp(f)] + Var[exp(f)].
        predictive_mean <- mapply(function(m, v) {
          integral(exp, m, v)
        }, latent_mean, latent_variance)
        predictive_variance <- mapply(function(m, v) {
          integral(function(z) exp(2 * z), m, v)
        }, latent_mean, latent_variance) - predictive_mean^2 + predictive_mean
        max(
          abs(scores$log_score - log_score),
          abs(scores$predictive_mean - predictive_mean) / predictive_mean,
          abs(scores$predictive_variance - predictive_variance) /
            predictive_variance
        )
      }
    ),
    validation_reference(
      "rolling-origin-refit",
      "Predictive scores",
      "rolling-origin forecasts and forecast metrics (RMSE, coverage, log predictive density)",
      "expanding-window refits with fit_gp() and predict_gp(), and the metrics computed directly",
      "brute force",
      tolerance = 1e-10,
      run = function() {
        data <- validation_data(40, seed = 29L)
        kernel <- matern32_kernel(1, 1.5)
        backtest <- rolling_origin_gp(data$x, data$y, kernel,
                                      initial_window = 25, horizon = 3,
                                      step = 4, noise_variance = 0.04)
        reference <- do.call(rbind, lapply(c(25, 29, 33, 37), function(origin) {
          model <- fit_gp(data$x[seq_len(origin)], data$y[seq_len(origin)],
                          kernel, noise_variance = 0.04)
          target <- origin + 1:3
          prediction <- predict_gp(model, data$x[target])
          cbind(prediction$mean, prediction$observation_variance,
                data$y[target])
        }))
        metrics <- forecast_metrics(backtest)
        lower <- reference[, 1] - stats::qnorm(0.975) * sqrt(reference[, 2])
        upper <- reference[, 1] + stats::qnorm(0.975) * sqrt(reference[, 2])
        max(
          validation_scaled_error(backtest$mean, reference[, 1]),
          validation_scaled_error(backtest$observation_variance,
                                  reference[, 2]),
          abs(metrics$rmse - sqrt(mean((reference[, 3] - reference[, 1])^2))),
          abs(metrics$empirical_coverage -
                mean(reference[, 3] >= lower & reference[, 3] <= upper)),
          abs(metrics$mean_log_predictive_density -
                mean(stats::dnorm(reference[, 3], reference[, 1],
                                  sqrt(reference[, 2]), log = TRUE)))
        )
      }
    ),

    # Kernels ---------------------------------------------------------------
    validation_reference(
      "kernels-kernlab",
      "Kernels",
      "covariance matrices of the RBF, Matern-1/2 (exponential), and linear kernels on 20 three-dimensional inputs",
      "kernlab::kernelMatrix() with rbfdot, laplacedot, and vanilladot",
      "implementation",
      tolerance = 1e-8,
      package = "kernlab",
      recorded_version = "0.9-33",
      run = function() {
        x <- validation_with_seed(30L, function() {
          matrix(stats::runif(60, -2, 2), 20, 3)
        })
        reference <- function(kernel) {
          unclass(kernlab::kernelMatrix(kernel, x))[, ]
        }
        max(
          validation_scaled_error(
            evaluate_kernel(rbf_kernel(1, 0.8), x),
            reference(kernlab::rbfdot(sigma = 1 / (2 * 0.8^2)))
          ),
          validation_scaled_error(
            evaluate_kernel(matern12_kernel(1, 0.8), x),
            reference(kernlab::laplacedot(sigma = 1 / 0.8))
          ),
          validation_scaled_error(
            evaluate_kernel(linear_kernel(1), x),
            reference(kernlab::vanilladot())
          )
        )
      }
    ),
    validation_reference(
      "kernels-gradient",
      "Kernels",
      "derivatives of covariance matrices in every log parameter: rational quadratic x periodic, changepoint, and spectral-mixture kernels",
      "Richardson-extrapolated central differences of evaluate_kernel()",
      "brute force",
      tolerance = 1e-6,
      run = function() {
        # Irregular inputs: on an even grid every lag can be a multiple of
        # a component's half period, where its frequency derivative is zero.
        x <- validation_with_seed(31L, function() sort(stats::runif(12, -2, 3)))
        kernels <- list(
          product_kernel(rational_quadratic_kernel(1.2, 0.8, 1.5),
                         periodic_kernel(1, 1.3, 2)),
          changepoint_kernel(rbf_kernel(1, 1), matern52_kernel(0.5, 0.7),
                             location = 0.4, steepness = 2),
          spectral_mixture_kernel(2, weights = c(1, 0.5),
                                  frequencies = c(0.3, 1.1),
                                  spectral_variances = c(0.02, 0.05))
        )
        errors <- vapply(kernels, function(kernel) {
          values <- kernel_parameters(kernel, flatten = TRUE)
          # kernel_gradient() differentiates in each parameter's optimizer
          # coordinate: the log of a positive parameter, a real parameter
          # itself, or the logit of a bounded one.
          constraints <- gaussianprocesses:::.flatten_kernel_constraints(kernel)
          gradient <- kernel_gradient(kernel, x)
          max(vapply(names(gradient), function(name) {
            constraint <- constraints[[name]]
            at <- function(coordinate) {
              updated <- values
              updated[[name]] <- constraint$from_coordinate(coordinate)
              evaluate_kernel(update_kernel_parameters(kernel, updated), x)
            }
            coordinate <- constraint$to_coordinate(values[[name]])
            difference <- function(h) {
              (at(coordinate + h) - at(coordinate - h)) / (2 * h)
            }
            reference <- (4 * difference(5e-4) - difference(1e-3)) / 3
            validation_scaled_error(gradient[[name]], reference)
          }, numeric(1)))
        }, numeric(1))
        max(errors)
      }
    ),
    validation_reference(
      "kernels-spectral-density",
      "Kernels",
      "spectral-mixture covariance at 30 lags",
      "the Fourier integral of its spectral density, a mixture of Gaussians, with stats::integrate()",
      "brute force",
      tolerance = 1e-10,
      run = function() {
        weights <- c(1, 0.5, 0.2)
        frequencies <- c(0, 0.4, 1.3)
        spectral_variances <- c(0.01, 0.03, 0.1)
        kernel <- spectral_mixture_kernel(3, weights, frequencies,
                                          spectral_variances)
        lags <- seq(0, 6, length.out = 30)
        density <- function(s) {
          rowSums(vapply(1:3, function(q) {
            weights[[q]] / 2 * (
              stats::dnorm(s, frequencies[[q]], sqrt(spectral_variances[[q]])) +
                stats::dnorm(s, -frequencies[[q]], sqrt(spectral_variances[[q]]))
            )
          }, numeric(length(s))))
        }
        reference <- vapply(lags, function(lag) {
          stats::integrate(function(s) density(s) * cos(2 * pi * s * lag),
                           -Inf, Inf, rel.tol = 1e-12,
                           subdivisions = 1000L)$value
        }, numeric(1))
        validation_scaled_error(drop(evaluate_kernel(kernel, lags, 0)),
                                reference)
      }
    )
  )
}


# Runner ----------------------------------------------------------------------

# Runs references and returns one row per reference. A reference whose
# package is not installed is reported as "not run" without running it.
run_validation_references <- function(references = validation_references()) {
  rows <- lapply(references, function(reference) {
    available <- is.na(reference$package) ||
      requireNamespace(reference$package, quietly = TRUE)
    seconds <- NA_real_
    discrepancy <- NA_real_

    if (available) {
      seconds <- system.time(discrepancy <- reference$run())[["elapsed"]]
    }

    data.frame(
      id = reference$id,
      method = reference$method,
      quantity = reference$quantity,
      reference = reference$reference,
      type = reference$type,
      package = reference$package,
      recorded_version = reference$recorded_version,
      installed_version = if (!is.na(reference$package) && available) {
        as.character(utils::packageVersion(reference$package))
      } else {
        NA_character_
      },
      tolerance = reference$tolerance,
      discrepancy = discrepancy,
      status = if (!available) {
        "not run"
      } else if (discrepancy <= reference$tolerance) {
        "pass"
      } else {
        "fail"
      },
      seconds = seconds,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
