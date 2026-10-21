# How should the log marginal likelihood gradient form its trace term?
#
# Run manually after installing the package, for example
#   Rscript inst/benchmarks/gradient-numerics.R
# This file is intentionally excluded from CI and needs the callr package for
# the memory measurements. Run it on an otherwise idle machine.
#
# The gradient with respect to a log parameter eta is
#   (1/2) (alpha' D alpha - tr(C^{-1} D)),   D = dC / d eta.
# Two ways to form the trace are compared:
#   chol2inv  C^{-1} from the Cholesky factor once (about n^3 flops), then
#             sum(C^{-1} * D) per parameter (n^2 flops each). This is what
#             log_marginal_likelihood_gradient() does.
#   solves    two triangular solves per parameter, tr(R^{-1} (R^{-T} D)),
#             which never forms C^{-1} (about 2 n^3 flops per parameter).
#
# Accuracy needs a reference that is itself accurate when C is
# ill-conditioned. Finite differences of the log marginal likelihood lose
# accuracy as fast as the methods being tested, so the main reference is
# exact: a periodic kernel on a regular grid over one period gives a
# circulant C whose eigenvalues are sums of modified Bessel functions, and
# for y made of a few Fourier modes the gradient has a closed form with
# relative accuracy near machine precision. Richardson-extrapolated central
# differences are also reported, for kernels without a closed form.

library(gaussianprocesses)

kernel_derivatives <- function(model) {
  gaussianprocesses:::.kernel_value_and_gradient(
    model$kernel,
    model$x,
    model$x
  )$gradient
}

gradient_by_chol2inv <- function(model) {
  log_marginal_likelihood_gradient(model)
}

gradient_by_solves <- function(model) {
  factor <- model$cholesky
  alpha <- model$alpha
  trace_of <- function(derivative) {
    sum(diag(backsolve(factor, backsolve(factor, derivative, transpose = TRUE))))
  }
  gradient_term <- function(derivative) {
    0.5 * (drop(crossprod(alpha, derivative %*% alpha)) - trace_of(derivative))
  }

  derivatives <- kernel_derivatives(model)
  derivatives$noise_variance <- diag(model$noise_variance, model$n_observations)
  vapply(derivatives, gradient_term, numeric(1))
}


# Exact reference: periodic kernel on a regular grid ------------------------
#
# k(d) = s exp(-2 sin^2(pi d) / l^2) = s sum_k e^{-x} I_k(x) e^{2 pi i k d}
# with x = 1 / l^2, so on the grid t_a = a / n the eigenvalue of the
# circulant kernel matrix for Fourier mode m is
#   lambda_m = s n sum_{k = m mod n} e^{-x} I_k(x).
# The derivative matrices share the eigenvectors, with eigenvalues lambda_m
# (log variance), sigma^2 (log noise variance), and
#   s n sum_{k = m mod n} (-2 x) e^{-x} ((I_{k-1} + I_{k+1}) / 2 - I_k)
# (log length scale).

circulant_case <- function(n, variance, length_scale, target_condition) {
  x <- 1 / length_scale^2
  # High orders underflow (below about 1e-300) and R warns that precision is
  # lost; such terms contribute nothing next to the noise variance.
  scaled_bessel <- function(order) {
    suppressWarnings(besselI(x, abs(order), expon.scaled = TRUE))
  }
  aliases <- function(m) m + n * (-3:3)
  modes <- 0:(n - 1)

  lambda <- vapply(
    modes,
    function(m) variance * n * sum(scaled_bessel(aliases(m))),
    numeric(1)
  )
  length_scale_eigen <- vapply(
    modes,
    function(m) {
      k <- aliases(m)
      variance * n * sum(
        -2 * x * ((scaled_bessel(k - 1) + scaled_bessel(k + 1)) / 2 -
          scaled_bessel(k))
      )
    },
    numeric(1)
  )
  noise <- max(lambda) / target_condition

  # y from a few Fourier modes with fixed coefficients.
  grid <- (0:(n - 1)) / n
  used <- c(1, 2, 3, 5, 8)
  cosine <- c(0.9, -0.4, 0.7, 0.3, -0.2)
  sine <- c(-0.5, 0.8, 0.1, -0.6, 0.4)
  y <- rowSums(
    vapply(
      seq_along(used),
      function(i) {
        cosine[i] * cos(2 * pi * used[i] * grid) +
          sine[i] * sin(2 * pi * used[i] * grid)
      },
      numeric(n)
    )
  )

  exact_gradient <- function(derivative_eigen) {
    eigen_used <- derivative_eigen[used + 1]
    covariance_used <- lambda[used + 1] + noise
    quadratic <- sum((cosine^2 + sine^2) * (n / 2) * eigen_used / covariance_used^2)
    trace <- sum(derivative_eigen / (lambda + noise))
    c(gradient = 0.5 * (quadratic - trace), scale = 0.5 * (abs(quadratic) + abs(trace)))
  }

  reference <- rbind(
    variance = exact_gradient(lambda),
    length_scale = exact_gradient(length_scale_eigen),
    noise_variance = exact_gradient(rep(noise, n))
  )

  model <- fit_gp(
    grid,
    y,
    periodic_kernel(variance, length_scale, period = 1),
    noise_variance = noise
  )

  list(
    model = model,
    reference = reference,
    condition = (max(lambda) + noise) / (min(lambda) + noise)
  )
}


accuracy_rows <- list()

for (target in c(1e4, 1e8, 1e12)) {
  case <- circulant_case(
    n = 400,
    variance = 1.3,
    length_scale = 0.7,
    target_condition = target
  )
  stopifnot(case$model$jitter == 0)
  parameters <- rownames(case$reference)

  for (method in c("chol2inv", "solves")) {
    gradient <- switch(
      method,
      chol2inv = gradient_by_chol2inv(case$model),
      solves = gradient_by_solves(case$model)
    )[parameters]
    error <- abs(gradient - case$reference[, "gradient"])

    accuracy_rows[[length(accuracy_rows) + 1L]] <- data.frame(
      reference = "exact (circulant)",
      condition = signif(case$condition, 2),
      method = method,
      max_relative_error = signif(max(error / abs(case$reference[, "gradient"])), 2),
      max_error_over_term_scale = signif(max(error / case$reference[, "scale"]), 2)
    )
  }
}


# Richardson reference for kernels without a closed form --------------------

richardson_gradient <- function(model, step = 1e-2) {
  parameters <- c(
    log(kernel_parameters(model$kernel, flatten = TRUE)),
    noise_variance = log(model$noise_variance)
  )
  lml_at <- function(theta) {
    natural <- exp(theta)
    kernel <- update_kernel_parameters(model$kernel, natural[names(natural) != "noise_variance"])
    log_marginal_likelihood(
      fit_gp(model$x, model$y, kernel, noise_variance = natural[["noise_variance"]])
    )
  }
  central <- function(index, h) {
    shift <- replace(numeric(length(parameters)), index, h)
    (lml_at(parameters + shift) - lml_at(parameters - shift)) / (2 * h)
  }

  estimates <- vapply(
    seq_along(parameters),
    function(index) {
      coarse <- central(index, step)
      fine <- central(index, step / 2)
      c(value = (4 * fine - coarse) / 3, uncertainty = abs((4 * fine - coarse) / 3 - fine))
    },
    numeric(2)
  )
  colnames(estimates) <- names(parameters)
  estimates
}

ard_x <- cbind(seq(0, 1, length.out = 150), sin(1:150), cos(3 * (1:150)))
ard_y <- sin(4 * ard_x[, 1]) + 0.5 * ard_x[, 2]

for (noise in c(1e-2, 1e-6)) {
  model <- fit_gp(
    ard_x,
    ard_y,
    sum_kernel(rbf_kernel(1, c(0.3, 1.2, 2)), periodic_kernel(0.4, 0.8, period = 1.5)),
    noise_variance = noise
  )
  stopifnot(model$jitter == 0)
  reference <- richardson_gradient(model)
  condition <- 1 / rcond(gaussianprocesses:::.model_observation_covariance(model))

  for (method in c("chol2inv", "solves")) {
    gradient <- switch(
      method,
      chol2inv = gradient_by_chol2inv(model),
      solves = gradient_by_solves(model)
    )[colnames(reference)]

    accuracy_rows[[length(accuracy_rows) + 1L]] <- data.frame(
      reference = sprintf(
        "Richardson (uncertainty %.0e)",
        max(reference["uncertainty", ] / pmax(1, abs(reference["value", ])))
      ),
      condition = signif(condition, 2),
      method = method,
      max_relative_error = signif(
        max(abs(gradient - reference["value", ]) / pmax(1, abs(reference["value", ]))),
        2
      ),
      max_error_over_term_scale = NA_real_
    )
  }
}

cat("Accuracy\n")
print(do.call(rbind, accuracy_rows), row.names = FALSE)


# Cost ------------------------------------------------------------------------
#
# Time per gradient for ARD RBF kernels in d input dimensions, so p = d + 2
# parameters including the noise variance. The kernel derivative matrices,
# which both methods need, are included. Peak memory is measured in a fresh
# R process per case, as R's allocation high-water mark from gc().

# Runs in a fresh R process, so it defines everything it uses.
cost_case <- function(n, dimensions, method, repeats) {
  library(gaussianprocesses)

  gradient_by_solves <- function(model) {
    factor <- model$cholesky
    alpha <- model$alpha
    trace_of <- function(derivative) {
      sum(diag(backsolve(factor, backsolve(factor, derivative, transpose = TRUE))))
    }
    derivatives <- gaussianprocesses:::.kernel_value_and_gradient(
      model$kernel,
      model$x,
      model$x
    )$gradient
    derivatives$noise_variance <- diag(model$noise_variance, model$n_observations)
    vapply(
      derivatives,
      function(derivative) {
        0.5 * (drop(crossprod(alpha, derivative %*% alpha)) - trace_of(derivative))
      },
      numeric(1)
    )
  }

  set.seed(1)
  x <- matrix(stats::runif(n * dimensions), n, dimensions)
  model <- fit_gp(
    x,
    sin(rowSums(x)),
    rbf_kernel(1, rep(0.7, dimensions)),
    noise_variance = 0.01
  )
  gradient <- switch(
    method,
    chol2inv = log_marginal_likelihood_gradient,
    solves = gradient_by_solves
  )

  seconds <- vapply(
    seq_len(repeats),
    function(i) system.time(gradient(model))[["elapsed"]],
    numeric(1)
  )
  baseline <- sum(gc(reset = TRUE)[, 2L])
  gradient(model)
  peak_mb <- sum(gc()[, 6L]) - baseline

  list(seconds = stats::median(seconds), peak_mb = peak_mb)
}

cost_rows <- list()

for (n in c(250L, 500L, 1000L)) {
  for (dimensions in c(1L, 3L, 7L)) {
    for (method in c("chol2inv", "solves")) {
      result <- callr::r(
        cost_case,
        list(n = n, dimensions = dimensions, method = method, repeats = 3L)
      )

      cost_rows[[length(cost_rows) + 1L]] <- data.frame(
        n = n,
        parameters = dimensions + 2L,
        method = method,
        seconds = round(result$seconds, 3),
        peak_mb = round(result$peak_mb, 1)
      )
    }
  }
}

cat("\nCost per gradient\n")
print(do.call(rbind, cost_rows), row.names = FALSE)
