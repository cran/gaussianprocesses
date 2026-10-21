# Gauss-Hermite quadrature.

test_that("Gauss-Hermite rules match published tables", {
  # Abramowitz and Stegun (1964), formula 25.4.46 and table 25.10, n = 5.
  rule <- .gauss_hermite(5)

  expect_lt(
    max(abs(rule$nodes - c(-2.020182870456086, -0.958572464613819, 0,
                           0.958572464613819, 2.020182870456086))),
    1e-14
  )
  expect_lt(
    max(abs(rule$weights - c(0.0199532420590459, 0.393619323152241,
                             0.945308720482942, 0.393619323152241,
                             0.0199532420590459))),
    1e-14
  )
  expect_equal(.gauss_hermite(1), list(nodes = 0, weights = sqrt(pi)))
  expect_identical(.gauss_hermite(5L), rule)
  expect_error(.gauss_hermite(0), class = "gaussianprocesses_validation_error")
})


test_that("Gauss-Hermite rules integrate polynomials of degree 2n - 1", {
  for (order in c(2, 7, 20, 40)) {
    rule <- .gauss_hermite(order)

    # int t^k exp(-t^2) dt is Gamma((k + 1) / 2) for even k and 0 otherwise,
    # compared relative to sum_j w_j |t_j|^k, the scale of rounding in the
    # sum.
    for (degree in seq(0, 2 * order - 1)) {
      moment <- sum(rule$weights * rule$nodes^degree)
      exact <- if (degree %% 2 == 0) gamma((degree + 1) / 2) else 0

      expect_lt(
        abs(moment - exact),
        1e-12 * sum(rule$weights * abs(rule$nodes)^degree),
        label = sprintf("order %d, degree %d", order, degree)
      )
    }

    # The smallest weights keep their relative accuracy: the recurrence for
    # the weights at the outermost node, w = 1 / sum_k p_k(t)^2, is the
    # reciprocal of the Christoffel function; check it against the closed
    # form w = 2^(n-1) n! sqrt(pi) / (n^2 H_(n-1)(t)^2) in logs.
    outermost <- rule$nodes[order]
    hermite <- c(1, 2 * outermost)

    for (k in seq_len(max(order - 2, 0))) {
      hermite <- c(
        hermite,
        2 * outermost * hermite[k + 1] - 2 * k * hermite[k]
      )
    }

    log_weight <- (order - 1) * log(2) + lfactorial(order) + 0.5 * log(pi) -
      2 * log(order) - 2 * log(abs(hermite[order]))
    expect_equal(log(rule$weights[order]), log_weight, tolerance = 1e-12)
  }
})


test_that("quadrature reproduces closed-form Gaussian expectations", {
  mean <- c(-2, -0.5, 0, 0.7, 1.5, 3)

  for (order in c(20, 30, 40)) {
    for (variance in c(0.01, 0.25, 1)) {
      variance <- rep(variance, length(mean))

      # Probit: E[Phi(f)] = Phi(mu / sqrt(1 + sigma^2)).
      expect_lt(
        max(abs(
          .gaussian_expectation(function(f, index) stats::pnorm(f), mean,
                                variance, order) -
            stats::pnorm(mean / sqrt(1 + variance))
        )),
        1e-10
      )

      # Lognormal: E[exp(k f)] = exp(k mu + k^2 sigma^2 / 2).
      for (power in 1:3) {
        expect_lt(
          max(abs(
            .gaussian_expectation(function(f, index) exp(power * f), mean,
                                  variance, order) /
              exp(power * mean + power^2 * variance / 2) - 1
          )),
          1e-10
        )
      }
    }
  }

  # A wider latent distribution needs the upper end of the range.
  expect_lt(
    max(abs(
      .gaussian_expectation(function(f, index) stats::pnorm(f), mean,
                            rep(2.25, 6), 40) -
        stats::pnorm(mean / sqrt(3.25))
    )),
    1e-10
  )

  # The integrand sees which element each row belongs to; zero variance
  # evaluates it at the mean.
  expect_equal(
    .gaussian_expectation(function(f, index) f * index, c(1, 2, 3),
                          c(0.5, 0, 2), 10),
    c(1, 4, 9)
  )
  expect_equal(
    .gaussian_log_expectation(function(f, index) -f^2, c(0.3, -1), c(0, 0),
                              5),
    c(-0.09, -1)
  )
})


test_that("log expectations are accurate, adaptive or not", {
  mean <- c(-1, 0.5, 2)
  variance <- c(0.3, 1, 2)

  # log E[exp(2 f)] = 2 mu + 2 sigma^2, beyond the range of exp().
  extreme <- c(-500, 0.5, 400)
  expect_equal(
    .gaussian_log_expectation(function(f, index) 2 * f, extreme, variance, 30),
    2 * extreme + 2 * variance,
    tolerance = 1e-12
  )

  # With h quadratic, the adaptive rule at the mode is exact at any order:
  # h(f) = -(f - c)^2 / (2 s), so h + log N has curvature -(1 / s + 1 / v)
  # and its mode is the precision-weighted mean.
  centre <- c(3, -2, 0.5)
  spread <- c(0.1, 2, 0.5)
  h <- function(f, index) -(f - centre[index])^2 / (2 * spread[index])
  precision <- 1 / spread + 1 / variance
  mode <- (centre / spread + mean / variance) / precision
  exact <- 0.5 * log(spread / (spread + variance)) -
    (centre - mean)^2 / (2 * (spread + variance))

  for (order in c(1, 3, 10)) {
    expect_equal(
      .gaussian_log_expectation(h, mean, variance, order, center = mode,
                                scale = 1 / sqrt(precision)),
      exact,
      tolerance = 1e-12
    )
  }
})
