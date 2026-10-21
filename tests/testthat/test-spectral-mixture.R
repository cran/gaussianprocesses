# Spectral-mixture kernels.

spectral_inputs <- list(
  one = matrix(c(0, 0.3, 1, 1, 2.5, 4.2), ncol = 1),
  two = cbind(c(0, 0.3, 1, 1, 2.5, 4.2), c(1, -0.5, 0.2, 0.2, 2, 0.7))
)

# Mixtures with 1, 2, and 3 components; ARD values in two dimensions.
spectral_cases <- function(dimension) {
  per_dimension <- function(values) {
    if (dimension == 1L) values[, 1] else values
  }

  list(
    q1 = spectral_mixture_kernel(
      1,
      weights = 1.3,
      frequencies = per_dimension(rbind(c(0.2, 0.05))),
      spectral_variances = per_dimension(rbind(c(0.02, 0.04)))
    ),
    q2 = spectral_mixture_kernel(
      2,
      weights = c(1, 0.4),
      frequencies = per_dimension(rbind(c(0, 0.1), c(0.35, 0.2))),
      spectral_variances = per_dimension(rbind(c(0.01, 0.03), c(0.005, 0.02)))
    ),
    q3 = spectral_mixture_kernel(
      3,
      weights = c(0.8, 0.5, 0.3),
      frequencies = c(0.05, 0.25, 0.6),
      spectral_variances = c(0.02, 0.01, 0.05)
    )
  )
}

# The defining formula, with weights, frequencies, and spectral variances
# per component and dimension.
spectral_reference <- function(kernel, x, y) {
  components <- if (identical(kernel$type, "sum")) {
    kernel$children
  } else {
    list(kernel)
  }

  Reduce(
    `+`,
    lapply(
      components,
      function(component) {
        p <- component$parameters
        frequency <- rep_len(p$frequency, ncol(x))
        spectral_variance <- rep_len(p$spectral_variance, ncol(x))
        p$weight * Reduce(
          `*`,
          lapply(
            seq_len(ncol(x)),
            function(d) {
              tau <- outer(x[, d], y[, d], "-")
              exp(-2 * pi^2 * tau^2 * spectral_variance[[d]]) *
                cos(2 * pi * tau * frequency[[d]])
            }
          )
        )
      }
    )
  )
}


test_that("the kernel is the spectral-mixture formula", {
  for (dimension in 1:2) {
    x <- spectral_inputs[[dimension]]
    y <- x[c(2, 5), , drop = FALSE] + 0.15
    cases <- spectral_cases(dimension)

    for (name in names(cases)) {
      kernel <- cases[[name]]
      label <- paste(name, "D =", dimension)
      covariance <- evaluate_kernel(kernel, x)
      weights <- sum(kernel_parameters(kernel, flatten = TRUE)[
        grepl("weight", names(kernel_parameters(kernel, flatten = TRUE)))
      ])

      expect_equal(
        evaluate_kernel(kernel, x, y),
        spectral_reference(kernel, x, y),
        tolerance = 1e-14,
        info = label
      )
      expect_equal(diag(covariance), rep(weights, nrow(x)), info = label)
      expect_equal(kernel_diagonal(kernel, x), diag(covariance), info = label)
      # Rows 3 and 4 repeat, so the matrix is singular.
      expect_gte(
        min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values),
        -1e-10 * weights,
        label = paste(label, "smallest eigenvalue")
      )
    }
  }

  # cos is even, so a frequency and its negative give the same kernel.
  expect_equal(
    evaluate_kernel(spectral_mixture_kernel(1, 1, -0.3, 0.02),
                    spectral_inputs$one),
    evaluate_kernel(spectral_mixture_kernel(1, 1, 0.3, 0.02),
                    spectral_inputs$one),
    tolerance = 1e-15
  )
})


test_that("gradients and input derivatives match finite differences", {
  for (dimension in 1:2) {
    x <- spectral_inputs[[dimension]]
    y <- x[c(2, 5), , drop = FALSE] + 0.15
    cases <- spectral_cases(dimension)

    for (name in names(cases)) {
      kernel <- cases[[name]]
      label <- paste(name, "D =", dimension)
      analytical <- kernel_gradient(kernel, x, y)
      numerical <- coordinate_fd_kernel_gradient(kernel, x, y)

      expect_identical(names(analytical), names(numerical), info = label)

      for (path in names(analytical)) {
        expect_lte(
          max(abs(analytical[[path]] - numerical[[path]])) /
            max(1, max(abs(analytical[[path]]))),
          1e-7,
          label = paste(label, path)
        )
      }

      expect_input_derivative_contract(kernel, x, y, label)
      expect_input_parameter_gradient_contract(kernel, x, y, label)
    }
  }
})


test_that("a component with frequency 0 is an RBF kernel", {
  for (length_scale in list(0.7, c(0.4, 1.6))) {
    x <- spectral_inputs[[length(length_scale)]]
    spectral <- spectral_mixture_kernel(
      1,
      weights = 2.3,
      frequencies = 0,
      spectral_variances = rbind(1 / (4 * pi^2 * length_scale^2))
    )

    expect_lt(
      max(abs(
        evaluate_kernel(spectral, x) -
          evaluate_kernel(rbf_kernel(2.3, length_scale), x)
      )),
      1e-12
    )
  }
})


test_that("the constructor shapes and validates its parameters", {
  kernel <- spectral_mixture_kernel(2, weights = 0.5, frequencies = c(0, 0.3))

  expect_identical(kernel$type, "sum")
  expect_identical(
    names(kernel_parameters(kernel, flatten = TRUE)),
    c(
      "kernel1.weight", "kernel1.frequency", "kernel1.spectral_variance",
      "kernel2.weight", "kernel2.frequency", "kernel2.spectral_variance"
    )
  )
  expect_identical(kernel$children[[2]]$parameters$frequency, 0.3)
  expect_identical(kernel$children[[1]]$parameters$weight, 0.5)

  single <- spectral_mixture_kernel()
  expect_identical(single$type, "spectral_component")
  expect_identical(
    single$parameters,
    list(weight = 1, frequency = 0, spectral_variance = 1)
  )
  expect_identical(
    gaussianprocesses:::.flatten_kernel_constraints(single)$frequency$type,
    "real"
  )

  for (call in list(
    function() spectral_mixture_kernel(0),
    function() spectral_mixture_kernel(1.5),
    function() spectral_mixture_kernel(2, weights = c(1, 2, 3)),
    function() spectral_mixture_kernel(2, weights = c(1, -1)),
    function() spectral_mixture_kernel(2, frequencies = matrix(0, 3, 2)),
    function() spectral_mixture_kernel(2, spectral_variances = 0),
    function() spectral_mixture_kernel(1, frequencies = NA),
    function() spectral_mixture_kernel(1, weights = "a")
  )) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }

  expect_identical(
    utils::capture.output(print(single)),
    "SpectralComponent(weight=1, frequency=0, spectral_variance=1)"
  )
})


sinusoids <- function(x) {
  2 * sin(2 * pi * 0.1 * x) + sin(2 * pi * 0.27 * x + 1)
}


test_that("initialisation recovers the dominant frequencies", {
  for (sampling in c("regular", "irregular")) {
    set.seed(29)
    x <- if (sampling == "regular") 0:199 else sort(runif(200, 0, 200))
    y <- sinusoids(x) + rnorm(200, sd = 0.5)
    resolution <- 1 / diff(range(x))

    kernel <- initialize_spectral_mixture(x, y, n_components = 2)
    parameters <- lapply(kernel$children, `[[`, "parameters")
    weights <- vapply(parameters, `[[`, numeric(1), "weight")
    frequencies <- vapply(parameters, `[[`, numeric(1), "frequency")[
      order(weights, decreasing = TRUE)
    ]

    # Within one frequency-resolution bin, the larger amplitude first.
    expect_lt(abs(frequencies[[1]] - 0.1), resolution, label = sampling)
    expect_lt(abs(frequencies[[2]] - 0.27), resolution, label = sampling)
    expect_true(all(weights > 0))
  }

  expect_error(
    initialize_spectral_mixture(cbind(1:20, 1:20), rnorm(20)),
    "one input dimension",
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    initialize_spectral_mixture(1:5, rnorm(5), n_components = 2),
    "at least 6 distinct inputs",
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    initialize_spectral_mixture(1:20, rnorm(20), n_components = 0),
    class = "gaussianprocesses_validation_error"
  )
  expect_identical(
    initialize_spectral_mixture(0:39, sinusoids(0:39), n_components = 1)$type,
    "spectral_component"
  )
})


test_that("optimize_gp() bounds frequencies by the Nyquist frequency", {
  set.seed(3)
  x <- 0:59
  y <- sinusoids(x) + rnorm(60, sd = 0.3)
  kernel <- spectral_mixture_kernel(2, weights = 1,
                                    frequencies = c(0.12, 0.9),
                                    spectral_variances = 1e-4)

  # Start 3 displaces kernel1.frequency; a few iterations suffice to check
  # the bounds.
  model <- optimize_gp(x, y, kernel, noise_variance = 0.1, n_starts = 3,
                       control = list(maxit = 5))
  starts <- model$optimization$starts
  resolution <- 1 / 59

  # The initial 0.9 is clamped to the Nyquist frequency 0.5, and later
  # starts move a frequency by one resolution bin.
  expect_identical(unname(starts[1, "kernel2.frequency"]), 0.5)
  displaced <- starts[, "kernel1.frequency"] != 0.12
  expect_equal(
    unname(abs(starts[displaced, "kernel1.frequency"] - 0.12)),
    resolution
  )

  frequencies <- vapply(
    model$optimization$start_results,
    function(result) {
      result$parameters[c("kernel1.frequency", "kernel2.frequency")]
    },
    numeric(2)
  )
  expect_true(all(frequencies >= 0 & frequencies <= 0.5))

  # Irregular inputs: half the reciprocal of the smallest spacing.
  limits <- gaussianprocesses:::.spectral_frequency_limits(
    select_dimensions(spectral_mixture_kernel(1), 2),
    cbind(1:4, c(0, 0.5, 0.75, 2)),
    1:2
  )
  expect_identical(names(limits), "frequency")
  expect_equal(limits$frequency$upper, 2)
  expect_equal(limits$frequency$resolution, 0.5)
})


test_that("initialisation and optimization find the sinusoids (smoke test)", {
  # Slow; skipped on CRAN to keep its check within about a minute.
  skip_on_cran()
  set.seed(29)
  x <- 0:59
  y <- sinusoids(x) + rnorm(60, sd = 0.3)
  kernel <- initialize_spectral_mixture(x, y, n_components = 2)

  # One start from the initialisation; the displaced starts of the previous
  # test take a few hundred fits each.
  model <- optimize_gp(x, y, kernel, noise_variance = 0.1, n_starts = 1)
  frequencies <- sort(vapply(
    model$kernel$children,
    function(component) component$parameters$frequency,
    numeric(1)
  ))

  expect_lt(max(abs(frequencies - c(0.1, 0.27))), 0.01)
  expect_true(model$optimization$converged)
})
