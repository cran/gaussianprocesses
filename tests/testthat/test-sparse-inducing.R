# Inducing-point gradients and optimization, greedy variance selection, and
# collapse detection.

test_that("gradients in the inducing inputs pass finite-difference checks", {
  set.seed(36)

  for (n_features in 1:3) {
    for (n_inducing in c(3L, 6L, 10L)) {
      x <- matrix(stats::runif(30 * n_features, 0, 4), 30, n_features)
      y <- sin(rowSums(x)) + stats::rnorm(30, sd = 0.2)
      kernel <- if (n_features == 2L) {
        sum_kernel(rbf_kernel(1, c(1, 1.5)), linear_kernel(0.2))
      } else {
        matern52_kernel(1.2, 1.1)
      }
      inducing <- x[sample(30, n_inducing), , drop = FALSE] + 0.05

      for (method in c("fitc", "vfe")) {
        fit_at <- function(values) {
          fit_sparse_gp(x, y, kernel, noise_variance = 0.05,
                        inducing_points = matrix(values, n_inducing,
                                                 n_features),
                        method = method)
        }
        analytical <- .sparse_gradient(fit_at(as.vector(inducing)),
                                       include_noise = FALSE,
                                       include_inducing = TRUE)
        analytical <- analytical[startsWith(names(analytical), "inducing.")]

        # Richardson differences with a step of 1e-3: their O(h^4) error is
        # far below the rounding error of smaller steps.
        numerical <- vapply(seq_along(inducing), function(k) {
          shifted <- function(step) {
            up <- as.vector(inducing)
            down <- up
            up[k] <- up[k] + step
            down[k] <- down[k] - step
            (fit_at(up)$log_marginal_likelihood -
                fit_at(down)$log_marginal_likelihood) / (2 * step)
          }
          (4 * shifted(5e-4) - shifted(1e-3)) / 3
        }, numeric(1))

        expect_identical(
          names(analytical),
          .inducing_parameter_names(n_inducing, n_features)
        )
        expect_lt(
          max(abs(analytical - numerical) / pmax(abs(numerical), 1e-3)),
          1e-6,
          label = sprintf("%s, D = %d, m = %d", method, n_features,
                          n_inducing)
        )
      }
    }
  }
})


test_that("greedy variance selection is the pivoted Cholesky order", {
  set.seed(1)

  for (case in 1:4) {
    n_features <- if (case %% 2 == 0) 2L else 1L
    x <- matrix(stats::runif(30 * n_features), 30, n_features)
    kernel <- if (case <= 2) {
      rbf_kernel(1, rep(0.3, n_features))
    } else {
      matern32_kernel(0.7, rep(0.4, n_features))
    }
    selected <- select_inducing_points(x, 8, method = "variance",
                                       kernel = kernel)
    pivot <- attr(
      suppressWarnings(chol(evaluate_kernel(kernel, x), pivot = TRUE)),
      "pivot"
    )

    expect_identical(attr(selected, "indices"), pivot[1:8],
                     label = paste("case", case))
    expect_identical(attr(selected, "method"), "variance")
    expect_equal(matrix(selected, ncol = n_features),
                 x[pivot[1:8], , drop = FALSE])
  }

  # The kernel is evaluated one column at a time, n values per chosen point,
  # never as the full n x n matrix.
  evaluated <- 0
  original <- evaluate_kernel
  local_mocked_bindings(evaluate_kernel = function(kernel, x, y = NULL) {
    value <- original(kernel, x, y)
    evaluated <<- evaluated + length(value)
    value
  })
  select_inducing_points(stats::runif(500), 20, method = "variance",
                         kernel = rbf_kernel(length_scale = 0.05))
  expect_identical(evaluated, 500 * 20)
})


test_that("greedy selection stops when the prior variance is explained", {
  # The linear kernel on one input has rank one: one point explains it.
  expect_warning(
    selected <- select_inducing_points(seq(1, 2, length.out = 20), 3,
                                       method = "variance",
                                       kernel = linear_kernel()),
    class = "gaussianprocesses_inducing_warning"
  )
  expect_identical(nrow(selected), 1L)

  expect_error(
    select_inducing_points(1:10, 3, method = "variance"),
    class = "gaussianprocesses_validation_error"
  )

  # fit_sparse_gp() passes its kernel to the selection.
  set.seed(2)
  x <- sort(stats::runif(200, 0, 10))
  y <- sin(x) + stats::rnorm(200, sd = 0.1)
  kernel <- rbf_kernel(1, 1)
  variance <- fit_sparse_gp(x, y, kernel, noise_variance = 0.01,
                            n_inducing = 12, selection = "variance",
                            method = "vfe")
  farthest <- fit_sparse_gp(x, y, kernel, noise_variance = 0.01,
                            n_inducing = 12, method = "vfe")
  expect_identical(variance$selection_method, "variance")
  expect_lt(variance$trace_term, farthest$trace_term)
})


test_that("inducing inputs are optimized jointly with the hyperparameters", {
  set.seed(2)
  x <- sort(stats::runif(150, -4, 4))
  y <- sin(2 * x) + stats::rnorm(150, sd = sqrt(0.1))

  fixed <- optimize_sparse_gp(x, y, rbf_kernel(), noise_variance = 0.5,
                              n_inducing = 6, selection = "variance",
                              n_starts = 1)
  free <- optimize_sparse_gp(x, y, rbf_kernel(), noise_variance = 0.5,
                             n_inducing = 6, selection = "variance",
                             optimize_inducing = TRUE, n_starts = 1)
  exact <- log_marginal_likelihood(
    fit_gp(x, y, free$kernel, noise_variance = free$noise_variance)
  )

  expect_true(free$optimization$optimize_inducing)
  expect_false(fixed$optimization$optimize_inducing)
  expect_identical(free$selection_method, "optimized")
  expect_null(free$inducing_indices)
  expect_named(
    free$optimization$optimized_parameters,
    c("variance", "length_scale", "noise_variance",
      .inducing_parameter_names(6, 1))
  )
  expect_equal(
    unname(free$inducing_points[, 1]),
    unname(free$optimization$optimized_parameters[
      .inducing_parameter_names(6, 1)
    ])
  )
  # Moving the inducing inputs can only raise the bound, which stays below
  # the exact marginal likelihood.
  expect_gt(free$log_marginal_likelihood, fixed$log_marginal_likelihood)
  expect_lt(free$log_marginal_likelihood, exact)

  expect_error(
    optimize_sparse_gp(x, y, matern12_kernel(), noise_variance = 0.5,
                       n_inducing = 6, optimize_inducing = TRUE),
    class = "gaussianprocesses_smoothness_error"
  )
  expect_error(
    optimize_sparse_gp(x, y, rbf_kernel(), n_inducing = 6,
                       optimize_inducing = NA),
    class = "gaussianprocesses_validation_error"
  )
})


test_that("nearly identical inducing points are detected", {
  x <- seq(0, 5, length.out = 60)
  y <- sin(x)
  inducing <- c(0.5, 1.5, 1.5 + 1e-9, 3, 4.5)

  expect_warning(
    model <- fit_sparse_gp(x, y, rbf_kernel(), noise_variance = 0.01,
                           inducing_points = inducing, method = "vfe"),
    class = "gaussianprocesses_inducing_collapse_warning"
  )
  expect_identical(model$inducing_collapse$first, 2L)
  expect_identical(model$inducing_collapse$second, 3L)
  expect_gt(model$inducing_collapse$correlation, 1 - 1e-6)
  expect_gt(model$kuu_jitter, 0)
  expect_true(is.finite(model$log_marginal_likelihood))

  expect_no_warning(
    separated <- fit_sparse_gp(x, y, rbf_kernel(), noise_variance = 0.01,
                               inducing_points = c(0.5, 1.5, 3, 4.5))
  )
  expect_null(separated$inducing_collapse)
})
