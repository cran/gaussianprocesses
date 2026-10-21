# Changepoint kernels.

changepoint_x <- cbind(
  c(1, 3.5, 4.9, 5, 5, 5.2, 6.4, 9),
  c(0.3, -1, 0.8, 0.2, 0.2, -0.4, 1.1, 0)
)
changepoint_y <- changepoint_x[c(1, 4, 7), ] + 0.25

changepoint_cases <- function() {
  list(
    simple = changepoint_kernel(
      rbf_kernel(1, 2),
      rbf_kernel(0.3, 0.5),
      location = 5,
      steepness = 3
    ),
    ard_children_second_column = changepoint_kernel(
      matern52_kernel(1.2, c(1.5, 0.8)),
      rational_quadratic_kernel(0.6, c(0.7, 1.3), alpha = 1.5),
      location = 0.1,
      steepness = 2,
      column = 2
    ),
    nested = changepoint_kernel(
      changepoint_kernel(
        rbf_kernel(1, 1.5),
        periodic_kernel(0.5, 1, period = 2),
        location = 3,
        steepness = 4
      ),
      linear_kernel(0.2),
      location = 6,
      steepness = 1.5
    ),
    selected = sum_kernel(
      select_dimensions(
        changepoint_kernel(rbf_kernel(), matern32_kernel(0.4, 0.6),
                           location = 5.1, steepness = 6),
        1
      ),
      select_dimensions(rbf_kernel(0.2, 1), 2)
    )
  )
}

# The defining formula, with dense matrices.
changepoint_reference <- function(kernel, x, y) {
  rise <- function(points) {
    stats::plogis(
      kernel$parameters$steepness *
        (points[, kernel$column] - kernel$parameters$location)
    )
  }

  outer(1 - rise(x), 1 - rise(y)) *
    evaluate_kernel(kernel$children[[1]], x, y) +
    outer(rise(x), rise(y)) * evaluate_kernel(kernel$children[[2]], x, y)
}


test_that("the changepoint kernel is the weighted sum of its two kernels", {
  cases <- changepoint_cases()

  for (name in c("simple", "ard_children_second_column", "nested")) {
    kernel <- cases[[name]]

    expect_equal(
      evaluate_kernel(kernel, changepoint_x, changepoint_y),
      changepoint_reference(kernel, changepoint_x, changepoint_y),
      tolerance = 1e-14,
      info = name
    )
  }

  for (name in names(cases)) {
    kernel <- cases[[name]]
    covariance <- evaluate_kernel(kernel, changepoint_x)

    expect_equal(covariance, t(covariance), tolerance = 1e-14, info = name)
    expect_equal(
      kernel_diagonal(kernel, changepoint_x),
      diag(covariance),
      tolerance = 1e-14,
      info = name
    )
    # Rows 4 and 5 repeat, so the matrix is singular.
    expect_gte(
      min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values),
      -1e-10 * max(diag(covariance)),
      label = paste(name, "smallest eigenvalue")
    )
  }
})


test_that("hyperparameter gradients match finite differences", {
  for (case in names(changepoint_cases())) {
    kernel <- changepoint_cases()[[case]]

    for (y in list(NULL, changepoint_y)) {
      analytical <- kernel_gradient(kernel, changepoint_x, y)
      numerical <- coordinate_fd_kernel_gradient(kernel, changepoint_x, y)

      expect_identical(names(analytical), names(numerical), info = case)

      for (path in names(analytical)) {
        expect_lte(
          max(abs(analytical[[path]] - numerical[[path]])) /
            max(1, max(abs(analytical[[path]]))),
          1e-7,
          label = paste(case, path)
        )
      }
    }
  }
})


test_that("input derivatives and their parameter derivatives are correct", {
  for (case in names(changepoint_cases())) {
    kernel <- changepoint_cases()[[case]]

    expect_input_derivative_contract(kernel, changepoint_x, changepoint_y,
                                     case)
    expect_input_parameter_gradient_contract(kernel, changepoint_x,
                                             changepoint_y, case)
  }

  expect_error(
    kernel_input_gradient(
      changepoint_kernel(rbf_kernel(), matern12_kernel(), location = 5),
      changepoint_x
    ),
    class = "gaussianprocesses_smoothness_error"
  )

  # Derivative observations of a changepoint process.
  kernel <- changepoint_cases()$simple
  x <- c(1, 3, 4.5, 5.5, 7, 9, 2, 8)
  types <- c(0, 0, 0, 0, 0, 0, 1, 1)
  y <- c(sin(x[1:6]), cos(x[7:8]))
  model <- fit_gp(x, y, kernel, noise_variance = 0.01, derivative = types)
  analytical <- log_marginal_likelihood_gradient(model)
  step <- 1e-6
  numerical <- vapply(
    c("location", "steepness"),
    function(path) {
      at <- function(direction) {
        shifted <- kernel_parameters(kernel, flatten = TRUE)[path]
        shifted[] <- if (path == "location") {
          shifted + direction * step
        } else {
          shifted * exp(direction * step)
        }
        log_marginal_likelihood(
          fit_gp(x, y, update_kernel_parameters(kernel, shifted),
                 noise_variance = 0.01, derivative = types)
        )
      }
      (at(1) - at(-1)) / (2 * step)
    },
    numeric(1)
  )
  expect_equal(analytical[c("location", "steepness")], numerical,
               tolerance = 1e-6)
})


test_that("a steep change separates the two regimes", {
  before <- rbf_kernel(1, 2)
  after <- matern52_kernel(0.4, 0.7)
  kernel <- changepoint_kernel(before, after, location = 5, steepness = 500)
  # Every point is at least 0.1 from the change, where the sigmoid is within
  # exp(-50) of 0 or 1.
  x <- c(1, 2.5, 4.9, 5.1, 6, 8.5)
  left <- x < 5
  covariance <- evaluate_kernel(kernel, x)
  tolerance <- 1e-12

  expect_lt(max(abs(covariance[left, !left])), tolerance)
  expect_lt(
    max(abs(covariance[left, left] - evaluate_kernel(before, x[left]))),
    tolerance
  )
  expect_lt(
    max(abs(covariance[!left, !left] - evaluate_kernel(after, x[!left]))),
    tolerance
  )

  # A gentle change mixes them.
  gentle <- evaluate_kernel(
    changepoint_kernel(before, after, location = 5, steepness = 0.5),
    x
  )
  expect_gt(max(abs(gentle[left, !left])), 0.1)
})


test_that("parameter paths, updates, and printing", {
  kernel <- changepoint_cases()$nested

  expect_identical(
    names(kernel_parameters(kernel, flatten = TRUE)),
    c(
      "location", "steepness",
      "before.location", "before.steepness",
      "before.before.variance", "before.before.length_scale",
      "before.after.variance", "before.after.length_scale",
      "before.after.period",
      "after.variance"
    )
  )
  expect_identical(
    names(kernel_parameters(kernel)),
    c("location", "steepness", "before", "after")
  )
  expect_identical(
    gaussianprocesses:::.flatten_kernel_constraints(kernel)$location$type,
    "real"
  )
  expect_identical(
    gaussianprocesses:::.flatten_kernel_constraints(kernel)$steepness$type,
    "positive"
  )

  updated <- update_kernel_parameters(
    kernel,
    c(location = -2, "before.steepness" = 7, "after.variance" = 0.5)
  )
  expect_identical(updated$parameters$location, -2)
  expect_identical(updated$children[[1]]$parameters$steepness, 7)
  expect_identical(updated$children[[2]]$parameters$variance, 0.5)
  expect_error(
    update_kernel_parameters(kernel, c(steepness = -1)),
    class = "gaussianprocesses_validation_error"
  )

  expect_identical(
    utils::capture.output(print(changepoint_cases()$simple)),
    c(
      "ChangepointKernel(location=5, steepness=3, column=1,",
      "  RBF(variance=1, length_scale=2)",
      "  RBF(variance=0.3, length_scale=0.5)",
      ")"
    )
  )
})


test_that("invalid changepoints raise classed errors", {
  for (call in list(
    function() changepoint_kernel(list(), rbf_kernel()),
    function() changepoint_kernel(rbf_kernel(), "a"),
    function() changepoint_kernel(rbf_kernel(), rbf_kernel(), location = NA),
    function() changepoint_kernel(rbf_kernel(), rbf_kernel(), steepness = 0),
    function() changepoint_kernel(rbf_kernel(), rbf_kernel(), column = 0),
    function() changepoint_kernel(rbf_kernel(), rbf_kernel(), column = 1.5),
    function() changepoint_kernel(rbf_kernel(), rbf_kernel(), column = 1:2)
  )) {
    expect_error(call(), class = "gaussianprocesses_validation_error")
  }

  outside <- changepoint_kernel(rbf_kernel(), rbf_kernel(), column = 3)
  expect_error(
    evaluate_kernel(outside, changepoint_x),
    "acts along column 3, but its inputs have 2 column(s)",
    fixed = TRUE,
    class = "gaussianprocesses_validation_error"
  )
  expect_error(
    fit_gp(changepoint_x, changepoint_x[, 1], outside),
    "acts along column 3",
    class = "gaussianprocesses_validation_error"
  )
})


test_that("optimize_gp() bounds locations and starts them across quantiles", {
  x <- cbind(seq(0, 10, length.out = 40), rep(c(-1, 1), 20))
  y <- sin(x[, 1]) + 0.1 * x[, 2]
  kernel <- sum_kernel(
    select_dimensions(
      changepoint_kernel(rbf_kernel(), rbf_kernel(), location = 50),
      1
    ),
    select_dimensions(linear_kernel(), 2)
  )

  model <- optimize_gp(x, y, kernel, noise_variance = 0.05, n_starts = 4,
                       control = list(maxit = 5))
  starts <- model$optimization$starts[, "kernel1.location"]

  # The first start is clamped into the range of column 1 and the others are
  # its quartiles.
  expect_identical(unname(starts), c(10, 2.5, 5, 7.5))
  locations <- vapply(
    model$optimization$start_results,
    function(result) result$parameters[["kernel1.location"]],
    numeric(1)
  )
  expect_true(all(locations >= 0 & locations <= 10))
  expect_identical(
    unname(model$optimization$coordinates[["kernel1.location"]]),
    "identity"
  )
})


test_that("optimization recovers a simulated changepoint (smoke test)", {
  # Slow; skipped on CRAN to keep its check within about a minute.
  skip_on_cran()
  # A smoke test on one seeded data set, not a coverage claim: over seeds 1 to
  # 8 the estimates were within 0.3 of the true location.
  truth <- changepoint_kernel(
    rbf_kernel(1, 2),
    rbf_kernel(1, 0.3),
    location = 6,
    steepness = 5
  )
  x <- seq(0, 10, length.out = 80)
  simulated <- simulate_gp_data(x, truth, noise_variance = 0.01, seed = 1)

  model <- optimize_gp(
    x,
    simulated$observed,
    changepoint_kernel(rbf_kernel(), rbf_kernel(), location = 2,
                       steepness = 2),
    noise_variance = 0.05,
    n_starts = 4
  )

  expect_lt(abs(model$kernel$parameters$location - 6), 0.5)
  expect_true(model$optimization$converged)
})


test_that("time-series kernels and models accept a changepoint", {
  kernel <- time_series_kernel(period = 12, changepoint = 30,
                               changepoint_steepness = 2)

  expect_identical(kernel$type, "changepoint")
  expect_identical(kernel$parameters, list(location = 30, steepness = 2))
  expect_identical(kernel$children[[1]], time_series_kernel(period = 12))
  expect_identical(kernel$children[[2]], time_series_kernel(period = 12))
  expect_identical(
    names(kernel_parameters(kernel, flatten = TRUE))[1:4],
    c("location", "steepness", "before.kernel1.variance",
      "before.kernel2.variance")
  )

  time <- as.Date("2026-01-01") + 0:47
  break_index <- gp_time_index(c(time[1], as.Date("2026-01-31")))[2]
  expect_identical(break_index, 30)

  y <- ifelse(seq_along(time) <= 30, sin(seq_along(time) / 3), 2)
  model <- fit_time_series_gp(
    time,
    y,
    time_series_kernel(period = 12, changepoint = break_index),
    noise_variance = 0.05
  )
  forecast <- forecast_gp(model, max(time) + 1:3)
  expect_true(all(forecast$region == "forecast"))
  expect_length(forecast$mean, 3L)
})
