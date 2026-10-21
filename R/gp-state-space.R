# State-space inference for one-dimensional Matern kernels ------------------
#
# A Matern process with nu = p + 1/2 on one input is the first component of
# the state x = (f, f', ..., f^(p)) of the linear stochastic differential
# equation dx = F x dt + L dW (Hartikainen and Sarkka, 2010). F is the
# companion matrix of (s + lambda)^(p + 1) with lambda = sqrt(2 nu) / ell,
# L = (0, ..., 0, 1)', and W is white noise of spectral density q. The
# stationary covariance P solves the Lyapunov equation
#   F P + P F' + L q L' = 0,
# and k(tau) = H exp(F tau) P H' for tau >= 0 with H = (1, 0, ..., 0).
# Because N = F + lambda I is nilpotent,
#   exp(F t) = exp(-lambda t) (I + t N + t^2 N^2 / 2 + ...),
# a closed form with p + 1 terms. Sums of kernels stack their states, and
# scaled kernels scale P and q.
#
# Over sorted inputs t_1 <= ... <= t_n the process is the linear Gaussian
# state-space model
#   x_k = A_k x_(k-1) + q_k,  A_k = exp(F (t_k - t_(k-1))),  q_k ~ N(0, P - A_k P A_k'),
#   y_k = H x_k + e_k,        e_k ~ N(0, sigma_k^2),
# started from N(0, P). Kalman filtering and Rauch-Tung-Striebel smoothing
# give the exact Gaussian-process posterior and log marginal likelihood in
# O(n) time; tied inputs have A = I and Q = 0.


# State-space representation of a Matern leaf of order p (nu = p + 1/2).
.matern_state_space <- function(order) {
  force(order)

  function(parameters) {
    variance <- parameters$variance
    length_scale <- parameters$length_scale

    if (length(length_scale) != 1L) {
      .gp_abort(
        "State-space inference needs a single length scale (one input).",
        c("gaussianprocesses_state_space_error",
          "gaussianprocesses_validation_error")
      )
    }

    dimension <- order + 1L
    lambda <- sqrt(2 * order + 1) / length_scale
    feedback <- matrix(0, dimension, dimension)

    if (dimension > 1L) {
      feedback[cbind(seq_len(dimension - 1L), 2:dimension)] <- 1
    }

    feedback[dimension, ] <- -choose(dimension, 0:(dimension - 1L)) *
      lambda^(dimension:1)

    stationary <- switch(
      dimension,
      matrix(variance, 1L, 1L),
      diag(c(variance, lambda^2 * variance)),
      {
        kappa <- lambda^2 * variance / 3
        matrix(
          c(variance, 0, -kappa, 0, kappa, 0, -kappa, 0, lambda^4 * variance),
          3L,
          3L
        )
      }
    )

    list(
      dimension = dimension,
      decay = lambda,
      feedback = feedback,
      stationary = stationary,
      spectral_density = variance *
        c(2 * lambda, 4 * lambda^3, 16 * lambda^5 / 3)[[dimension]],
      noise_effect = c(rep(0, dimension - 1L), 1),
      observation = c(1, rep(0, dimension - 1L))
    )
  }
}


.state_space_error <- function(message) {
  .gp_abort(
    message,
    c("gaussianprocesses_state_space_error",
      "gaussianprocesses_validation_error")
  )
}


# Leaf blocks of a kernel's state-space representation, with the scale of
# enclosing scale_kernel()s applied.
.state_space_blocks <- function(kernel, scale = 1) {
  if (length(kernel$children) == 0L) {
    entry <- .kernel_entry(kernel$type)

    if (is.null(entry$state_space)) {
      .state_space_error(sprintf(
        paste0(
          "The %s kernel has no exact state-space representation; use ",
          "method = \"exact\". Matern-1/2, 3/2, and 5/2 kernels, their ",
          "sums, and scaled versions have one."
        ),
        .kernel_display_name(kernel$type)
      ))
    }

    block <- entry$state_space(kernel$parameters)
    block$stationary <- scale * block$stationary
    block$spectral_density <- scale * block$spectral_density

    return(list(block))
  }

  switch(
    kernel$type,
    sum = do.call(c, lapply(kernel$children, .state_space_blocks,
                            scale = scale)),
    scale = .state_space_blocks(kernel$children[[1L]],
                                scale * kernel$parameters$scale),
    select = {
      if (!identical(as.integer(kernel$columns), 1L)) {
        .state_space_error(
          "State-space inference has one input; select_dimensions() must select column 1."
        )
      }

      .state_space_blocks(kernel$children[[1L]], scale)
    },
    .state_space_error(sprintf(
      paste0(
        "A %s kernel has no exact state-space representation in this ",
        "package; use method = \"exact\". Sums and scaled kernels of ",
        "Matern-1/2, 3/2, and 5/2 kernels have one."
      ),
      kernel$type
    ))
  )
}


# The stacked representation: feedback F, stationary covariance P, process
# noise L q L', observation H, and the transition exp(F dt) in closed form.
.state_space_representation <- function(kernel) {
  blocks <- .state_space_blocks(kernel)
  sizes <- vapply(blocks, function(block) block$dimension, integer(1))
  dimension <- sum(sizes)
  offsets <- cumsum(c(0L, sizes))[seq_along(sizes)]
  index <- lapply(seq_along(blocks), function(b) offsets[[b]] + seq_len(sizes[[b]]))

  block_diagonal <- function(make) {
    matrix_value <- matrix(0, dimension, dimension)

    for (b in seq_along(blocks)) {
      matrix_value[index[[b]], index[[b]]] <- make(blocks[[b]])
    }

    matrix_value
  }

  # Powers of N = F + lambda I for each block.
  nilpotent_powers <- lapply(blocks, function(block) {
    nilpotent <- block$feedback + block$decay * diag(block$dimension)
    powers <- list(diag(block$dimension))

    for (power in seq_len(block$dimension - 1L)) {
      powers[[power + 1L]] <- powers[[power]] %*% nilpotent
    }

    powers
  })

  list(
    dimension = dimension,
    feedback = block_diagonal(function(block) block$feedback),
    stationary = block_diagonal(function(block) block$stationary),
    process_noise = block_diagonal(function(block) {
      block$spectral_density * tcrossprod(block$noise_effect)
    }),
    observation = unlist(lapply(blocks, function(block) block$observation)),
    transition = function(step) {
      result <- matrix(0, dimension, dimension)

      for (b in seq_along(blocks)) {
        powers <- nilpotent_powers[[b]]
        exponential <- powers[[1L]]

        for (power in seq_along(powers)[-1L]) {
          exponential <- exponential +
            powers[[power]] * step^(power - 1L) / factorial(power - 1L)
        }

        result[index[[b]], index[[b]]] <- exp(-blocks[[b]]$decay * step) *
          exponential
      }

      result
    }
  )
}


# Kalman filter over sorted `times`. Entries with `observed` FALSE are only
# predicted. Updates use the Joseph form and symmetrization.
.kalman_filter <- function(representation, times, y, noise, observed) {
  n <- length(times)
  d <- representation$dimension
  h <- representation$observation
  stationary <- representation$stationary
  identity <- diag(d)
  predicted_mean <- matrix(0, d, n)
  filtered_mean <- matrix(0, d, n)
  predicted_covariance <- array(0, c(d, d, n))
  filtered_covariance <- array(0, c(d, d, n))
  transitions <- array(0, c(d, d, n))
  mean <- numeric(d)
  covariance <- stationary
  log_likelihood <- 0
  transitions[, , 1L] <- identity

  for (k in seq_len(n)) {
    if (k > 1L) {
      transition <- representation$transition(times[[k]] - times[[k - 1L]])
      mean <- drop(transition %*% mean)
      covariance <- tcrossprod(transition %*% (covariance - stationary),
                               transition) + stationary
      covariance <- 0.5 * (covariance + t(covariance))
      transitions[, , k] <- transition
    }

    predicted_mean[, k] <- mean
    predicted_covariance[, , k] <- covariance

    if (observed[[k]]) {
      projected <- drop(covariance %*% h)
      innovation_variance <- sum(h * projected) + noise[[k]]

      if (!isTRUE(innovation_variance > 0)) {
        .gp_abort(
          "A Kalman innovation variance is not positive.",
          "gaussianprocesses_numerical_error"
        )
      }

      innovation <- y[[k]] - sum(h * mean)
      gain <- projected / innovation_variance
      mean <- mean + gain * innovation
      joseph <- identity - outer(gain, h)
      covariance <- tcrossprod(joseph %*% covariance, joseph) +
        noise[[k]] * tcrossprod(gain)
      covariance <- 0.5 * (covariance + t(covariance))
      log_likelihood <- log_likelihood - 0.5 * (
        log(2 * pi * innovation_variance) +
          innovation^2 / innovation_variance
      )
    }

    filtered_mean[, k] <- mean
    filtered_covariance[, , k] <- covariance
  }

  list(
    predicted_mean = predicted_mean,
    predicted_covariance = predicted_covariance,
    filtered_mean = filtered_mean,
    filtered_covariance = filtered_covariance,
    transitions = transitions,
    log_likelihood = log_likelihood
  )
}


# The smoother gain G_k = P_k^f A_(k+1)' (P_(k+1)^p)^{-1}.
.smoother_gain <- function(filter, k) {
  t(solve(
    filter$predicted_covariance[, , k + 1L],
    filter$transitions[, , k + 1L] %*% filter$filtered_covariance[, , k]
  ))
}


# Rauch-Tung-Striebel smoother; also returns the gains, which give the
# posterior cross-covariances Cov(x_i, x_j) = G_i ... G_(j-1) P_j^s.
.rts_smoother <- function(filter) {
  n <- ncol(filter$filtered_mean)
  d <- nrow(filter$filtered_mean)
  smoothed_mean <- filter$filtered_mean
  smoothed_covariance <- filter$filtered_covariance
  gains <- array(0, c(d, d, n))

  for (k in rev(seq_len(n - 1L))) {
    gain <- .smoother_gain(filter, k)
    smoothed_mean[, k] <- filter$filtered_mean[, k] +
      drop(gain %*% (smoothed_mean[, k + 1L] - filter$predicted_mean[, k + 1L]))
    covariance <- filter$filtered_covariance[, , k] +
      tcrossprod(gain %*% (smoothed_covariance[, , k + 1L] -
                             filter$predicted_covariance[, , k + 1L]), gain)
    smoothed_covariance[, , k] <- 0.5 * (covariance + t(covariance))
    gains[, , k] <- gain
  }

  list(mean = smoothed_mean, covariance = smoothed_covariance, gains = gains)
}


# Filters training and query inputs together, query entries unobserved, in
# a stable sort; returns the filter and the positions of both in the sorted
# sequence.
.state_space_pass <- function(model, query = numeric(0)) {
  times <- c(drop(model$x), query)
  observed <- c(rep(TRUE, model$n_observations), rep(FALSE, length(query)))
  values <- c(model$y - model$prior_mean, rep(NA_real_, length(query)))
  noise <- c(model$training_noise_variance, rep(NA_real_, length(query)))
  sorted <- order(times, !observed)
  representation <- .state_space_representation(model$kernel)
  filter <- .kalman_filter(representation, times[sorted], values[sorted],
                           noise[sorted], observed[sorted])
  position <- integer(length(times))
  position[sorted] <- seq_along(sorted)

  list(
    representation = representation,
    filter = filter,
    training = position[seq_len(model$n_observations)],
    query = position[model$n_observations + seq_along(query)]
  )
}


.prepare_state_space_model <- function(model, name = "model") {
  .prepare_model(model, "gaussianprocesses_state_space_model", name)
}


# A state-space model at given hyperparameters on a numeric input. Inputs
# may be unsorted and tied. Without `smooth`, only the log marginal
# likelihood is computed, as the optimizer needs.
.fit_state_space_model <- function(x, y, kernel, noise_variance = 1e-6,
                                   mean = zero_mean(), smooth = TRUE) {
  data <- .validate_gp_training_data(x, y)

  if (ncol(data$x) != 1L) {
    .state_space_error("State-space inference needs one input dimension.")
  }

  .assert_kernel(kernel)
  .assert_kernel_dimensions(kernel, 1L)
  .assert_mean_function(mean)
  .assert_fixed_mean(mean, "State-space inference")
  representation <- .state_space_representation(kernel)
  training_noise <- .resolve_training_noise_variance(noise_variance,
                                                     nrow(data$x))

  if (any(training_noise <= 0)) {
    .gp_abort(
      "State-space inference needs positive noise variances.",
      "gaussianprocesses_validation_error"
    )
  }

  mean <- .resolve_mean(mean, data$x)
  model <- list(
    x = data$x,
    y = data$y,
    kernel = kernel,
    mean = mean,
    prior_mean = evaluate_mean(mean, data$x),
    noise_variance = as.numeric(noise_variance),
    training_noise_variance = training_noise,
    noise_structure = if (length(noise_variance) == 1L) {
      "homoscedastic"
    } else {
      "known_heteroscedastic"
    },
    method = "state_space",
    state_dimension = representation$dimension,
    n_observations = nrow(data$x),
    n_features = 1L,
    schema_version = .model_schema_version,
    package_version = .package_version()
  )
  pass <- .state_space_pass(model)
  model$log_marginal_likelihood <- pass$filter$log_likelihood

  if (!smooth) {
    return(structure(model, class = "gaussianprocesses_state_space_model"))
  }

  smoother <- .rts_smoother(pass$filter)
  h <- pass$representation$observation
  positions <- pass$training
  model$fitted_mean <- model$prior_mean +
    drop(crossprod(h, smoother$mean[, positions, drop = FALSE]))
  model$fitted_variance <- vapply(
    positions,
    function(k) sum(h * drop(smoother$covariance[, , k] %*% h)),
    numeric(1)
  )

  structure(model, class = "gaussianprocesses_state_space_model")
}


# Arguments of fit_gp() and optimize_gp() that state-space inference does
# not use are refused rather than ignored.
.assert_no_exact_arguments <- function(...) {
  extra <- list(...)

  if (length(extra) == 0L) {
    return(invisible(NULL))
  }

  labels <- names(extra)

  if (is.null(labels)) {
    labels <- character(length(extra))
  }

  labels <- ifelse(nzchar(labels), paste0("'", labels, "'"), "unnamed")
  .gp_abort(
    sprintf(
      "Argument(s) %s apply only to method = \"exact\".",
      paste(labels, collapse = ", ")
    ),
    "gaussianprocesses_validation_error"
  )
}


.fit_state_space_gp <- function(x, y, kernel, noise_variance = 1e-6,
                                mean = zero_mean(), ...) {
  .assert_no_exact_arguments(...)
  .fit_state_space_model(x, y, kernel, noise_variance, mean)
}


# Maximizes the Kalman-filter log marginal likelihood with the optimizer of
# optimize_gp(); L-BFGS-B approximates the gradient by finite differences.
.optimize_state_space_gp <- function(
    x,
    y,
    kernel,
    noise_variance = 1e-6,
    mean = zero_mean(),
    optimize_noise = TRUE,
    n_starts = 3L,
    start_spread = 1,
    lower = 1e-8,
    upper = 1e8,
    control = list(maxit = 500),
    ...) {
  .assert_no_exact_arguments(...)
  data <- .validate_gp_training_data(x, y)
  .assert_kernel(kernel)
  .state_space_representation(kernel)
  .assert_mean_function(mean)
  .assert_fixed_mean(mean, "State-space inference")
  training_noise <- .resolve_training_noise_variance(noise_variance,
                                                     nrow(data$x))

  if (!is.logical(optimize_noise) || length(optimize_noise) != 1L ||
      is.na(optimize_noise)) {
    .gp_abort(
      "'optimize_noise' must be TRUE or FALSE.",
      "gaussianprocesses_validation_error"
    )
  }

  n_starts <- .assert_positive_integer(n_starts, "n_starts")
  .assert_numeric_scalar(start_spread, "start_spread", lower = 0)
  bounds <- .validate_optimization_bounds(lower, upper)
  .assert_optimizer_control(control)
  noise_pattern <- if (optimize_noise && length(noise_variance) > 1L) {
    training_noise
  }

  if (any(training_noise <= 0)) {
    .gp_abort(
      "State-space inference needs positive noise variances.",
      "gaussianprocesses_validation_error"
    )
  }

  noise_parameter <- .noise_parameter(noise_pattern)
  extra_initial <- if (!optimize_noise) {
    numeric(0)
  } else if (is.null(noise_pattern)) {
    c(noise_variance = max(noise_variance, bounds[["lower"]]))
  } else {
    c(noise_scale = 1)
  }
  setup <- .optimization_setup(
    kernel,
    data$x,
    extra_initial = extra_initial,
    extra_constraints = .optimization_constraints(
      kernel,
      optimize_noise,
      noise_parameter
    )[names(extra_initial)],
    bounds = bounds,
    n_starts = n_starts,
    start_spread = start_spread
  )
  decode <- function(theta) {
    .decode_optimization_parameters(
      theta,
      kernel = kernel,
      optimize_noise = optimize_noise,
      fixed_noise_variance = noise_variance,
      constraints = setup$constraints,
      noise_pattern = noise_pattern
    )
  }

  evaluations <- c(objective = 0L, gradient = 0L, fits = 0L)
  problem <- list(
    objective = function(theta) {
      evaluations[["objective"]] <<- evaluations[["objective"]] + 1L
      evaluations[["fits"]] <<- evaluations[["fits"]] + 1L
      decoded <- decode(theta)
      value <- tryCatch(
        .fit_state_space_model(data$x, data$y, decoded$kernel,
                               decoded$noise_variance, mean,
                               smooth = FALSE)$log_marginal_likelihood,
        error = function(error) NA_real_
      )

      if (!is.finite(value)) {
        return(1e100)
      }

      -value
    },
    evaluations = function() evaluations,
    reset = function() {
      evaluations[] <<- 0L
      invisible(NULL)
    }
  )

  run <- .run_optimization_starts(problem, setup, control, "numerical")
  decoded <- decode(run$best_theta)
  fitted <- .fit_state_space_model(data$x, data$y, decoded$kernel,
                                   decoded$noise_variance, mean)
  fitted$optimization <- .optimization_record(
    setup,
    run,
    gradient = "numerical",
    bounds = bounds,
    n_starts = n_starts,
    start_spread = start_spread,
    specific = list(
      optimize_noise = optimize_noise,
      noise_parameterization = if (!optimize_noise) {
        "fixed"
      } else if (is.null(noise_pattern)) {
        "variance"
      } else {
        "scale"
      }
    ),
    results = list(log_marginal_likelihood = fitted$log_marginal_likelihood)
  )

  fitted
}


# Posterior prediction at numeric inputs, with the fields of predict_gp().
.predict_state_space <- function(
    model,
    x,
    interval_level = 0.95,
    include_covariance = FALSE,
    variance_tolerance = sqrt(.Machine$double.eps),
    observation_noise_variance = NULL) {
  prediction_x <- .validate_prediction_inputs(model, x)
  .validate_interval_level(interval_level)
  .assert_numeric_scalar(variance_tolerance, "variance_tolerance", lower = 0)

  if (!is.logical(include_covariance) || length(include_covariance) != 1L ||
      is.na(include_covariance)) {
    .gp_abort(
      "'include_covariance' must be TRUE or FALSE.",
      "gaussianprocesses_validation_error"
    )
  }

  query <- drop(prediction_x)
  noise <- .resolve_prediction_noise_variance(model, observation_noise_variance,
                                              length(query))
  pass <- .state_space_pass(model, query)
  smoother <- .rts_smoother(pass$filter)
  h <- pass$representation$observation
  positions <- pass$query
  prior_variance <- kernel_diagonal(model$kernel, prediction_x)

  latent_mean <- evaluate_mean(model$mean, prediction_x) +
    drop(crossprod(h, smoother$mean[, positions, drop = FALSE]))
  latent_variance <- .sanitize_posterior_variance(
    vapply(positions, function(k) {
      sum(h * drop(smoother$covariance[, , k] %*% h))
    }, numeric(1)),
    reference_scale = max(abs(prior_variance)),
    tolerance = variance_tolerance
  )
  latent_sd <- sqrt(latent_variance)
  observation_sd <- sqrt(latent_variance + noise)
  quantile <- stats::qnorm((1 + interval_level) / 2)

  result <- list(
    x = prediction_x,
    mean = latent_mean,
    latent_variance = latent_variance,
    observation_noise_variance = noise,
    observation_variance = latent_variance + noise,
    latent_sd = latent_sd,
    observation_sd = observation_sd,
    interval_level = interval_level,
    latent_interval = cbind(lower = latent_mean - quantile * latent_sd,
                            upper = latent_mean + quantile * latent_sd),
    prediction_interval = cbind(
      lower = latent_mean - quantile * observation_sd,
      upper = latent_mean + quantile * observation_sd
    )
  )

  if (include_covariance) {
    covariance <- .state_space_query_covariance(smoother, h, positions)
    diag(covariance) <- latent_variance
    result$latent_covariance <- covariance
    result$observation_covariance <- covariance + diag(noise, length(query))
  }

  structure(result, class = "gaussianprocesses_prediction")
}


# Posterior covariance between query positions from the smoother's
# cross-covariances Cov(x_i, x_j) = G_i ... G_(j-1) P_j^s for i < j. The
# gains are multiplied once per segment between consecutive queries, and
# each row carries the vector h' G_i ... G_(j-1), so m queries among N
# sequence entries cost O(N + m^2) small products.
.state_space_query_covariance <- function(smoother, h, positions) {
  m <- length(positions)
  ordered <- order(positions)
  sorted <- positions[ordered]
  segments <- vector("list", max(m - 1L, 0L))

  for (a in seq_len(m - 1L)) {
    product <- diag(length(h))

    for (k in seq.int(sorted[[a]], sorted[[a + 1L]] - 1L)) {
      product <- product %*% smoother$gains[, , k]
    }

    segments[[a]] <- product
  }

  projected <- lapply(sorted, function(j) {
    drop(smoother$covariance[, , j] %*% h)
  })
  covariance <- matrix(0, m, m)

  for (a in seq_len(m)) {
    covariance[a, a] <- sum(h * projected[[a]])
    row <- h

    for (b in a + seq_len(m - a)) {
      row <- drop(row %*% segments[[b - 1L]])
      covariance[a, b] <- sum(row * projected[[b]])
      covariance[b, a] <- covariance[a, b]
    }
  }

  result <- matrix(0, m, m)
  result[ordered, ordered] <- covariance
  result
}


# Square root R with R'R = covariance for small positive semidefinite
# matrices, robust to singular ones (tied inputs give zero covariance).
.semidefinite_root <- function(covariance) {
  decomposition <- eigen(0.5 * (covariance + t(covariance)), symmetric = TRUE)
  values <- pmax(decomposition$values, 0)
  t(decomposition$vectors %*% diag(sqrt(values), length(values)))
}


#' @rdname sample_gp_posterior
#' @export
sample_gp_posterior.gaussianprocesses_state_space_model <- function(
    model,
    x,
    n_draws = 1L,
    type = c("latent", "observation"),
    observation_noise_variance = NULL,
    seed = NULL,
    ...) {
  model <- .prepare_state_space_model(model)
  type <- match.arg(type)
  n_draws <- .assert_draw_count(n_draws)
  .validate_seed(seed)
  prediction_x <- .validate_prediction_inputs(model, x)
  query <- drop(prediction_x)

  if (identical(type, "latent") && !is.null(observation_noise_variance)) {
    .gp_abort(
      paste0(
        "'observation_noise_variance' is used only for observation draws; ",
        "use type = \"observation\"."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  noise <- if (identical(type, "observation")) {
    .resolve_prediction_noise_variance(model, observation_noise_variance,
                                       length(query))
  }
  pass <- .state_space_pass(model, query)
  filter <- pass$filter
  smoother <- .rts_smoother(filter)
  h <- pass$representation$observation
  positions <- pass$query
  n_steps <- ncol(filter$filtered_mean)
  prior_mean <- evaluate_mean(model$mean, prediction_x)

  # Forward filtering, backward sampling: x_N from the last filtered state,
  # then x_k | x_(k+1) ~ N(m_k + G_k (x_(k+1) - m_(k+1)^p),
  # P_k - G_k P_(k+1)^p G_k').
  draws <- .with_local_seed(seed, function() {
    standard <- function() {
      matrix(stats::rnorm(length(h) * n_draws), length(h))
    }
    # Only the current state is kept; steps before the first query are not
    # needed.
    query_at <- integer(n_steps)
    query_at[positions] <- seq_along(positions)
    latent <- matrix(0, length(query), n_draws)
    state <- filter$filtered_mean[, n_steps] +
      crossprod(.semidefinite_root(filter$filtered_covariance[, , n_steps]),
                standard())

    for (k in rev(seq.int(min(positions), n_steps))) {
      if (k < n_steps) {
        gain <- smoother$gains[, , k]
        conditional_mean <- filter$filtered_mean[, k] +
          gain %*% (state - filter$predicted_mean[, k + 1L])
        conditional_covariance <- filter$filtered_covariance[, , k] -
          tcrossprod(gain %*% filter$predicted_covariance[, , k + 1L], gain)
        state <- conditional_mean +
          crossprod(.semidefinite_root(conditional_covariance), standard())
      }

      if (query_at[[k]] > 0L) {
        latent[query_at[[k]], ] <- prior_mean[[query_at[[k]]]] +
          drop(crossprod(h, state))
      }
    }

    observation <- if (identical(type, "observation")) {
      latent + sqrt(noise) *
        matrix(stats::rnorm(length(query) * n_draws), length(query))
    }

    list(latent = latent, observation = observation)
  })

  .new_gp_samples(
    x = prediction_x,
    distribution = "posterior",
    model = "state space",
    latent = draws$latent,
    observation = draws$observation,
    mean = prior_mean +
      drop(crossprod(h, smoother$mean[, positions, drop = FALSE])),
    latent_covariance = NULL,
    observation_noise_variance = noise,
    decomposition = list(jitter = 0, relative_jitter = 0,
                         covariance_scale = NULL, attempts = 0L),
    seed = seed
  )
}


#' @export
print.gaussianprocesses_state_space_model <- function(x, ...) {
  x <- .prepare_state_space_model(x, "x")
  cat("Gaussian-process regression model (state space, Kalman smoothing)\n")
  cat("  observations: ", x$n_observations, "\n", sep = "")
  cat("  state dimension: ", x$state_dimension, "\n", sep = "")
  cat("  mean: ", .format_mean(x$mean), "\n", sep = "")
  cat("  noise structure: ", x$noise_structure, "\n", sep = "")
  cat(
    "  log marginal likelihood: ",
    format(signif(x$log_marginal_likelihood, 7), trim = TRUE),
    "\n",
    sep = ""
  )

  if (!is.null(x$optimization)) {
    cat(
      "  optimization: ",
      if (isTRUE(x$optimization$converged)) "converged" else "not converged",
      " (", x$optimization$method, ", ", x$optimization$gradient,
      " gradient)\n",
      sep = ""
    )
  }

  cat("  kernel:\n")
  cat(paste0("    ", .format_kernel(x$kernel)), sep = "\n")
  cat("\n")
  invisible(x)
}
