# Sparse Gaussian-process regression -----------------------------------------

.scale_inputs_for_inducing_selection <- function(x) {
  center <- colMeans(x)

  if (nrow(x) < 2L) {
    scale <- rep(1, ncol(x))
  } else {
    scale <- apply(x, 2L, stats::sd)
    scale[!is.finite(scale) | scale == 0] <- 1
  }

  sweep(
    sweep(x, 2L, center, FUN = "-"),
    2L,
    scale,
    FUN = "/"
  )
}


.select_farthest_inducing_indices <- function(x, n_inducing) {
  # Standardize with every row, but only offer the first occurrence of each
  # distinct input: repeated rows would make K_uu exactly singular.
  candidates <- which(!duplicated(x))
  scaled <- .scale_inputs_for_inducing_selection(x)[
    candidates,
    ,
    drop = FALSE
  ]

  first <- which.min(rowSums(scaled^2))
  selected <- first

  if (n_inducing == 1L) {
    return(candidates[selected])
  }

  minimum_distance <- rowSums(
    sweep(scaled, 2L, scaled[first, ], FUN = "-")^2
  )
  minimum_distance[first] <- -Inf

  while (length(selected) < n_inducing) {
    next_index <- which.max(minimum_distance)
    selected <- c(selected, next_index)

    distance_to_new <- rowSums(
      sweep(
        scaled,
        2L,
        scaled[next_index, ],
        FUN = "-"
      )^2
    )

    minimum_distance <- pmin(
      minimum_distance,
      distance_to_new
    )
    minimum_distance[selected] <- -Inf
  }

  candidates[selected]
}


.select_quantile_inducing_indices <- function(x, n_inducing) {
  if (ncol(x) != 1L) {
    .gp_abort(
      "'method = quantile' is available only for one-dimensional inputs.",
      "gaussianprocesses_validation_error"
    )
  }

  targets <- as.numeric(
    stats::quantile(
      x[, 1L],
      probs = seq(0, 1, length.out = n_inducing),
      names = FALSE,
      type = 8
    )
  )

  available <- which(!duplicated(x[, 1L]))
  selected <- integer(0)

  for (target in targets) {
    local_index <- which.min(
      abs(x[available, 1L] - target)
    )
    selected <- c(selected, available[[local_index]])
    available <- available[-local_index]
  }

  selected
}


# Greedy conditional-variance selection: an incremental pivoted Cholesky
# decomposition of K(X, X). Column t of the n x m factor is
#   (K(X, x_i) - L[, 1:(t-1)] L[i, 1:(t-1)]') / sqrt(d_i)
# for the pivot i with the largest residual variance d_i, and the residuals
# fall by its square. Only the chosen columns of K(X, X) are evaluated. The
# product uses the whole factor, whose later columns are still zero: taking
# the first t - 1 columns would copy them at every step, which was slower.
.select_variance_inducing_indices <- function(x, n_inducing, kernel) {
  n <- nrow(x)
  residual <- kernel_diagonal(kernel, x)
  threshold <- 1e-12 * max(abs(residual))
  factor <- matrix(0, n, n_inducing)
  indices <- integer(0)

  for (step in seq_len(n_inducing)) {
    pivot <- which.max(residual)

    if (residual[[pivot]] <= threshold) {
      warning(
        warningCondition(
          sprintf(
            paste0(
              "The prior variance at the inputs is explained by %d ",
              "inducing point(s); returning %d instead of %d."
            ),
            step - 1L,
            step - 1L,
            n_inducing
          ),
          class = "gaussianprocesses_inducing_warning"
        )
      )
      break
    }

    column <- drop(evaluate_kernel(kernel, x, x[pivot, , drop = FALSE]))

    if (step > 1L) {
      column <- column - drop(factor %*% factor[pivot, ])
    }

    factor[, step] <- column / sqrt(residual[[pivot]])
    residual <- residual - factor[, step]^2
    residual[[pivot]] <- 0
    indices <- c(indices, pivot)
  }

  indices
}


#' Select inducing points deterministically
#'
#' @param x Numeric vector or matrix of candidate inputs.
#' @param n_inducing Number of inducing points to select.
#' @param method Selection method. `"farthest"` performs deterministic
#'   farthest-point sampling after column standardization. `"quantile"`
#'   selects points near equally spaced empirical quantiles and is available
#'   only for one-dimensional inputs. `"variance"` performs greedy
#'   conditional-variance selection under `kernel`.
#' @param kernel Kernel specification, needed by `method = "variance"`.
#'
#' @details
#' Repeated input rows are offered as candidates only once, so the selected
#' inducing locations are always distinct. `n_inducing` therefore cannot
#' exceed the number of distinct rows in `x`.
#'
#' Greedy conditional-variance selection (Burt, Rasmussen, and van der Wilk,
#' 2019, 2020) adds, one at a time, the input whose prior variance is least
#' explained by the points chosen so far: the largest residual
#' \eqn{k(x, x) - Q(x, x)} with \eqn{Q} the Nyström approximation from the
#' chosen points. This is the pivot order of a pivoted Cholesky
#' decomposition of \eqn{K(X, X)}, computed incrementally without forming
#' \eqn{K(X, X)} in \eqn{O(n m^2)} time and \eqn{O(n m)} memory. It
#' minimizes the trace term \eqn{\mathrm{tr}(K_{ff} - Q_{ff})} greedily, so
#' it is the natural choice for VFE, whose bound loses half that trace over
#' the noise variance. If the prior variance is explained, to a relative
#' \eqn{10^{-12}}, before `n_inducing` points are chosen, the selection stops
#' with a warning of class `gaussianprocesses_inducing_warning` and returns
#' fewer points.
#'
#' @references Burt, D. R., Rasmussen, C. E., and van der Wilk, M. (2019).
#'   Rates of convergence for sparse variational Gaussian process regression.
#'   *Proceedings of the 36th International Conference on Machine Learning*,
#'   862--871.
#'
#'   Burt, D. R., Rasmussen, C. E., and van der Wilk, M. (2020). Convergence
#'   of sparse variational inference in Gaussian processes regression.
#'   *Journal of Machine Learning Research*, 21(131), 1--63.
#'
#' @return A numeric matrix with one inducing point per row. The selected
#'   training-row indices are stored in the `indices` attribute.
#' @examples
#' x <- seq(0, 10, length.out = 101)
#'
#' inducing <- select_inducing_points(x, n_inducing = 5, method = "quantile")
#' inducing[, 1]
#' attr(inducing, "indices")
#' @template experimental
#' @export
select_inducing_points <- function(
    x,
    n_inducing,
    method = c("farthest", "quantile", "variance"),
    kernel = NULL) {
  x_matrix <- .as_input_matrix(x, name = "x")
  n_inducing <- .assert_positive_integer(
    n_inducing,
    "n_inducing"
  )
  method <- match.arg(method)
  n_distinct <- sum(!duplicated(x_matrix))

  if (n_inducing > n_distinct) {
    .gp_abort(
      sprintf(
        paste0(
          "'n_inducing' cannot exceed the number of distinct candidate ",
          "inputs (%d)."
        ),
        n_distinct
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (identical(method, "variance")) {
    if (is.null(kernel)) {
      .gp_abort(
        "Variance selection needs a 'kernel'.",
        "gaussianprocesses_validation_error"
      )
    }

    .assert_kernel(kernel)
    .assert_kernel_dimensions(kernel, ncol(x_matrix))
  }

  indices <- switch(
    method,
    farthest = .select_farthest_inducing_indices(
      x_matrix,
      n_inducing
    ),
    quantile = .select_quantile_inducing_indices(
      x_matrix,
      n_inducing
    ),
    variance = .select_variance_inducing_indices(
      x_matrix,
      n_inducing,
      kernel
    )
  )

  result <- x_matrix[indices, , drop = FALSE]
  attr(result, "indices") <- indices
  attr(result, "method") <- method
  result
}


.triangular_whiten <- function(factor, rhs) {
  # Returns R^{-T} rhs for an upper-triangular Cholesky factor of
  # S = R'R, so that crossprod(result) equals rhs' S^{-1} rhs.
  forwardsolve(
    t(factor),
    rhs,
    upper.tri = FALSE,
    transpose = FALSE
  )
}


.validate_inducing_points <- function(
    inducing_points,
    n_features) {
  matrix <- .as_input_matrix(
    inducing_points,
    name = "inducing_points"
  )

  if (ncol(matrix) != n_features) {
    .gp_abort(
      sprintf(
        paste0(
          "'inducing_points' must have %d input dimension(s). ",
          "Received %d."
        ),
        n_features,
        ncol(matrix)
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (anyDuplicated(data.frame(matrix))) {
    .gp_abort(
      "'inducing_points' must contain unique rows.",
      "gaussianprocesses_validation_error"
    )
  }

  matrix
}


# Inducing points: explicit, or selected from the training inputs.
.sparse_inducing_points <- function(x, n_inducing, inducing_points,
                                    selection, kernel) {
  if (!is.null(inducing_points)) {
    return(list(
      points = .validate_inducing_points(inducing_points, ncol(x)),
      indices = NULL,
      method = "explicit"
    ))
  }

  if (is.null(n_inducing)) {
    n_inducing <- min(20L, sum(!duplicated(x)))
  }

  n_inducing <- .assert_positive_integer(n_inducing, "n_inducing")
  points <- select_inducing_points(x, n_inducing = n_inducing,
                                   method = selection, kernel = kernel)

  list(
    points = points,
    indices = attr(points, "indices"),
    method = attr(points, "method")
  )
}


.assert_vfe_noise <- function(method, training_noise) {
  if (identical(method, "vfe") && any(training_noise <= 0)) {
    .gp_abort(
      "VFE needs positive noise variances; its bound divides by them.",
      "gaussianprocesses_validation_error"
    )
  }

  invisible(training_noise)
}


.sparse_method <- function(model) {
  if (is.null(model$method)) "fitc" else model$method
}


.sparse_label <- function(method) {
  if (identical(method, "vfe")) "VFE" else "FITC"
}


# A sparse model at given hyperparameters and inducing points; validation is
# left to the callers. Both approximations use
#   C = Q_ff + Lambda,   Q_ff = K_fu K_uu^{-1} K_uf = V'V,   V = R_uu^{-T} K_uf,
# with Lambda = diag(k_ff - q_ff) + noise for FITC and Lambda = noise for
# VFE, whose objective subtracts tr(K_ff - Q_ff) / 2 for each unit of noise
# variance (Titsias, 2009).
.fit_sparse_model <- function(
    x,
    y,
    kernel,
    noise_variance,
    training_noise,
    mean,
    inducing,
    inducing_indices,
    selection_method,
    method,
    fitc_diagonal_floor,
    policy) {
  m <- nrow(inducing)
  prior_mean <- evaluate_mean(mean, x)
  centered_y <- y - prior_mean

  kuu <- evaluate_kernel(kernel, inducing)
  kuu_decomposition <- .stable_cholesky(
    kuu,
    initial_jitter = policy$initial_jitter,
    fallback_jitter = policy$fallback_jitter,
    jitter_multiplier = policy$jitter_multiplier,
    max_attempts = policy$max_attempts,
    tolerance = policy$symmetry_tolerance
  )

  kuf <- evaluate_kernel(kernel, inducing, x)

  # Every quantity below uses the same, possibly jittered, inducing
  # covariance R_uu' R_uu through V = R_uu^{-T} K_uf, so Q_ff = V'V.
  whitened_kuf <- .triangular_whiten(kuu_decomposition$factor, kuf)

  qff_diagonal <- colSums(whitened_kuf^2)
  prior_diagonal <- kernel_diagonal(kernel, x)

  residual_diagonal <- .sanitize_posterior_variance(
    prior_diagonal - qff_diagonal,
    reference_scale = max(abs(prior_diagonal)),
    tolerance = policy$symmetry_tolerance
  )

  if (identical(method, "vfe")) {
    lambda <- training_noise
    clipped_lambda <- rep(FALSE, length(lambda))
  } else {
    lambda_raw <- residual_diagonal + training_noise
    clipped_lambda <- lambda_raw < fitc_diagonal_floor
    lambda <- pmax(lambda_raw, fitc_diagonal_floor)
  }

  # Factorize B = I + V Lambda^{-1} V', whose eigenvalues are at least one,
  # instead of A = K_uu + K_uf Lambda^{-1} K_fu = R_uu' B R_uu. A inherits
  # the conditioning of K_uu, needed its own jitter, and then broke the
  # determinant identity whenever that jitter differed from K_uu's.
  scaled_whitened_kuf <- sweep(whitened_kuf, 2L, sqrt(lambda), FUN = "/")
  b_matrix <- diag(m) + tcrossprod(scaled_whitened_kuf)
  b_decomposition <- .stable_cholesky(
    b_matrix,
    initial_jitter = policy$initial_jitter,
    fallback_jitter = policy$fallback_jitter,
    jitter_multiplier = policy$jitter_multiplier,
    max_attempts = policy$max_attempts,
    tolerance = policy$symmetry_tolerance
  )

  # B^{-1} V Lambda^{-1} y; beta = A^{-1} K_uf Lambda^{-1} y = R_uu^{-1} of it.
  whitened_beta <- .cholesky_solve(
    b_decomposition$factor,
    drop(whitened_kuf %*% (centered_y / lambda))
  )
  beta <- backsolve(
    kuu_decomposition$factor,
    whitened_beta,
    upper.tri = TRUE,
    transpose = FALSE
  )

  training_alpha <- (
    centered_y - drop(crossprod(whitened_kuf, whitened_beta))
  ) / lambda

  # Matrix determinant lemma: |Q_ff + Lambda| = |Lambda| |B|.
  log_determinant <- sum(log(lambda)) +
    2 * sum(log(diag(b_decomposition$factor)))

  log_density <- -0.5 * drop(crossprod(centered_y, training_alpha)) -
    0.5 * log_determinant -
    0.5 * nrow(x) * log(2 * pi)

  # The inducing-point coverage diagnostic t = tr(K_ff - Q_ff) / sigma^2,
  # with one noise variance per observation; it is the VFE penalty times 2.
  trace_term <- if (all(training_noise > 0)) {
    sum(residual_diagonal / training_noise)
  } else {
    NA_real_
  }

  structure(
    list(
      x = x,
      y = y,
      centered_y = centered_y,
      kernel = kernel,
      mean = mean,
      prior_mean = prior_mean,
      noise_variance = as.numeric(noise_variance),
      training_noise_variance = training_noise,
      noise_structure = if (length(noise_variance) == 1L) {
        "homoscedastic"
      } else {
        "known_heteroscedastic"
      },
      method = method,
      inducing_points = inducing,
      inducing_indices = inducing_indices,
      selection_method = selection_method,
      n_inducing = m,
      kuu = kuu,
      kuu_cholesky = kuu_decomposition$factor,
      kuu_jitter = kuu_decomposition$jitter,
      kuf = kuf,
      qff_diagonal = qff_diagonal,
      residual_diagonal = residual_diagonal,
      lambda = lambda,
      fitc_diagonal_floor = fitc_diagonal_floor,
      fitc_floor_count = sum(clipped_lambda),
      b_cholesky = b_decomposition$factor,
      # Upper-triangular factor of A, since A = (R_B R_uu)' (R_B R_uu).
      a_cholesky = b_decomposition$factor %*% kuu_decomposition$factor,
      # Jitter applied to the inducing posterior system B.
      a_jitter = b_decomposition$jitter,
      beta = beta,
      training_alpha = training_alpha,
      trace_term = trace_term,
      inducing_collapse = .inducing_collapse(kuu),
      log_marginal_likelihood = if (identical(method, "vfe")) {
        log_density - 0.5 * trace_term
      } else {
        log_density
      },
      numerical_policy = c(
        policy,
        list(fitc_diagonal_floor = fitc_diagonal_floor)
      ),
      n_observations = nrow(x),
      n_features = ncol(x),
      approximation = .sparse_label(method),
      schema_version = .model_schema_version,
      package_version = .package_version()
    ),
    class = "gaussianprocesses_sparse_model"
  )
}


# Derivatives of k(x_i, x_i) in the optimizer coordinates of the kernel
# parameters, one column per parameter, from the diagonals of small blocks.
.kernel_diagonal_gradient <- function(kernel, x, block_size = 256L) {
  rows_of <- split(seq_len(nrow(x)), ceiling(seq_len(nrow(x)) / block_size))
  blocks <- lapply(rows_of, function(rows) {
    block <- x[rows, , drop = FALSE]
    gradient <- .kernel_value_and_gradient(kernel, block, block)$gradient
    matrix(
      vapply(gradient, diag, numeric(length(rows))),
      nrow = length(rows),
      dimnames = list(NULL, names(gradient))
    )
  })

  do.call(rbind, blocks)
}


# Gradient of the FITC objective or the VFE bound in the optimizer
# coordinates of the kernel parameters and, with `include_noise`, of the
# noise (`noise_parameter`: its variance, or the scale of a noise pattern).
#
# With C = Q + Lambda, alpha = C^{-1} r, and M = alpha alpha' - C^{-1}, the
# Gaussian term changes by tr(M (dQ + dLambda)) / 2, and with
# P = K_uu^{-1} K_uf,
#   dQ = dK_fu P + P' dK_uf - P' dK_uu P,
#   tr(M dQ) = 2 sum(PM * dK_uf) - sum(P M P' * dK_uu).
# In whitened form, V = R_uu^{-T} K_uf, W = V Lambda^{-1}, and
# B = I + V Lambda^{-1} V', V C^{-1} = B^{-1} W, so
#   P C^{-1} = R_uu^{-1} B^{-1} W,   P C^{-1} P' = R_uu^{-1} (I - B^{-1}) R_uu^{-T},
# and no n x n matrix is formed. The diagonal of M is
# alpha_i^2 - 1 / lambda_i + |R_B^{-T} v_i|^2 / lambda_i^2. FITC adds
# tr(M dLambda) / 2 with dLambda_ii = dk_ii - dq_ii (zero where the floor
# applies); VFE adds the derivative of -sum((k_ii - q_ii) / lambda_i) / 2.
# Terms over observations are accumulated in blocks of columns.
#
# Altogether the objective changes by sum(G_uf * dK_uf) + sum(G_uu * dK_uu)
# + sum_i w_i dk_ii, with G_uf = P M - 2 P W and G_uu = -P M P' / 2 + P W P'.
# With `include_inducing`, the derivatives in the inducing inputs Z follow:
# Z_jd changes only row j of K_uf and row and column j of K_uu, so
#   dF / dZ_jd = sum_i G_uf[j, i] dk(z_j, x_i) / dz_jd
#                + 2 sum_l G_uu[j, l] dk(z_j, z_l) / dz_jd,
# with the input derivatives of kernel_input_gradient(); they are named
# "inducing.<j>.<d>" in the column-major order of Z.
.sparse_gradient <- function(
    model,
    include_noise,
    noise_parameter = "noise_variance",
    include_inducing = FALSE) {
  method <- .sparse_method(model)
  kernel <- model$kernel
  inducing <- model$inducing_points
  x <- model$x
  m <- model$n_inducing
  n <- model$n_observations
  lambda <- model$lambda
  noise <- model$training_noise_variance
  alpha <- model$training_alpha
  r_uu <- model$kuu_cholesky
  r_b <- model$b_cholesky
  inverse_b <- chol2inv(r_b)

  # R^{-1} X R^{-T} for symmetric X.
  sandwich <- function(matrix) {
    left <- backsolve(r_uu, matrix)
    t(backsolve(r_uu, t(left)))
  }

  r_inverse_b_inverse <- backsolve(r_uu, inverse_b)
  p_alpha <- backsolve(
    r_uu,
    forwardsolve(t(r_uu), drop(model$kuf %*% alpha))
  )
  p_m_p <- tcrossprod(p_alpha) - sandwich(diag(m) - inverse_b)

  kuu_gradient <- .kernel_value_and_gradient(kernel, inducing,
                                             inducing)$gradient
  parameter_names <- names(kuu_gradient)
  diagonal_gradient <- .kernel_diagonal_gradient(kernel, x)
  active <- if (identical(method, "fitc")) {
    model$residual_diagonal + noise >= model$fitc_diagonal_floor
  } else {
    rep(TRUE, n)
  }

  # The diagonal terms enter only as sum_i w_i (dk_ii - dq_ii), with weights
  # w_i = M_ii / 2 on the active FITC entries and -1 / (2 lambda_i) for VFE.
  # With W = diag(w),
  #   sum_i w_i dq_ii = 2 sum(dK_uf * (P W)) - sum(dK_uu * (P W P')),
  # so P W P' is accumulated once for all parameters instead of forming
  # P_i' dK_uu P_i for every observation and parameter.
  low_rank <- stats::setNames(numeric(length(parameter_names)),
                              parameter_names)
  weighted_cross <- low_rank
  weighted_inducing <- matrix(0, m, m)
  weights <- numeric(n)
  noise_gradient <- 0
  block_size <- max(1L, floor(.prediction_block_elements / m))
  n_features <- ncol(inducing)
  inducing_gradient <- matrix(0, m, n_features)

  for (start in seq.int(1L, n, by = block_size)) {
    rows <- start:min(start + block_size - 1L, n)
    kuf_gradient <- .kernel_value_and_gradient(
      kernel,
      inducing,
      x[rows, , drop = FALSE]
    )$gradient
    whitened <- .triangular_whiten(r_uu, model$kuf[, rows, drop = FALSE])
    projected <- backsolve(r_uu, whitened)
    scaled <- sweep(whitened, 2L, lambda[rows], FUN = "/")
    p_m <- outer(p_alpha, alpha[rows]) - r_inverse_b_inverse %*% scaled
    precision_diagonal <- 1 / lambda[rows] -
      colSums(.triangular_whiten(r_b, whitened)^2) / lambda[rows]^2
    m_diagonal <- alpha[rows]^2 - precision_diagonal

    block_weights <- if (identical(method, "fitc")) {
      0.5 * m_diagonal * active[rows]
    } else {
      -0.5 / lambda[rows]
    }
    weights[rows] <- block_weights
    weighted_projected <- sweep(projected, 2L, block_weights, FUN = "*")
    weighted_inducing <- weighted_inducing +
      tcrossprod(weighted_projected, projected)

    for (name in parameter_names) {
      d_kuf <- kuf_gradient[[name]]
      low_rank[[name]] <- low_rank[[name]] + sum(p_m * d_kuf)
      weighted_cross[[name]] <- weighted_cross[[name]] +
        2 * sum(d_kuf * weighted_projected)
    }

    if (include_inducing) {
      cross_weight <- p_m - 2 * weighted_projected
      input_gradient <- kernel_input_gradient(
        kernel,
        inducing,
        x[rows, , drop = FALSE]
      )

      for (dimension in seq_len(n_features)) {
        inducing_gradient[, dimension] <- inducing_gradient[, dimension] +
          rowSums(cross_weight * input_gradient[[dimension]])
      }
    }

    noise_gradient <- noise_gradient + if (identical(method, "fitc")) {
      0.5 * sum((m_diagonal * noise[rows])[active[rows]])
    } else {
      0.5 * sum(m_diagonal * lambda[rows]) +
        0.5 * sum(model$residual_diagonal[rows] / lambda[rows])
    }
  }

  gradient <- vapply(
    parameter_names,
    function(name) {
      weighted_dq <- weighted_cross[[name]] -
        sum(kuu_gradient[[name]] * weighted_inducing)
      low_rank[[name]] - 0.5 * sum(p_m_p * kuu_gradient[[name]]) +
        sum(weights * diagonal_gradient[, name]) - weighted_dq
    },
    numeric(1)
  )

  if (include_noise) {
    gradient <- c(gradient, stats::setNames(noise_gradient, noise_parameter))
  }

  if (include_inducing) {
    inducing_weight <- -0.5 * p_m_p + weighted_inducing
    input_gradient <- kernel_input_gradient(kernel, inducing, inducing)

    for (dimension in seq_len(n_features)) {
      inducing_gradient[, dimension] <- inducing_gradient[, dimension] +
        2 * rowSums(inducing_weight * input_gradient[[dimension]])
    }

    gradient <- c(
      gradient,
      stats::setNames(
        as.vector(inducing_gradient),
        .inducing_parameter_names(m, n_features)
      )
    )
  }

  gradient
}


.inducing_parameter_names <- function(n_inducing, n_features) {
  paste0(
    "inducing.",
    rep(seq_len(n_inducing), times = n_features),
    ".",
    rep(seq_len(n_features), each = n_inducing)
  )
}


# Pairs of inducing points that are nearly identical under the kernel, with
# correlation k(z_j, z_l) / sqrt(k(z_j, z_j) k(z_l, z_l)) above
# 1 - `tolerance`; NULL if there are none. Such points add almost nothing to
# the approximation and make K_uu nearly singular.
.inducing_collapse <- function(kuu, tolerance = 1e-6) {
  scale <- sqrt(pmax(diag(kuu), .Machine$double.xmin))
  correlation <- kuu / outer(scale, scale)
  pairs <- which(
    upper.tri(correlation) & correlation > 1 - tolerance,
    arr.ind = TRUE
  )

  if (nrow(pairs) == 0L) {
    return(NULL)
  }

  data.frame(
    first = pairs[, "row"],
    second = pairs[, "col"],
    correlation = correlation[pairs]
  )
}


.warn_inducing_collapse <- function(model) {
  collapse <- model$inducing_collapse

  if (!is.null(collapse)) {
    warning(
      warningCondition(
        sprintf(
          paste0(
            "%d pair(s) of inducing points are nearly identical under the ",
            "kernel (largest correlation %s, first pair %d and %d). They add ",
            "almost nothing to the approximation and make K_uu nearly ",
            "singular; the factorization used jitter %s. Remove duplicates or ",
            "select the inducing points afresh."
          ),
          nrow(collapse),
          format(signif(max(collapse$correlation), 10)),
          collapse$first[[1L]],
          collapse$second[[1L]],
          format(signif(model$kuu_jitter, 3))
        ),
        class = "gaussianprocesses_inducing_collapse_warning"
      )
    )
  }

  invisible(model)
}


#' Fit a sparse Gaussian process (FITC or VFE)
#'
#' Fits a sparse approximation with fixed inducing locations: the Fully
#' Independent Training Conditional (FITC) approximation or the collapsed
#' variational free energy (VFE) of Titsias (2009). Exact GP regression
#' remains the reference model; these provide a lower-cost path when the
#' number of inducing points is much smaller than the number of observations.
#'
#' @param x Numeric vector or matrix of training inputs.
#' @param y Numeric response vector.
#' @param kernel Gaussian-process kernel specification.
#' @param noise_variance Non-negative scalar or one non-negative value per
#'   training observation.
#' @param mean Gaussian-process mean specification with fixed
#'   coefficients.
#' @param n_inducing Optional number of inducing points selected from `x`.
#' @param inducing_points Optional explicit inducing locations. If supplied,
#'   `n_inducing` is ignored.
#' @param selection Inducing-point selection method used when explicit points
#'   are not supplied: `"farthest"`, `"quantile"`, or `"variance"` (greedy
#'   conditional-variance selection under `kernel`; see
#'   [select_inducing_points()]).
#' @param method The approximation: `"fitc"` (the default) or `"vfe"`. VFE
#'   needs positive noise variances.
#' @param fitc_diagonal_floor Strictly positive numerical floor applied to the
#'   FITC diagonal term before inversion.
#' @param initial_jitter Non-negative jitter tried first in the
#'   inducing-system Cholesky factorizations, relative to the scale of each
#'   factorized matrix (the mean of its diagonal).
#' @param fallback_jitter Positive relative jitter tried after the first
#'   failed factorization.
#' @param jitter_multiplier Multiplicative jitter escalation factor.
#' @param max_attempts Maximum Cholesky attempts.
#' @param symmetry_tolerance Relative tolerance for checking that
#'   covariance matrices are symmetric.
#'
#' @details
#' Let \eqn{m} be the number of inducing points and \eqn{n} the number of
#' training observations. FITC replaces the full training covariance with a
#' low-rank inducing covariance plus a diagonal correction. The dominant
#' training cost is \eqn{O(n m^2 + m^3)} rather than the exact GP's
#' \eqn{O(n^3)}, with \eqn{m \ll n}. The implementation stores an
#' \eqn{m \times n} cross covariance rather than an \eqn{n \times n} training
#' covariance.
#'
#' FITC preserves each training point's prior marginal variance through the
#' diagonal correction, unlike simpler deterministic training conditional
#' approximations. It is still an approximation: inducing locations and their
#' number can affect posterior means, variances, and the approximate marginal
#' likelihood.
#'
#' If the inducing covariance needs numerical jitter, the jittered matrix is
#' used consistently in the low-rank term, the posterior, and the approximate
#' log marginal likelihood. The inducing posterior is factorized in whitened
#' form, whose eigenvalues are bounded below by one.
#'
#' @section VFE:
#' With \eqn{Q_{ff} = K_{fu} K_{uu}^{-1} K_{uf}} and noise
#' \eqn{\Sigma = \mathrm{diag}(\sigma_i^2)}, the collapsed variational
#' bound of Titsias (2009) is
#' \deqn{\mathcal F = \log N(y \mid m(X), Q_{ff} + \Sigma) -
#' \frac12 \sum_i \frac{k(x_i, x_i) - [Q_{ff}]_{ii}}{\sigma_i^2}
#' \le \log p(y).}
#' It is a lower bound on the exact log marginal likelihood, reported as
#' `log_marginal_likelihood`, so maximizing it over the hyperparameters
#' cannot overfit the way FITC can: FITC's diagonal correction
#' \eqn{\mathrm{diag}(K_{ff} - Q_{ff})} can absorb noise, and FITC tends to
#' underestimate the noise variance (Bauer, van der Wilk, and Rasmussen,
#' 2016). VFE predicts with the deterministic training conditional, which
#' has the same form as the FITC prediction with \eqn{\Lambda = \Sigma}.
#'
#' Both models report the trace term
#' \eqn{t = \sum_i (k(x_i, x_i) - [Q_{ff}]_{ii}) / \sigma_i^2} in
#' `trace_term`: the prior variance that the inducing points do not explain,
#' in units of noise variance. A trace term much larger than 1 means that
#' the inducing points do not cover the data well: VFE's bound is loose by at
#' least \eqn{t/2} and FITC leans on its diagonal correction. More or better
#' placed inducing points reduce it.
#'
#' @references Titsias, M. (2009). Variational learning of inducing variables
#'   in sparse Gaussian processes. *Proceedings of the 12th International
#'   Conference on Artificial Intelligence and Statistics*, 567--574.
#'
#'   Bauer, M., van der Wilk, M., and Rasmussen, C. E. (2016). Understanding
#'   probabilistic sparse Gaussian process approximations. *Advances in
#'   Neural Information Processing Systems*, 29.
#'
#' @return An object of class `gaussianprocesses_sparse_model`.
#' @seealso [gaussianprocesses_model] for the structure, versioning, and
#'   persistence of fitted models.
#' @examples
#' x <- seq(-3, 3, length.out = 200)
#' y <- sin(x) + 0.1 * cos(7 * x)
#'
#' model <- fit_sparse_gp(
#'   x,
#'   y,
#'   kernel = rbf_kernel(length_scale = 0.8),
#'   noise_variance = 0.01,
#'   n_inducing = 15
#' )
#' model
#' @template experimental
#' @export
fit_sparse_gp <- function(
    x,
    y,
    kernel,
    noise_variance = 1e-6,
    mean = zero_mean(),
    n_inducing = NULL,
    inducing_points = NULL,
    selection = c("farthest", "quantile", "variance"),
    method = c("fitc", "vfe"),
    fitc_diagonal_floor = 1e-10,
    initial_jitter = 0,
    fallback_jitter = 1e-10,
    jitter_multiplier = 10,
    max_attempts = 8L,
    symmetry_tolerance = sqrt(.Machine$double.eps)) {
  data <- .validate_gp_training_data(x, y)
  .assert_kernel(kernel)
  .assert_kernel_dimensions(kernel, ncol(data$x))
  .assert_mean_function(mean)
  .assert_fixed_mean(mean, "fit_sparse_gp()")
  mean <- .resolve_mean(mean, data$x)
  selection <- match.arg(selection)
  method <- match.arg(method)

  .assert_numeric_scalar(
    fitc_diagonal_floor,
    "fitc_diagonal_floor",
    lower = 0,
    lower_open = TRUE
  )

  training_noise <- .resolve_training_noise_variance(
    noise_variance,
    n_observations = nrow(data$x)
  )
  .assert_vfe_noise(method, training_noise)

  inducing <- .sparse_inducing_points(data$x, n_inducing, inducing_points,
                                      selection, kernel)

  model <- .fit_sparse_model(
    data$x,
    data$y,
    kernel,
    noise_variance,
    training_noise,
    mean,
    inducing$points,
    inducing$indices,
    inducing$method,
    method,
    fitc_diagonal_floor,
    list(
      initial_jitter = initial_jitter,
      fallback_jitter = fallback_jitter,
      jitter_multiplier = jitter_multiplier,
      max_attempts = max_attempts,
      symmetry_tolerance = symmetry_tolerance
    )
  )

  .warn_inducing_collapse(model)
}


#' Predict from a sparse Gaussian-process model
#'
#' Computes the posterior of the latent function and of new observations
#' under a FITC or VFE model ([fit_sparse_gp()]).
#'
#' @param model A fitted sparse GP model.
#' @param x Numeric vector or matrix of prediction inputs.
#' @param interval_level Probability level for intervals.
#' @param include_covariance If `TRUE`, construct the full latent covariance
#'   matrix. Otherwise only marginal variances are computed.
#' @param variance_tolerance Relative tolerance for posterior variances.
#' @param observation_noise_variance Optional prediction-time observation-noise
#'   variance. Required when the sparse model was fitted with observation-
#'   specific training noise.
#' @param block_size Number of prediction inputs processed together when
#'   `include_covariance = FALSE`. `NULL` chooses the largest block whose
#'   inducing-by-block cross-covariance holds at most \eqn{2^{20}} values, as
#'   in [predict_gp()]; memory then grows linearly with the number of
#'   prediction inputs. Results do not depend on it.
#'
#' @details
#' With \eqn{A = K_{uu} + K_{uf} \Lambda^{-1} K_{fu}}, where
#' \eqn{\Lambda} is the FITC diagonal or, for VFE, the noise, the latent
#' posterior at new inputs is
#' \deqn{\mu_* = m(x_*) + K_{*u} A^{-1} K_{uf} \Lambda^{-1} (y - m(X)),
#' \qquad \Sigma_{**} = K_{**} - Q_{**} + K_{*u} A^{-1} K_{u*},}
#' computed through the Cholesky factors of \eqn{K_{uu}} and
#' \eqn{B = I + V \Lambda^{-1} V^\top}. For VFE these are the predictive
#' equations of the deterministic training conditional.
#'
#' @return An object of class `gaussianprocesses_sparse_prediction`, with the
#'   shared prediction fields of [predict_gp()] and `approximation`.
#' @examples
#' x <- seq(-3, 3, length.out = 200)
#' y <- sin(x) + 0.1 * cos(7 * x)
#' kernel <- rbf_kernel(length_scale = 0.8)
#'
#' sparse <- fit_sparse_gp(x, y, kernel, noise_variance = 0.01, n_inducing = 15)
#' exact <- fit_gp(x, y, kernel, noise_variance = 0.01)
#'
#' x_new <- c(-1, 0, 1)
#' cbind(
#'   exact = predict_gp(exact, x_new)$mean,
#'   sparse = predict_sparse_gp(sparse, x_new)$mean
#' )
#' @template experimental
#' @export
predict_sparse_gp <- function(
    model,
    x,
    interval_level = 0.95,
    include_covariance = FALSE,
    variance_tolerance = sqrt(.Machine$double.eps),
    observation_noise_variance = NULL,
    block_size = NULL) {
  model <- .prepare_sparse_gp_model(model)
  prediction_x <- .validate_prediction_inputs(model, x)
  .validate_interval_level(interval_level)

  if (!is.logical(include_covariance) ||
      length(include_covariance) != 1L ||
      is.na(include_covariance)) {
    .gp_abort(
      "'include_covariance' must be TRUE or FALSE.",
      "gaussianprocesses_validation_error"
    )
  }

  prediction_noise <- .resolve_prediction_noise_variance(
    model,
    observation_noise_variance,
    n_predictions = nrow(prediction_x)
  )

  n_predictions <- nrow(prediction_x)
  block_size <- .prediction_block_size(
    block_size,
    n_training = model$n_inducing,
    n_predictions = n_predictions
  )

  # Q_ss = W'W and K_su A^{-1} K_us = Z'Z, with W = R_uu^{-T} K_us and
  # Z = R_B^{-T} W, matching the whitened factorization used in fitting.
  # Marginal quantities are accumulated in blocks of prediction inputs.
  whiten <- function(rows) {
    kus <- evaluate_kernel(
      model$kernel,
      model$inducing_points,
      prediction_x[rows, , drop = FALSE]
    )
    whitened <- .triangular_whiten(model$kuu_cholesky, kus)

    list(
      kus = kus,
      whitened = whitened,
      posterior = .triangular_whiten(model$b_cholesky, whitened)
    )
  }

  mean_offset <- numeric(n_predictions)
  qss_diagonal <- numeric(n_predictions)
  posterior_inducing_diagonal <- numeric(n_predictions)

  for (start in seq.int(1L, n_predictions, by = block_size)) {
    rows <- start:min(start + block_size - 1L, n_predictions)
    block <- whiten(rows)
    mean_offset[rows] <- drop(crossprod(block$kus, model$beta))
    qss_diagonal[rows] <- colSums(block$whitened^2)
    posterior_inducing_diagonal[rows] <- colSums(block$posterior^2)
  }

  posterior_mean <- evaluate_mean(model$mean, prediction_x) + mean_offset
  prior_diagonal <- kernel_diagonal(
    model$kernel,
    prediction_x
  )

  latent_variance <- .sanitize_posterior_variance(
    prior_diagonal -
      qss_diagonal +
      posterior_inducing_diagonal,
    reference_scale = max(abs(prior_diagonal)),
    tolerance = variance_tolerance
  )

  observation_variance <- latent_variance +
    prediction_noise
  latent_sd <- sqrt(latent_variance)
  observation_sd <- sqrt(observation_variance)
  quantile <- stats::qnorm((1 + interval_level) / 2)

  latent_interval <- cbind(
    lower = posterior_mean - quantile * latent_sd,
    upper = posterior_mean + quantile * latent_sd
  )
  prediction_interval <- cbind(
    lower = posterior_mean - quantile * observation_sd,
    upper = posterior_mean + quantile * observation_sd
  )

  result <- list(
    x = prediction_x,
    mean = posterior_mean,
    latent_variance = latent_variance,
    observation_noise_variance = prediction_noise,
    observation_variance = observation_variance,
    latent_sd = latent_sd,
    observation_sd = observation_sd,
    interval_level = interval_level,
    latent_interval = latent_interval,
    prediction_interval = prediction_interval,
    approximation = model$approximation
  )

  if (include_covariance) {
    prior_covariance <- evaluate_kernel(
      model$kernel,
      prediction_x
    )
    block <- whiten(seq_len(n_predictions))
    qss <- crossprod(block$whitened)
    posterior_inducing <- crossprod(block$posterior)

    latent_covariance <- prior_covariance -
      qss +
      posterior_inducing
    latent_covariance <- 0.5 * (
      latent_covariance + t(latent_covariance)
    )
    diag(latent_covariance) <- latent_variance

    result$latent_covariance <- latent_covariance
    result$observation_covariance <-
      latent_covariance +
      diag(prediction_noise, nrow(prediction_x))
  }

  structure(
    result,
    class = "gaussianprocesses_sparse_prediction"
  )
}


#' Benchmark sparse FITC against exact GP regression
#'
#' Runs both models with identical fixed kernel and noise parameters, then
#' compares retained model size, elapsed fitting/prediction time, posterior mean
#' RMSE, and latent-variance MAE. The exact GP is treated as the numerical
#' reference, not as ground truth.
#'
#' @param x Numeric vector or matrix of training inputs.
#' @param y Numeric response vector.
#' @param kernel Fixed kernel specification shared by both models.
#' @param n_inducing Number of inducing points for FITC.
#' @param prediction_x Optional prediction inputs; defaults to `x`.
#' @param noise_variance Non-negative observation-noise variance: a scalar,
#'   or one value per training observation. Both models use the same noise.
#' @param observation_noise_variance Optional noise variance at
#'   `prediction_x`, passed to both predictions. It is required when
#'   `noise_variance` has one value per observation, as in [predict_gp()].
#' @param mean Mean specification shared by both models.
#' @param selection Inducing-point selection method.
#'
#' @details
#' The benchmark reports retained R object size, not peak process memory.
#' Elapsed timings are simple wall-clock measurements intended for local
#' comparison, not publication-grade benchmarking.
#'
#' @return A list containing a comparison data frame and both prediction
#'   objects.
#' @examples
#' x <- seq(-3, 3, length.out = 200)
#' y <- sin(x) + 0.1 * cos(7 * x)
#'
#' benchmark <- benchmark_sparse_gp(
#'   x,
#'   y,
#'   kernel = rbf_kernel(length_scale = 0.8),
#'   n_inducing = 15,
#'   noise_variance = 0.01
#' )
#' benchmark$comparison
#' @template experimental
#' @export
benchmark_sparse_gp <- function(
    x,
    y,
    kernel,
    n_inducing,
    prediction_x = x,
    noise_variance = 1e-6,
    mean = zero_mean(),
    selection = c("farthest", "quantile", "variance"),
    observation_noise_variance = NULL) {
  selection <- match.arg(selection)
  .resolve_training_noise_variance(
    noise_variance,
    n_observations = NROW(x)
  )

  exact_fit_time <- system.time({
    exact_model <- fit_gp(
      x,
      y,
      kernel = kernel,
      noise_variance = noise_variance,
      mean = mean
    )
  })[["elapsed"]]

  sparse_fit_time <- system.time({
    sparse_model <- fit_sparse_gp(
      x,
      y,
      kernel = kernel,
      noise_variance = noise_variance,
      mean = mean,
      n_inducing = n_inducing,
      selection = selection
    )
  })[["elapsed"]]

  exact_predict_time <- system.time({
    exact_prediction <- predict_gp(
      exact_model,
      prediction_x,
      observation_noise_variance = observation_noise_variance
    )
  })[["elapsed"]]

  sparse_predict_time <- system.time({
    sparse_prediction <- predict_sparse_gp(
      sparse_model,
      prediction_x,
      observation_noise_variance = observation_noise_variance
    )
  })[["elapsed"]]

  mean_rmse <- sqrt(
    mean(
      (sparse_prediction$mean -
        exact_prediction$mean)^2
    )
  )
  variance_mae <- mean(
    abs(
      sparse_prediction$latent_variance -
        exact_prediction$latent_variance
    )
  )

  comparison <- data.frame(
    model = c("exact", "FITC"),
    n_training = c(
      exact_model$n_observations,
      sparse_model$n_observations
    ),
    n_inducing = c(
      exact_model$n_observations,
      sparse_model$n_inducing
    ),
    fit_elapsed_seconds = c(
      exact_fit_time,
      sparse_fit_time
    ),
    predict_elapsed_seconds = c(
      exact_predict_time,
      sparse_predict_time
    ),
    retained_model_bytes = c(
      as.numeric(utils::object.size(exact_model)),
      as.numeric(utils::object.size(sparse_model))
    ),
    posterior_mean_rmse_vs_exact = c(0, mean_rmse),
    latent_variance_mae_vs_exact = c(0, variance_mae),
    check.names = FALSE
  )

  list(
    comparison = comparison,
    exact_prediction = exact_prediction,
    sparse_prediction = sparse_prediction
  )
}


#' @export
print.gaussianprocesses_sparse_model <- function(x, ...) {
  x <- .prepare_sparse_gp_model(x, "x")
  vfe <- identical(.sparse_method(x), "vfe")
  cat(
    "Sparse Gaussian-process regression model (", x$approximation, ")\n",
    sep = ""
  )
  cat("  observations: ", x$n_observations, "\n", sep = "")
  cat("  inducing points: ", x$n_inducing, "\n", sep = "")
  cat(
    "  inducing selection: ",
    x$selection_method,
    "\n",
    sep = ""
  )

  if (!vfe) {
    cat(
      "  FITC floor applications: ",
      x$fitc_floor_count,
      "\n",
      sep = ""
    )
  }

  cat(
    "  trace term tr(K - Q) / noise: ",
    format(signif(x$trace_term, 4), trim = TRUE),
    "\n",
    sep = ""
  )
  cat(
    if (vfe) {
      "  variational lower bound on the log marginal likelihood: "
    } else {
      "  approximate log marginal likelihood: "
    },
    format(
      signif(x$log_marginal_likelihood, 7),
      trim = TRUE
    ),
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

  invisible(x)
}


#' @export
print.gaussianprocesses_sparse_prediction <- function(x, ...) {
  cat(
    "Sparse Gaussian-process prediction (",
    if (is.null(x$approximation)) "FITC" else x$approximation,
    ")\n",
    sep = ""
  )
  cat("  prediction points: ", length(x$mean), "\n", sep = "")
  cat(
    "  interval level: ",
    format(100 * x$interval_level, trim = TRUE),
    "%\n",
    sep = ""
  )
  invisible(x)
}
