# Spectral-mixture kernels ----------------------------------------------------
#
# A spectral-mixture kernel (Wilson and Adams, 2013) is a sum of spectral
# components. A component with weight w, frequencies mu_d, and spectral
# variances v_d is
#   k(x, y) = w prod_d g_d(tau_d),  tau_d = x_d - y_d,
#   g_d(tau) = exp(-2 pi^2 v_d tau^2) cos(2 pi mu_d tau),
# the covariance of a Gaussian spectral density at +-mu with variance v,
# symmetrized. Each factor is the real part of
#   h_d(tau) = exp(phi_d),  phi_d = -2 pi^2 v_d tau^2 + 2 pi i mu_d tau,
# so with phi' = -4 pi^2 v tau + 2 pi i mu and phi'' = -4 pi^2 v,
#   g' = Re(phi' h),  g'' = Re((phi'' + phi'^2) h).
# Parameters enter through phi only, which keeps their derivatives short.
#
# k is multilinear in the triples (g_d, g'_d, g''_d): every entry of its
# derivative tuple contains exactly one member of each triple. A derivative
# with respect to a parameter of dimension c is therefore the same tuple
# with the triple of c replaced by its derivative.


# The triple (g, g', g'') of every input dimension and, with `derivatives =
# TRUE`, the derivatives of the triple with respect to log v_d and mu_d.
# With `inputs = FALSE` only g and its parameter derivatives are computed.
.spectral_factors <- function(
    x,
    y,
    parameters,
    derivatives = FALSE,
    inputs = TRUE) {
  n_features <- ncol(x)
  frequency <- .resolve_length_scale(parameters$frequency, n_features)
  spectral_variance <- .resolve_length_scale(
    parameters$spectral_variance,
    n_features
  )

  lapply(
    seq_len(n_features),
    function(d) {
      tau <- outer(x[, d], y[, d], FUN = "-")
      v <- spectral_variance[[d]]
      mu <- frequency[[d]]
      if (!inputs) {
        # Real arithmetic: g = e cos(2 pi mu tau), dg / d log v =
        # -2 pi^2 v tau^2 g, and dg / d mu = -2 pi tau e sin(2 pi mu tau).
        envelope <- exp(-2 * pi^2 * v * tau^2)
        angle <- 2 * pi * mu * tau
        factor <- list(value = envelope * cos(angle))

        if (derivatives) {
          factor$spectral_variance <- list(
            value = -2 * pi^2 * v * tau^2 * factor$value
          )
          factor$frequency <- list(
            value = -2 * pi * tau * envelope * sin(angle)
          )
        }

        return(factor)
      }

      h <- exp(-2 * pi^2 * v * tau^2 + 1i * (2 * pi * mu * tau))

      slope <- -4 * pi^2 * v * tau + 1i * (2 * pi * mu)
      curvature <- -4 * pi^2 * v + slope^2
      factor <- list(
        value = Re(h),
        first = Re(slope * h),
        second = Re(curvature * h)
      )

      if (!derivatives) {
        return(factor)
      }

      # d / d log v: d phi = -2 pi^2 v tau^2, d phi' = -4 pi^2 v tau,
      # d phi'' = -4 pi^2 v. d / d mu: d phi = 2 pi i tau, d phi' = 2 pi i,
      # d phi'' = 0.
      differentiate <- function(d_phase, d_slope, d_second) {
        d_h <- d_phase * h
        list(
          value = Re(d_h),
          first = Re(d_slope * h + slope * d_h),
          second = Re((d_second + 2 * slope * d_slope) * h + curvature * d_h)
        )
      }

      factor$spectral_variance <- differentiate(
        -2 * pi^2 * v * tau^2,
        -4 * pi^2 * v * tau,
        -4 * pi^2 * v
      )
      factor$frequency <- differentiate(
        1i * (2 * pi * tau),
        1i * (2 * pi),
        0
      )
      factor
    }
  )
}


# Product of the matrices of `values` at positions other than `skip`, or
# `ones` when there are none.
.product_except <- function(values, skip, ones) {
  kept <- values[-skip]

  if (length(kept) == 0L) {
    return(ones)
  }

  Reduce(`*`, kept)
}


# The derivative tuple (value, gradient_x, gradient_y, cross) of
# w prod_d g_d from the triples of every dimension.
.spectral_tuple <- function(weight, factors, order) {
  dimensions <- seq_along(factors)
  values <- lapply(factors, `[[`, "value")
  ones <- matrix(1, nrow(values[[1L]]), ncol(values[[1L]]))

  gradient_x <- lapply(
    dimensions,
    function(d) {
      weight * factors[[d]]$first * .product_except(values, d, ones)
    }
  )
  tuple <- list(
    value = weight * Reduce(`*`, values),
    gradient_x = gradient_x,
    # k depends on x - y, so d / dy = -d / dx.
    gradient_y = lapply(gradient_x, `-`)
  )

  if (order == 2L) {
    tuple$cross <- lapply(
      dimensions,
      function(d) {
        lapply(
          dimensions,
          function(e) {
            if (d == e) {
              return(
                -weight * factors[[d]]$second *
                  .product_except(values, d, ones)
              )
            }

            -weight * factors[[d]]$first * factors[[e]]$first *
              .product_except(values, c(d, e), ones)
          }
        )
      }
    )
  }

  tuple
}


# Derivatives of the tuple with respect to one parameter of each dimension,
# summed over dimensions when the parameter is shared.
.spectral_parameter_tuples <- function(parameters, factors, name, order) {
  per_dimension <- lapply(
    seq_along(factors),
    function(c_index) {
      changed <- factors
      changed[[c_index]] <- factors[[c_index]][[name]]
      .spectral_tuple(parameters$weight, changed, order)
    }
  )

  if (length(parameters[[name]]) == 1L) {
    return(Reduce(.add_derivatives, per_dimension))
  }

  per_dimension
}


# Hyperparameter derivatives of the covariance: w dg_c prod_{d != c} g_d for
# a parameter of dimension c, summed over c when it is shared.
.spectral_component_gradient <- function(x, y, parameters, covariance) {
  factors <- .spectral_factors(
    x,
    y,
    parameters,
    derivatives = TRUE,
    inputs = FALSE
  )
  values <- lapply(factors, `[[`, "value")
  ones <- matrix(1, nrow(x), nrow(y))

  per_parameter <- function(name) {
    per_dimension <- lapply(
      seq_along(factors),
      function(c_index) {
        parameters$weight * factors[[c_index]][[name]]$value *
          .product_except(values, c_index, ones)
      }
    )

    if (length(parameters[[name]]) == 1L) {
      return(Reduce(`+`, per_dimension))
    }

    per_dimension
  }

  list(
    weight = covariance,
    frequency = per_parameter("frequency"),
    spectral_variance = per_parameter("spectral_variance")
  )
}


.spectral_input_parameter_gradient <- function(x, y, parameters) {
  factors <- .spectral_factors(x, y, parameters, derivatives = TRUE)
  input_terms <- function(tuple) tuple[c("gradient_x", "cross")]
  per_parameter <- function(name) {
    tuples <- .spectral_parameter_tuples(parameters, factors, name, 2L)

    if (!is.null(tuples$value)) {
      return(input_terms(tuples))
    }

    lapply(tuples, input_terms)
  }

  list(
    weight = input_terms(.spectral_tuple(parameters$weight, factors, 2L)),
    frequency = per_parameter("frequency"),
    spectral_variance = per_parameter("spectral_variance")
  )
}


.spectral_component_entry <- function() {
  list(
    type = "spectral_component",
    label = "SpectralComponent",
    parameters = list(
      weight = .kernel_parameter(1),
      # Real-valued: cos is even, so -mu and mu give the same kernel, and
      # mu = 0 is an RBF component.
      frequency = .kernel_parameter(
        0,
        shape = "ard",
        constraint = .real_constraint()
      ),
      spectral_variance = .kernel_parameter(1, shape = "ard")
    ),
    stationary = TRUE,
    smoothness = Inf,
    evaluate = function(x, y, parameters) {
      factors <- .spectral_factors(x, y, parameters, inputs = FALSE)
      parameters$weight * Reduce(`*`, lapply(factors, `[[`, "value"))
    },
    diagonal = function(x, parameters) rep(parameters$weight, nrow(x)),
    gradient = .spectral_component_gradient,
    input_gradient = function(x, y, parameters) {
      .spectral_tuple(
        parameters$weight,
        .spectral_factors(x, y, parameters),
        1L
      )$gradient_x
    },
    input_hessian = function(x, y, parameters) {
      .spectral_tuple(
        parameters$weight,
        .spectral_factors(x, y, parameters),
        2L
      )$cross
    },
    input_parameter_gradient = .spectral_input_parameter_gradient,
    state_space = NULL
  )
}


# One spectral component; spectral_mixture_kernel() combines them.
.spectral_component_kernel <- function(
    weight = 1,
    frequency = 0,
    spectral_variance = 1) {
  .new_leaf_kernel(
    "spectral_component",
    list(
      weight = weight,
      frequency = frequency,
      spectral_variance = spectral_variance
    )
  )
}


# Values for each component from a scalar, a vector with one value per
# component, or a matrix with one row per component and one column per input
# dimension.
.component_rows <- function(values, n_components, name) {
  if (is.matrix(values)) {
    if (nrow(values) != n_components) {
      .gp_abort(
        sprintf(
          "'%s' must have one row per component (%d); it has %d.",
          name,
          n_components,
          nrow(values)
        ),
        "gaussianprocesses_validation_error"
      )
    }

    return(lapply(seq_len(n_components), function(q) values[q, ]))
  }

  if (!is.numeric(values) || !(length(values) %in% c(1L, n_components))) {
    .gp_abort(
      sprintf(
        paste0(
          "'%s' must be one value, one value per component (%d), or a ",
          "matrix with one row per component."
        ),
        name,
        n_components
      ),
      "gaussianprocesses_validation_error"
    )
  }

  as.list(rep_len(as.numeric(values), n_components))
}


#' Spectral-mixture kernels
#'
#' A spectral-mixture kernel (Wilson and Adams, 2013) is a sum of
#' \eqn{Q} components with Gaussian spectral densities. With enough
#' components it approximates any stationary covariance.
#'
#' @param n_components Number of components \eqn{Q}.
#' @param weights Positive component weights \eqn{w_q}: one value, recycled,
#'   or one per component.
#' @param frequencies Frequencies \eqn{\mu_q}, in cycles per input unit: one
#'   value, one per component, or a matrix with one row per component and one
#'   column per input dimension.
#' @param spectral_variances Positive spectral variances \eqn{v_q}, in
#'   squared cycles per input unit, shaped like `frequencies`.
#'
#' @details
#' By Bochner's theorem a stationary kernel is the Fourier transform of a
#' non-negative spectral density. The spectral mixture uses
#' \eqn{S(s) = \sum_q \tfrac{w_q}{2} [N(s; \mu_q, V_q) + N(s; -\mu_q, V_q)]},
#' whose transform is
#' \deqn{k(\tau) = \sum_{q=1}^{Q} w_q \prod_{d=1}^{D}
#' \exp(-2 \pi^2 \tau_d^2 v_{qd}) \cos(2 \pi \tau_d \mu_{qd}).}
#' So \eqn{k(0) = \sum_q w_q}. A component with \eqn{\mu_q = 0} is an RBF
#' kernel with variance \eqn{w_q} and length scale
#' \eqn{\ell = 1 / (2 \pi \sqrt{v_q})}; a component with a small spectral
#' variance is a long-lived oscillation with period \eqn{1 / \mu_q}.
#'
#' Each component is a kernel with parameters `weight`, `frequency`, and
#' `spectral_variance`, and the mixture is their [sum_kernel()], so parameter
#' paths are `kernel1.weight`, `kernel1.frequency`, and so on; a single
#' component has no prefix. Frequencies are real-valued and optimized on the
#' identity scale, the others on the log scale. Because \eqn{\cos} is even,
#' \eqn{-\mu} and \eqn{\mu} give the same kernel. [optimize_gp()] bounds
#' every frequency by 0 and the Nyquist frequency \eqn{1 / (2 \Delta)} of its
#' input column, with \eqn{\Delta} the smallest spacing between distinct
#' inputs, and displaces it in later starts by one frequency-resolution
#' bin, the reciprocal of the input range.
#'
#' The likelihood is strongly multimodal: components can swap, merge, or
#' settle on harmonics. Start from the data's spectrum with
#' [initialize_spectral_mixture()], and use several starts; the examples use
#' `n_starts = 5`.
#'
#' @return A kernel specification: one spectral component, or the sum of
#'   `n_components` of them.
#' @references Wilson, A. G. and Adams, R. P. (2013). Gaussian process
#'   kernels for pattern discovery and extrapolation. *Proceedings of the
#'   30th International Conference on Machine Learning*, 1067--1075.
#' @seealso [initialize_spectral_mixture()].
#' @examples
#' kernel <- spectral_mixture_kernel(
#'   n_components = 2,
#'   weights = c(1, 0.5),
#'   frequencies = c(0, 0.25),
#'   spectral_variances = c(0.01, 0.001)
#' )
#' kernel
#'
#' # k(0) is the sum of the weights.
#' evaluate_kernel(kernel, 0)
#'
#' # A component with frequency 0 is an RBF kernel.
#' all.equal(
#'   evaluate_kernel(spectral_mixture_kernel(1, 2, 0, 1 / (4 * pi^2)), 0:3),
#'   evaluate_kernel(rbf_kernel(variance = 2, length_scale = 1), 0:3)
#' )
#' @template experimental
#' @export
spectral_mixture_kernel <- function(
    n_components = 1L,
    weights = 1 / n_components,
    frequencies = 0,
    spectral_variances = 1) {
  if (!is.numeric(n_components) || length(n_components) != 1L ||
      !is.finite(n_components) || n_components < 1 ||
      n_components != round(n_components)) {
    .gp_abort(
      "'n_components' must be a positive integer.",
      "gaussianprocesses_validation_error"
    )
  }

  n_components <- as.integer(n_components)
  weights <- .component_rows(weights, n_components, "weights")
  frequencies <- .component_rows(frequencies, n_components, "frequencies")
  spectral_variances <- .component_rows(
    spectral_variances,
    n_components,
    "spectral_variances"
  )

  components <- lapply(
    seq_len(n_components),
    function(q) {
      .spectral_component_kernel(
        weight = weights[[q]],
        frequency = frequencies[[q]],
        spectral_variance = spectral_variances[[q]]
      )
    }
  )

  if (n_components == 1L) {
    return(components[[1L]])
  }

  do.call(sum_kernel, components)
}


# Frequency limits -------------------------------------------------------------

# Smallest spacing between distinct values, or NA with fewer than two.
.minimum_spacing <- function(values) {
  distinct <- sort(unique(values))

  if (length(distinct) < 2L) {
    return(NA_real_)
  }

  min(diff(distinct))
}


# Nyquist frequency 1 / (2 Delta) and frequency resolution 1 / range of each
# spectral frequency of a kernel, named by parameter path. `columns` maps the
# columns a node receives to columns of `x`.
.spectral_frequency_limits <- function(kernel, x, columns, prefix = "") {
  if (length(kernel$children) == 0L) {
    if (!identical(kernel$type, "spectral_component")) {
      return(list())
    }

    n_values <- length(kernel$parameters$frequency)
    paths <- .parameter_paths("frequency", n_values, prefix)
    column_sets <- if (n_values == 1L) {
      list(columns)
    } else {
      as.list(columns)
    }

    limits <- lapply(
      column_sets,
      function(set) {
        spacing <- vapply(set, function(j) .minimum_spacing(x[, j]), 1)
        span <- vapply(set, function(j) diff(range(x[, j])), 1)

        list(
          upper = if (all(is.na(spacing))) {
            Inf
          } else {
            1 / (2 * min(spacing, na.rm = TRUE))
          },
          resolution = if (all(span == 0)) 1 else 1 / max(span)
        )
      }
    )
    names(limits) <- paths
    return(limits)
  }

  if (identical(kernel$type, "select")) {
    columns <- columns[kernel$columns]
  }

  do.call(
    c,
    unname(
      Map(
        function(child, child_prefix) {
          .spectral_frequency_limits(child, x, columns, child_prefix)
        },
        kernel$children,
        .child_prefixes(kernel, prefix)
      )
    )
  )
}


# Initialisation from the empirical spectrum -----------------------------------

# Lomb-Scargle periodogram of irregularly sampled `y` at `frequencies`
# (cycles per input unit), with the time offset tau that makes the sine and
# cosine terms orthogonal.
.lomb_scargle <- function(x, y, frequencies) {
  vapply(
    frequencies,
    function(frequency) {
      omega <- 2 * pi * frequency
      tau <- atan2(sum(sin(2 * omega * x)), sum(cos(2 * omega * x))) /
        (2 * omega)
      cosine <- cos(omega * (x - tau))
      sine <- sin(omega * (x - tau))

      0.5 * (sum(y * cosine)^2 / sum(cosine^2) + sum(y * sine)^2 / sum(sine^2))
    },
    numeric(1)
  )
}


# Empirical spectrum of centred `y`: the periodogram for regular inputs and
# the Lomb-Scargle periodogram, on a grid four times finer than the
# resolution and up to half the reciprocal median spacing, otherwise.
.empirical_spectrum <- function(x, y) {
  spacing <- diff(x)
  span <- x[[length(x)]] - x[[1L]]

  if (all(abs(spacing - spacing[[1L]]) <= 1e-8 * abs(spacing[[1L]]))) {
    periodogram <- stats::spec.pgram(
      stats::ts(y, deltat = spacing[[1L]]),
      taper = 0,
      detrend = FALSE,
      demean = FALSE,
      fast = FALSE,
      plot = FALSE
    )

    return(
      list(
        frequency = periodogram$freq,
        power = periodogram$spec,
        resolution = 1 / (length(x) * spacing[[1L]]),
        sampling = "regular"
      )
    )
  }

  frequency <- seq(
    1 / span,
    1 / (2 * stats::median(spacing)),
    by = 1 / (4 * span)
  )

  list(
    frequency = frequency,
    power = .lomb_scargle(x, y, frequency),
    resolution = 1 / span,
    sampling = "irregular"
  )
}


# Weighted EM for a mixture of `n_components` Gaussians and a uniform
# background on the spectrum, which absorbs the flat spectrum of the noise.
# Components start at the largest local maxima of the power.
.fit_spectral_mixture <- function(spectrum, n_components) {
  frequency <- spectrum$frequency
  power <- spectrum$power
  resolution <- spectrum$resolution
  n_grid <- length(frequency)

  peaks <- which(
    power > c(-Inf, power[-n_grid]) & power >= c(power[-1L], -Inf)
  )
  peaks <- peaks[order(power[peaks], decreasing = TRUE)]
  # With fewer peaks than components, the strongest other frequencies.
  others <- setdiff(order(power, decreasing = TRUE), peaks)
  start <- c(peaks, others)[seq_len(n_components)]

  mean <- frequency[start]
  variance <- rep(resolution^2, n_components)
  variance_floor <- (resolution / 4)^2
  share <- c(0.5, 0.5 * power[start] / sum(power[start]))
  background <- 1 / (diff(range(frequency)) + resolution)
  weight <- power / sum(power)
  previous <- -Inf

  for (iteration in seq_len(1000L)) {
    density <- cbind(
      share[[1L]] * background,
      vapply(
        seq_len(n_components),
        function(q) {
          share[[q + 1L]] *
            stats::dnorm(frequency, mean[[q]], sqrt(variance[[q]]))
        },
        numeric(n_grid)
      )
    )
    total <- pmax(rowSums(density), .Machine$double.xmin)
    log_likelihood <- sum(weight * log(total))
    responsibility <- density / total
    mass <- colSums(weight * responsibility)
    share <- mass / sum(mass)

    component <- responsibility[, -1L, drop = FALSE]
    component_mass <- pmax(mass[-1L], .Machine$double.xmin)
    mean <- colSums(weight * component * frequency) / component_mass
    variance <- pmax(
      colSums(weight * component * outer(frequency, mean, "-")^2) /
        component_mass,
      variance_floor
    )

    if (abs(log_likelihood - previous) <= 1e-10 * abs(log_likelihood)) {
      break
    }

    previous <- log_likelihood
  }

  list(share = share[-1L], mean = mean, variance = variance)
}


#' Initialise a spectral-mixture kernel from the data's spectrum
#'
#' Builds a [spectral_mixture_kernel()] whose components sit at the main
#' peaks of the empirical spectrum of one-dimensional data, a starting point
#' for [optimize_gp()].
#'
#' @param x Numeric vector of inputs, regularly or irregularly spaced.
#' @param y Numeric vector of responses.
#' @param n_components Number of components \eqn{Q}.
#'
#' @details
#' The spectrum of the centred responses is the periodogram
#' ([stats::spec.pgram()], without taper) when the inputs are regularly
#' spaced, and the Lomb-Scargle periodogram otherwise, computed on a grid
#' four times finer than the frequency resolution \eqn{1 / T}, for an input
#' range \eqn{T}, up to half the reciprocal median spacing. A mixture of
#' \eqn{Q} Gaussians and a uniform background is fitted to the spectrum, as
#' a density over frequency weighted by the power, by EM, starting from the
#' \eqn{Q} largest local maxima. The background absorbs the flat spectrum of
#' observation noise, which would otherwise widen the components.
#'
#' Component \eqn{q} takes the mean and variance of Gaussian \eqn{q} as its
#' frequency and spectral variance, and the share of the variance of
#' \eqn{y} that the Gaussian explains as its weight; the background share is
#' left to the noise. Spectral variances are at least
#' \eqn{(1 / (4 T))^2}, the order of a spectral peak's width.
#'
#' The likelihood of a spectral mixture is strongly multimodal, and this
#' initialisation is only a starting point: optimize with several starts,
#' for example `n_starts = 5`, which displace each frequency by one
#' resolution bin in turn.
#'
#' @return A spectral-mixture kernel specification with `n_components`
#'   components.
#' @seealso [spectral_mixture_kernel()].
#' @examples
#' set.seed(1)
#' x <- 0:199
#' y <- 2 * sin(2 * pi * 0.1 * x) + sin(2 * pi * 0.27 * x) +
#'   rnorm(200, sd = 0.5)
#'
#' kernel <- initialize_spectral_mixture(x, y, n_components = 2)
#' kernel
#' @template experimental
#' @export
initialize_spectral_mixture <- function(x, y, n_components = 2L) {
  data <- .validate_gp_training_data(x, y)

  if (ncol(data$x) != 1L) {
    .gp_abort(
      "initialize_spectral_mixture() needs one input dimension.",
      "gaussianprocesses_validation_error"
    )
  }

  if (!is.numeric(n_components) || length(n_components) != 1L ||
      !is.finite(n_components) || n_components < 1 ||
      n_components != round(n_components)) {
    .gp_abort(
      "'n_components' must be a positive integer.",
      "gaussianprocesses_validation_error"
    )
  }

  order <- order(data$x[, 1L])
  inputs <- data$x[order, 1L]
  responses <- data$y[order] - mean(data$y)

  if (length(unique(inputs)) < 2L * n_components + 2L) {
    .gp_abort(
      sprintf(
        "A spectrum for %d component(s) needs at least %d distinct inputs.",
        n_components,
        2L * n_components + 2L
      ),
      "gaussianprocesses_validation_error"
    )
  }

  spectrum <- .empirical_spectrum(inputs, responses)
  mixture <- .fit_spectral_mixture(spectrum, as.integer(n_components))

  spectral_mixture_kernel(
    n_components = n_components,
    weights = pmax(mixture$share, 1e-8) * stats::var(responses),
    frequencies = mixture$mean,
    spectral_variances = mixture$variance
  )
}
