# Covariance derivatives with respect to the inputs ---------------------------
#
# For a kernel k(x, y), registry entries supply
#   input_gradient(x, y, parameters): list over d of dk(x_i, y_j) / dx_d
#   input_hessian(x, y, parameters):  list over d of lists over e of
#                                     d2k(x_i, y_j) / dx_d dy_e
# and `smoothness`, the number of mean-square derivatives of the process.
# These are the covariances of derivative processes: cov(df(x)/dx_d, f(y)) is
# dk / dx_d and cov(df(x)/dx_d, df(y)/dy_e) is d2k / dx_d dy_e, so they exist
# only for kernels with smoothness at least 1.
#
# For a radial kernel k = kappa(r) with r^2 = sum_d (x_d - y_d)^2 / l_d^2 and
# u_d = (x_d - y_d) / l_d^2,
#   dk / dx_d          = a u_d,                       a = kappa'(r) / r,
#   d2k / dx_d dy_e    = -b u_d u_e - a delta_de / l_d^2,
#                                    b = (kappa''(r) - kappa'(r) / r) / r^2.
# Each kernel supplies a and b in closed form, so r = 0 needs no limit. The
# derivations are in vignette("v07-kernel-derivatives").
#
# Derivative observations (R/derivative-observations.R) also need the
# hyperparameter derivatives of these matrices. Entries supply them as
#   input_parameter_gradient(x, y, parameters): a named list with, for every
#     parameter, list(gradient_x, cross): the derivatives of the two lists
#     above with respect to the parameter's optimizer coordinate. An ARD
#     parameter gives one such list per input dimension, or one list when a
#     single value is shared, as `gradient` does.


.radial_input_derivatives <- function(x, y, parameters, coefficients, order) {
  length_scale <- .resolve_length_scale(parameters$length_scale, ncol(x))
  differences <- lapply(
    seq_len(ncol(x)),
    function(dimension) outer(x[, dimension], y[, dimension], FUN = "-")
  )
  scaled <- Map(
    function(difference, scale) difference / scale^2,
    differences,
    length_scale
  )
  squared_distance <- Reduce(
    `+`,
    Map(function(difference, u) difference * u, differences, scaled)
  )
  coefficient <- coefficients(squared_distance, parameters)

  if (order == 1L) {
    return(lapply(scaled, function(u) coefficient$a * u))
  }

  lapply(
    seq_along(scaled),
    function(d) {
      lapply(
        seq_along(scaled),
        function(e) {
          cross <- -coefficient$b * scaled[[d]] * scaled[[e]]

          if (d == e) {
            cross <- cross - coefficient$a / length_scale[[d]]^2
          }

          cross
        }
      )
    }
  )
}


# input_gradient and input_hessian fields of a radial kernel's registry entry.
.radial_input_gradient <- function(coefficients) {
  function(x, y, parameters) {
    .radial_input_derivatives(x, y, parameters, coefficients, 1L)
  }
}


.radial_input_hessian <- function(coefficients) {
  function(x, y, parameters) {
    .radial_input_derivatives(x, y, parameters, coefficients, 2L)
  }
}


# input_parameter_gradient fields of radial entries.
.radial_input_parameter_gradient <- function(coefficients) {
  function(x, y, parameters) {
    .radial_input_parameter_derivatives(x, y, parameters, coefficients)
  }
}


# a, b, and c = (db / dr) / r for each radial kernel, as functions of r^2.
# Since da / d(r^2) = b / 2 and db / d(r^2) = c / 2, these give the
# length-scale derivatives. Kernels with further parameters add
# `parameters`, the derivatives of a and b in their optimizer coordinates.
.rbf_input_coefficients <- function(squared_distance, parameters) {
  covariance <- parameters$variance * exp(-0.5 * squared_distance)
  list(a = -covariance, b = covariance, c = -covariance)
}


.matern32_input_coefficients <- function(squared_distance, parameters) {
  distance <- sqrt(squared_distance)
  decay <- parameters$variance * exp(-sqrt(3) * distance)
  # b = 3 sqrt(3) variance exp(-sqrt(3) r) / r multiplies u_d u_e = O(r^2),
  # so its product is continuous and zero at r = 0. Likewise c = O(r^-3)
  # only multiplies s_c u_d u_e = O(r^4).
  b <- 3 * sqrt(3) * decay / distance
  c <- -3 * sqrt(3) * decay * (1 + sqrt(3) * distance) / distance^3
  b[distance == 0] <- 0
  c[distance == 0] <- 0

  list(a = -3 * decay, b = b, c = c)
}


.matern52_input_coefficients <- function(squared_distance, parameters) {
  distance <- sqrt(squared_distance)
  decay <- parameters$variance * exp(-sqrt(5) * distance)
  # c = O(1 / r) only multiplies s_c u_d u_e = O(r^4).
  c <- -(25 * sqrt(5) / 3) * decay / distance
  c[distance == 0] <- 0

  list(
    a = -(5 / 3) * (1 + sqrt(5) * distance) * decay,
    b = (25 / 3) * decay,
    c = c
  )
}


# With w = 1 + r^2 / (2 alpha), a = -variance w^(-alpha - 1) and
# b = variance ((alpha + 1) / alpha) w^(-alpha - 2), so
#   d a / d log(alpha) = a (-alpha log w + (alpha + 1) r^2 / (2 alpha w)),
#   d b / d log(alpha) = b (-1 / (alpha + 1) - alpha log w
#                           + (alpha + 2) r^2 / (2 alpha w)).
.rational_quadratic_input_coefficients <- function(
    squared_distance,
    parameters) {
  alpha <- parameters$alpha
  base <- 1 + squared_distance / (2 * alpha)
  log_base <- log1p(squared_distance / (2 * alpha))
  a <- -parameters$variance * base^(-alpha - 1)
  b <- parameters$variance * ((alpha + 1) / alpha) * base^(-alpha - 2)

  list(
    a = a,
    b = b,
    c = -parameters$variance * (alpha + 1) * (alpha + 2) / alpha^2 *
      base^(-alpha - 3),
    parameters = list(
      alpha = list(
        a = a * (-alpha * log_base +
          (alpha + 1) * squared_distance / (2 * alpha * base)),
        b = b * (-1 / (alpha + 1) - alpha * log_base +
          (alpha + 2) * squared_distance / (2 * alpha * base))
      )
    )
  )
}


# Elementwise sum of input-derivative lists list(gradient_x, cross).
.sum_input_derivative_terms <- function(terms) {
  Reduce(
    function(first, second) {
      list(
        gradient_x = Map(`+`, first$gradient_x, second$gradient_x),
        cross = Map(
          function(a, b) Map(`+`, a, b),
          first$cross,
          second$cross
        )
      )
    },
    terms
  )
}


# Hyperparameter derivatives of the input derivatives of a radial kernel.
# With s_c = (x_c - y_c)^2 / l_c^2,
#   d(a u_d) / d log l_c = -s_c b u_d - 2 a u_d delta_dc,
#   d(-b u_d u_e - a delta_de / l_d^2) / d log l_c
#     = s_c c u_d u_e + 2 b u_d u_e (delta_dc + delta_ec)
#       + s_c b delta_de / l_d^2 + 2 a delta_de delta_dc / l_d^2.
# A shared length scale differentiates as the sum over c. Every term is
# proportional to the variance, so its log-derivative is the term itself.
.radial_input_parameter_derivatives <- function(
    x,
    y,
    parameters,
    coefficients) {
  n_features <- ncol(x)
  dimensions <- seq_len(n_features)
  length_scale <- .resolve_length_scale(parameters$length_scale, n_features)
  differences <- lapply(
    dimensions,
    function(dimension) outer(x[, dimension], y[, dimension], FUN = "-")
  )
  scaled <- Map(
    function(difference, scale) difference / scale^2,
    differences,
    length_scale
  )
  components <- Map(`*`, differences, scaled)
  coefficient <- coefficients(Reduce(`+`, components), parameters)
  a <- coefficient$a
  b <- coefficient$b

  # gradient_x and cross for coefficients a and b.
  input_terms <- function(a, b) {
    list(
      gradient_x = lapply(scaled, function(u) a * u),
      cross = lapply(
        dimensions,
        function(d) {
          lapply(
            dimensions,
            function(e) {
              term <- -b * scaled[[d]] * scaled[[e]]

              if (d == e) {
                term <- term - a / length_scale[[d]]^2
              }

              term
            }
          )
        }
      )
    )
  }

  length_scale_terms <- lapply(
    dimensions,
    function(c_index) {
      s <- components[[c_index]]
      s_c <- s * coefficient$c
      s_b <- s * b

      list(
        gradient_x = lapply(
          dimensions,
          function(d) {
            term <- -s_b * scaled[[d]]

            if (d == c_index) {
              term <- term - 2 * a * scaled[[d]]
            }

            term
          }
        ),
        cross = lapply(
          dimensions,
          function(d) {
            lapply(
              dimensions,
              function(e) {
                factor <- s_c + 2 * b * ((d == c_index) + (e == c_index))
                term <- factor * scaled[[d]] * scaled[[e]]

                if (d == e) {
                  term <- term + s_b / length_scale[[d]]^2

                  if (d == c_index) {
                    term <- term + 2 * a / length_scale[[d]]^2
                  }
                }

                term
              }
            )
          }
        )
      )
    }
  )

  derivatives <- list(
    variance = input_terms(a, b),
    length_scale = if (length(parameters$length_scale) == 1L) {
      .sum_input_derivative_terms(length_scale_terms)
    } else {
      length_scale_terms
    }
  )

  for (name in names(coefficient$parameters)) {
    extra <- coefficient$parameters[[name]]
    derivatives[[name]] <- input_terms(extra$a, extra$b)
  }

  derivatives
}


# k = variance * exp(-2 sum_d sin^2(a_d) / l_d^2) with a_d = pi (x_d - y_d) / p.
# With g_d = -(2 pi / (p l_d^2)) sin(2 a_d),
#   dk / dx_d       = k g_d,
#   d2k / dx_d dy_e = -k g_d g_e + delta_de k (4 pi^2 / (p^2 l_d^2)) cos(2 a_d).
.periodic_input_derivatives <- function(x, y, parameters, order) {
  period <- parameters$period
  length_scale <- .resolve_length_scale(parameters$length_scale, ncol(x))
  angles <- lapply(
    seq_len(ncol(x)),
    function(dimension) {
      pi * outer(x[, dimension], y[, dimension], FUN = "-") / period
    }
  )
  covariance <- parameters$variance * exp(
    -2 * Reduce(
      `+`,
      Map(function(angle, scale) sin(angle)^2 / scale^2, angles, length_scale)
    )
  )
  factors <- Map(
    function(angle, scale) -(2 * pi / (period * scale^2)) * sin(2 * angle),
    angles,
    length_scale
  )

  if (order == 1L) {
    return(lapply(factors, function(factor) covariance * factor))
  }

  lapply(
    seq_along(factors),
    function(d) {
      lapply(
        seq_along(factors),
        function(e) {
          cross <- -covariance * factors[[d]] * factors[[e]]

          if (d == e) {
            cross <- cross + covariance *
              (4 * pi^2 / (period^2 * length_scale[[d]]^2)) *
              cos(2 * angles[[d]])
          }

          cross
        }
      )
    }
  )
}


# Hyperparameter derivatives of the periodic input derivatives. With
# q_c = sin^2(a_c) / l_c^2 and h_d = (4 pi^2 / (p^2 l_d^2)) cos(2 a_d),
#   dk / d log l_c = 4 k q_c,  dg_d / d log l_c = -2 g_d delta_dc,
#   dh_d / d log l_c = -2 h_d delta_dc,
#   dk / d log p = k P,  P = 2 sum_d a_d sin(2 a_d) / l_d^2,
#   dg_d / d log p = (2 pi / (p l_d^2)) (sin(2 a_d) + 2 a_d cos(2 a_d)),
#   dh_d / d log p = (4 pi^2 / (p^2 l_d^2)) (2 a_d sin(2 a_d) - 2 cos(2 a_d)),
# and the product rule gives the derivatives of k g_d and
# -k g_d g_e + delta_de k h_d.
.periodic_input_parameter_derivatives <- function(x, y, parameters) {
  period <- parameters$period
  dimensions <- seq_len(ncol(x))
  length_scale <- .resolve_length_scale(parameters$length_scale, ncol(x))
  angles <- lapply(
    dimensions,
    function(dimension) {
      pi * outer(x[, dimension], y[, dimension], FUN = "-") / period
    }
  )
  covariance <- parameters$variance * exp(
    -2 * Reduce(
      `+`,
      Map(function(angle, scale) sin(angle)^2 / scale^2, angles, length_scale)
    )
  )
  g <- Map(
    function(angle, scale) -(2 * pi / (period * scale^2)) * sin(2 * angle),
    angles,
    length_scale
  )
  h <- Map(
    function(angle, scale) (4 * pi^2 / (period^2 * scale^2)) * cos(2 * angle),
    angles,
    length_scale
  )

  # Derivatives of gradient_x and cross from the derivatives dk, dg, and dh
  # of k, g, and h with respect to one parameter.
  input_terms <- function(dk, dg, dh) {
    list(
      gradient_x = Map(
        function(g_d, dg_d) dk * g_d + covariance * dg_d,
        g,
        dg
      ),
      cross = lapply(
        dimensions,
        function(d) {
          lapply(
            dimensions,
            function(e) {
              term <- -dk * g[[d]] * g[[e]] -
                covariance * (dg[[d]] * g[[e]] + g[[d]] * dg[[e]])

              if (d == e) {
                term <- term + dk * h[[d]] + covariance * dh[[d]]
              }

              term
            }
          )
        }
      )
    )
  }

  zero <- matrix(0, nrow(x), nrow(y))
  no_change <- rep(list(zero), ncol(x))

  length_scale_terms <- lapply(
    dimensions,
    function(c_index) {
      changed <- function(values) {
        result <- no_change
        result[[c_index]] <- -2 * values[[c_index]]
        result
      }

      input_terms(
        4 * covariance * sin(angles[[c_index]])^2 / length_scale[[c_index]]^2,
        changed(g),
        changed(h)
      )
    }
  )

  period_factor <- 2 * Reduce(
    `+`,
    Map(
      function(angle, scale) angle * sin(2 * angle) / scale^2,
      angles,
      length_scale
    )
  )

  list(
    variance = input_terms(covariance, no_change, no_change),
    length_scale = if (length(parameters$length_scale) == 1L) {
      .sum_input_derivative_terms(length_scale_terms)
    } else {
      length_scale_terms
    },
    period = input_terms(
      covariance * period_factor,
      Map(
        function(angle, scale) {
          (2 * pi / (period * scale^2)) *
            (sin(2 * angle) + 2 * angle * cos(2 * angle))
        },
        angles,
        length_scale
      ),
      Map(
        function(angle, scale) {
          (4 * pi^2 / (period^2 * scale^2)) *
            (2 * angle * sin(2 * angle) - 2 * cos(2 * angle))
        },
        angles,
        length_scale
      )
    )
  )
}


# k = variance * sum_d x_d y_d: dk / dx_d = variance * y_d and
# d2k / dx_d dy_e = variance * delta_de.
.linear_input_derivatives <- function(x, y, parameters, order) {
  ones <- rep(1, nrow(x))

  if (order == 1L) {
    return(
      lapply(
        seq_len(ncol(x)),
        function(d) parameters$variance * outer(ones, y[, d])
      )
    )
  }

  if (order == 3L) {
    # Hyperparameter derivatives: only the variance, which scales both.
    return(
      list(
        variance = list(
          gradient_x = .linear_input_derivatives(x, y, parameters, 1L),
          cross = .linear_input_derivatives(x, y, parameters, 2L)
        )
      )
    )
  }

  constant <- matrix(parameters$variance, nrow(x), nrow(y))
  zero <- matrix(0, nrow(x), nrow(y))

  dimensions <- seq_len(ncol(x))

  lapply(
    dimensions,
    function(d) lapply(dimensions, function(e) if (d == e) constant else zero)
  )
}


# Smoothness of a kernel tree: the number of mean-square derivatives of its
# process. Sums, products, and scaled kernels are as smooth as their least
# smooth component.
.kernel_smoothness <- function(kernel) {
  if (length(kernel$children) == 0L) {
    return(.kernel_entry(kernel$type)$smoothness)
  }

  min(vapply(kernel$children, .kernel_smoothness, numeric(1)))
}


# Leaf kernels that are not mean-square differentiable.
.rough_kernel_labels <- function(kernel) {
  if (length(kernel$children) == 0L) {
    entry <- .kernel_entry(kernel$type)
    return(if (entry$smoothness < 1) entry$label else character(0))
  }

  unique(unlist(lapply(kernel$children, .rough_kernel_labels)))
}


.assert_differentiable_kernel <- function(kernel) {
  if (.kernel_smoothness(kernel) < 1) {
    .gp_abort(
      sprintf(
        paste0(
          "The kernel is not mean-square differentiable, so its derivatives ",
          "with respect to the inputs do not exist: %s. Use smoother ",
          "components, and observation noise instead of a white-noise kernel."
        ),
        paste(.rough_kernel_labels(kernel), collapse = ", ")
      ),
      c(
        "gaussianprocesses_smoothness_error",
        "gaussianprocesses_validation_error"
      )
    )
  }

  invisible(kernel)
}


.leaf_input_derivatives <- function(kernel, x, y, order) {
  entry <- .kernel_entry(kernel$type)
  parameters <- kernel$parameters
  value <- .evaluate_leaf_kernel(kernel, x, y)

  derivatives <- list(
    value = value,
    gradient_x = entry$input_gradient(x, y, parameters),
    gradient_y = lapply(entry$input_gradient(y, x, parameters), t)
  )

  if (order == 2L) {
    derivatives$cross <- entry$input_hessian(x, y, parameters)
  }

  derivatives
}


.map_derivatives <- function(derivatives, transform) {
  derivatives$value <- transform(derivatives$value)
  derivatives$gradient_x <- lapply(derivatives$gradient_x, transform)
  derivatives$gradient_y <- lapply(derivatives$gradient_y, transform)

  if (!is.null(derivatives$cross)) {
    derivatives$cross <- lapply(
      derivatives$cross,
      function(row) lapply(row, transform)
    )
  }

  derivatives
}


.add_derivatives <- function(first, second) {
  sum_lists <- function(a, b) Map(`+`, a, b)

  result <- list(
    value = first$value + second$value,
    gradient_x = sum_lists(first$gradient_x, second$gradient_x),
    gradient_y = sum_lists(first$gradient_y, second$gradient_y)
  )

  if (!is.null(first$cross)) {
    result$cross <- Map(sum_lists, first$cross, second$cross)
  }

  result
}


# Product rule for k = k_a k_b:
#   d(k_a k_b) / dx_d = da/dx_d b + a db/dx_d,
#   d2(k_a k_b) / dx_d dy_e = d2a b + da/dx_d db/dy_e + da/dy_e db/dx_d + a d2b.
.multiply_derivatives <- function(first, second) {
  a <- first$value
  b <- second$value

  result <- list(
    value = a * b,
    gradient_x = Map(
      function(da, db) da * b + a * db,
      first$gradient_x,
      second$gradient_x
    ),
    gradient_y = Map(
      function(da, db) da * b + a * db,
      first$gradient_y,
      second$gradient_y
    )
  )

  if (!is.null(first$cross)) {
    dimensions <- seq_along(first$gradient_x)
    result$cross <- lapply(
      dimensions,
      function(d) {
        lapply(
          dimensions,
          function(e) {
            first$cross[[d]][[e]] * b +
              first$gradient_x[[d]] * second$gradient_y[[e]] +
              first$gradient_y[[e]] * second$gradient_x[[d]] +
              a * second$cross[[d]][[e]]
          }
        )
      }
    )
  }

  result
}


# The covariance K(x, y) and its input derivatives: gradient_x[[d]] is
# dK / dx_d, gradient_y[[e]] is dK / dy_e, and, for order 2, cross[[d]][[e]]
# is d2K / dx_d dy_e. `x` and `y` are input matrices with matching columns.
.kernel_input_derivatives <- function(kernel, x, y, order = 1L) {
  if (length(kernel$children) == 0L) {
    return(.leaf_input_derivatives(kernel, x, y, order))
  }

  if (identical(kernel$type, "select")) {
    derivatives <- .kernel_input_derivatives(
      kernel$children[[1L]],
      .select_input_columns(kernel, x),
      .select_input_columns(kernel, y),
      order
    )

    return(.embed_input_derivatives(derivatives, kernel$columns, ncol(x)))
  }

  if (identical(kernel$type, "changepoint")) {
    return(.changepoint_input_derivatives(kernel, x, y, order))
  }

  children <- lapply(
    kernel$children,
    .kernel_input_derivatives,
    x = x,
    y = y,
    order = order
  )

  switch(
    kernel$type,
    sum = Reduce(.add_derivatives, children),
    product = Reduce(.multiply_derivatives, children),
    scale = .map_derivatives(
      children[[1L]],
      function(matrix) kernel$parameters$scale * matrix
    ),
    .gp_abort(
      sprintf("Unknown composite kernel type '%s'.", kernel$type),
      "gaussianprocesses_validation_error"
    )
  )
}


#' Derivatives of a covariance matrix with respect to the inputs
#'
#' Computes \eqn{\partial k(x_i, y_j) / \partial x_{i,d}} for every input
#' dimension \eqn{d}: the covariance between the derivative of the process
#' at \eqn{x_i} and its value at \eqn{y_j}.
#'
#' @param kernel A Gaussian-process kernel specification.
#' @param x Numeric vector or matrix of inputs; the derivatives are taken with
#'   respect to these.
#' @param y Optional numeric vector or matrix. If `NULL`, `x` is used.
#'
#' @details
#' Input derivatives exist only for mean-square differentiable kernels. The
#' Matérn-1/2 and white-noise kernels, and any sum or product that contains
#' one, are not, and raise an error of class
#' `gaussianprocesses_smoothness_error`. The Matérn-3/2 process is
#' differentiable once, the Matérn-5/2 process twice, and the RBF,
#' rational-quadratic, periodic, and linear processes infinitely often.
#'
#' Derivatives are closed-form expressions for every kernel, written so that
#' coincident inputs need no limit; sums, products, and scaled kernels follow
#' the sum and product rules. [predict_gradient_gp()] uses them to predict
#' the gradient of the latent function. The derivations are in
#' `vignette("v07-kernel-derivatives", package = "gaussianprocesses")`.
#'
#' @return A list with one numeric matrix per input dimension. Element `d`
#'   has `nrow(x)` rows and `nrow(y)` columns and holds
#'   \eqn{\partial k(x_i, y_j) / \partial x_{i,d}}.
#' @seealso [kernel_gradient()] for derivatives with respect to the
#'   hyperparameters.
#' @examples
#' kernel <- rbf_kernel(variance = 1, length_scale = 0.5)
#' x <- c(0, 0.25, 1)
#'
#' gradient <- kernel_input_gradient(kernel, x)
#' gradient[[1]]
#'
#' # Compare with a central difference in the first input.
#' h <- 1e-6
#' (evaluate_kernel(kernel, x + h, x) -
#'   evaluate_kernel(kernel, x - h, x)) / (2 * h)
#' @template experimental
#' @export
kernel_input_gradient <- function(kernel, x, y = NULL) {
  .assert_kernel(kernel)
  .assert_differentiable_kernel(kernel)
  inputs <- .prepare_kernel_inputs(x, y)

  .kernel_input_derivatives(kernel, inputs$x, inputs$y, order = 1L)$gradient_x
}
