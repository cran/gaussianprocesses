# Derivative observations ------------------------------------------------------
#
# Training observation i has a type t_i: 0 if it observes f(x_i), and d if it
# observes df(x_i) / dx_d. Differentiation is linear, so values and
# derivatives are jointly Gaussian with covariances
#   cov(f(x), f(y))                 = k(x, y),
#   cov(df(x) / dx_d, f(y))         = dk(x, y) / dx_d,
#   cov(f(x), df(y) / dy_e)         = dk(x, y) / dy_e,
#   cov(df(x) / dx_d, df(y) / dy_e) = d2k(x, y) / dx_d dy_e,
# which are the four parts of a derivative tuple from
# .kernel_input_derivatives() (R/kernel-input-derivatives.R). A typed
# covariance picks, for every pair of observations, the part their types
# select. Means contribute h(x)' beta for values and dh(x)' beta / dx_d for
# derivatives.
#
# The likelihood gradient needs d C / d eta for the typed C. Leaf kernels
# supply the hyperparameter derivatives of their input derivatives
# (input_parameter_gradient in the registry); sums, scaled kernels, and
# selections are linear, and products follow the product rule, which is
# bilinear in its two tuples.


.validate_derivative_observations <- function(derivative, x) {
  if (is.null(derivative)) {
    return(NULL)
  }

  if (!is.numeric(derivative) || !is.null(dim(derivative)) ||
      length(derivative) != nrow(x) || any(!is.finite(derivative)) ||
      any(derivative != round(derivative)) || any(derivative < 0) ||
      any(derivative > ncol(x))) {
    .gp_abort(
      sprintf(
        paste0(
          "'derivative' must give one integer per observation: 0 for a ",
          "function value, or d in 1..%d for the derivative with respect ",
          "to input column d."
        ),
        ncol(x)
      ),
      "gaussianprocesses_validation_error"
    )
  }

  if (all(derivative == 0)) {
    return(NULL)
  }

  as.integer(derivative)
}


# Labels of observation types, for diagnostics and printing.
.observation_type_labels <- function(types) {
  ifelse(types == 0L, "value", paste0("derivative x", types))
}


# The entries of a derivative tuple selected by the types of its rows and
# columns.
.select_typed <- function(tuple, x_types, y_types) {
  result <- tuple$value

  for (a in unique(x_types)) {
    rows <- x_types == a

    for (b in unique(y_types)) {
      if (a == 0L && b == 0L) {
        next
      }

      columns <- y_types == b
      block <- if (b == 0L) {
        tuple$gradient_x[[a]]
      } else if (a == 0L) {
        tuple$gradient_y[[b]]
      } else {
        tuple$cross[[a]][[b]]
      }
      result[rows, columns] <- block[rows, columns, drop = FALSE]
    }
  }

  result
}


# Covariance between typed points. `x_types` and `y_types` are integer
# vectors of types; untyped points have type 0.
.typed_covariance <- function(kernel, x, x_types, y, y_types) {
  derivative_rows <- any(x_types > 0L)
  derivative_columns <- any(y_types > 0L)

  if (!derivative_rows && !derivative_columns) {
    return(evaluate_kernel(kernel, x, y))
  }

  order <- if (derivative_rows && derivative_columns) 2L else 1L

  .select_typed(
    .kernel_input_derivatives(kernel, x, y, order = order),
    x_types,
    y_types
  )
}


# Prior variances of typed points.
.typed_diagonal <- function(kernel, x, types) {
  variance <- numeric(nrow(x))
  values <- types == 0L

  if (any(values)) {
    variance[values] <- kernel_diagonal(kernel, x[values, , drop = FALSE])
  }

  if (any(!values)) {
    derivative_variance <- .prior_gradient_variance(
      kernel,
      x[!values, , drop = FALSE]
    )
    variance[!values] <- derivative_variance[
      cbind(seq_len(sum(!values)), types[!values])
    ]
  }

  variance
}


# A derivative tuple and its hyperparameter derivatives ------------------------

.leaf_derivative_tuples <- function(kernel, x, y, prefix) {
  entry <- .kernel_entry(kernel$type)
  values <- .leaf_kernel_value_and_gradient(kernel, x, y, prefix)
  forward <- entry$input_parameter_gradient(x, y, kernel$parameters)
  backward <- entry$input_parameter_gradient(y, x, kernel$parameters)
  gradient <- list()

  for (name in names(kernel$parameters)) {
    n_values <- length(kernel$parameters[[name]])
    paths <- .parameter_paths(name, n_values, prefix)
    forward_terms <- forward[[name]]
    backward_terms <- backward[[name]]

    if (n_values == 1L) {
      forward_terms <- list(forward_terms)
      backward_terms <- list(backward_terms)
    }

    for (index in seq_along(paths)) {
      gradient[[paths[[index]]]] <- list(
        value = values$gradient[[paths[[index]]]],
        gradient_x = forward_terms[[index]]$gradient_x,
        gradient_y = lapply(backward_terms[[index]]$gradient_x, t),
        cross = forward_terms[[index]]$cross
      )
    }
  }

  list(
    tuple = .leaf_input_derivatives(kernel, x, y, 2L),
    gradient = gradient
  )
}


# The derivative tuple of a kernel and, named by parameter path, the
# derivatives of the tuple with respect to every optimizer coordinate.
.kernel_derivative_tuples <- function(kernel, x, y, prefix = "") {
  if (length(kernel$children) == 0L) {
    return(.leaf_derivative_tuples(kernel, x, y, prefix))
  }

  if (identical(kernel$type, "select")) {
    selected <- .kernel_derivative_tuples(
      kernel$children[[1L]],
      .select_input_columns(kernel, x),
      .select_input_columns(kernel, y),
      prefix
    )
    embed <- function(tuple) {
      .embed_input_derivatives(tuple, kernel$columns, ncol(x))
    }

    return(
      list(
        tuple = embed(selected$tuple),
        gradient = lapply(selected$gradient, embed)
      )
    )
  }

  if (identical(kernel$type, "changepoint")) {
    return(.changepoint_derivative_tuples(kernel, x, y, prefix))
  }

  children <- Map(
    function(child, child_prefix) {
      .kernel_derivative_tuples(child, x, y, child_prefix)
    },
    kernel$children,
    .child_prefixes(kernel, prefix)
  )
  tuples <- lapply(children, `[[`, "tuple")

  switch(
    kernel$type,
    sum = list(
      tuple = Reduce(.add_derivatives, tuples),
      gradient = do.call(c, lapply(children, `[[`, "gradient"))
    ),
    product = {
      gradient <- lapply(
        seq_along(children),
        function(index) {
          others <- Reduce(.multiply_derivatives, tuples[-index])
          lapply(
            children[[index]]$gradient,
            function(derivative) .multiply_derivatives(derivative, others)
          )
        }
      )

      list(
        tuple = Reduce(.multiply_derivatives, tuples),
        gradient = do.call(c, gradient)
      )
    },
    scale = {
      scale <- kernel$parameters$scale
      scaled <- function(tuple) {
        .map_derivatives(tuple, function(matrix) scale * matrix)
      }
      own <- list(scaled(tuples[[1L]]))
      names(own) <- .parameter_paths("scale", 1L, prefix)

      list(
        tuple = scaled(tuples[[1L]]),
        gradient = c(own, lapply(children[[1L]]$gradient, scaled))
      )
    },
    .gp_abort(
      sprintf("Unknown composite kernel type '%s'.", kernel$type),
      "gaussianprocesses_validation_error"
    )
  )
}


# d K / d eta for the typed training covariance, named by parameter path.
.typed_covariance_gradient <- function(kernel, x, types) {
  lapply(
    .kernel_derivative_tuples(kernel, x, x)$gradient,
    .select_typed,
    x_types = types,
    y_types = types
  )
}


# Typed means ------------------------------------------------------------------

# Basis rows h(x_i)' for values and dh(x_i)' / dx_d for derivatives.
.typed_mean_basis <- function(mean, x, types) {
  basis <- .mean_basis(mean, x)

  if (is.null(types)) {
    return(basis)
  }

  gradients <- .mean_basis_gradient(mean, x)

  for (d in setdiff(unique(types), 0L)) {
    rows <- types == d
    basis[rows, ] <- gradients[[d]][rows, , drop = FALSE]
  }

  basis
}


# Values of a mean with known coefficients at typed points.
.typed_mean_values <- function(mean, x, types) {
  if (is.null(types)) {
    return(evaluate_mean(mean, x))
  }

  treatment <- .mean_treatment(mean)

  if (identical(treatment$treatment, "none")) {
    return(rep(0, nrow(x)))
  }

  basis <- .typed_mean_basis(mean, x, types)
  drop(basis %*% .checked_coefficient_values(treatment$values, basis))
}


# Covariances of fitted models -------------------------------------------------

# The kernel part of the training covariance, K or its typed version.
.model_typed_kernel_covariance <- function(model) {
  if (is.null(model$derivative)) {
    return(evaluate_kernel(model$kernel, model$x))
  }

  .typed_covariance(
    model$kernel,
    model$x,
    model$derivative,
    model$x,
    model$derivative
  )
}


# Covariance between the training observations (rows) and function values
# at new inputs (columns).
.model_cross_covariance <- function(model, x) {
  if (is.null(model$derivative)) {
    return(evaluate_kernel(model$kernel, model$x, x))
  }

  .typed_covariance(
    model$kernel,
    model$x,
    model$derivative,
    x,
    integer(nrow(x))
  )
}


# Covariances between df(x) / dx_d at new inputs (rows) and the training
# observations (columns), one matrix per input dimension d.
.model_gradient_cross_covariance <- function(model, x) {
  types <- model$derivative

  if (is.null(types)) {
    return(
      .kernel_input_derivatives(model$kernel, x, model$x, order = 1L)$gradient_x
    )
  }

  derivatives <- .kernel_input_derivatives(model$kernel, x, model$x, 2L)

  lapply(
    seq_len(ncol(x)),
    function(d) .select_typed(derivatives, rep(d, nrow(x)), types)
  )
}


# d K / d eta for the training covariance of a model.
.model_covariance_gradient <- function(model) {
  if (is.null(model$derivative)) {
    return(.kernel_value_and_gradient(model$kernel, model$x, model$x)$gradient)
  }

  .typed_covariance_gradient(model$kernel, model$x, model$derivative)
}
