# Kernel utilities ------------------------------------------------------------

.as_input_matrix <- function(x, name = deparse(substitute(x))) {
  if (is.numeric(x) && is.null(dim(x))) {
    .assert_numeric_vector(x, name = name)
    return(matrix(x, ncol = 1L))
  }

  .assert_numeric_matrix(x, name = name)
  x
}


.prepare_kernel_inputs <- function(x, y = NULL) {
  x_matrix <- .as_input_matrix(x, name = "x")
  y_matrix <- if (is.null(y)) {
    x_matrix
  } else {
    .as_input_matrix(y, name = "y")
  }

  if (ncol(x_matrix) != ncol(y_matrix)) {
    .gp_abort(
      paste0(
        "'x' and 'y' must have the same number of columns. ",
        "Received ",
        ncol(x_matrix),
        " and ",
        ncol(y_matrix),
        "."
      ),
      "gaussianprocesses_validation_error"
    )
  }

  list(x = x_matrix, y = y_matrix)
}


.pairwise_squared_distance <- function(x, y) {
  # Accumulate squared coordinate differences rather than expanding
  # ||x||^2 + ||y||^2 - 2 x'y. The expansion cancels catastrophically for
  # inputs far from the origin and leaves round-off on coincident rows.
  squared <- matrix(0, nrow = nrow(x), ncol = nrow(y))

  for (dimension in seq_len(ncol(x))) {
    squared <- squared +
      outer(x[, dimension], y[, dimension], FUN = "-")^2
  }

  squared
}


.pairwise_identical_rows <- function(x, y) {
  identical_rows <- matrix(TRUE, nrow = nrow(x), ncol = nrow(y))

  for (dimension in seq_len(ncol(x))) {
    identical_rows <- identical_rows &
      outer(x[, dimension], y[, dimension], FUN = "==")
  }

  identical_rows
}


# One length scale per input dimension. Kernel parameters are validated by
# the kernel registry before they reach this helper.
.resolve_length_scale <- function(length_scale, n_features) {
  if (length(length_scale) == 1L) {
    return(rep(length_scale, n_features))
  }

  as.numeric(length_scale)
}


.scaled_pairwise_squared_distance <- function(x, y, length_scale) {
  resolved <- .resolve_length_scale(length_scale, ncol(x))
  scaled_x <- sweep(x, 2L, resolved, FUN = "/")
  scaled_y <- sweep(y, 2L, resolved, FUN = "/")

  .pairwise_squared_distance(scaled_x, scaled_y)
}


# Per-dimension squared differences of length-scale-scaled inputs,
# s_d = ((x_d - y_d) / l_d)^2, one matrix per input dimension. Their sum is
# exactly the scaled squared distance computed by
# .scaled_pairwise_squared_distance().
.scaled_squared_difference_components <- function(x, y, length_scale) {
  resolved <- .resolve_length_scale(length_scale, ncol(x))

  lapply(
    seq_len(ncol(x)),
    function(dimension) {
      outer(
        x[, dimension] / resolved[[dimension]],
        y[, dimension] / resolved[[dimension]],
        FUN = "-"
      )^2
    }
  )
}


.pairwise_periodic_scaled_distance <- function(
    x,
    y,
    length_scale,
    period) {
  resolved <- .resolve_length_scale(length_scale, ncol(x))
  result <- matrix(0, nrow = nrow(x), ncol = nrow(y))

  for (dimension in seq_len(ncol(x))) {
    differences <- outer(
      x[, dimension],
      y[, dimension],
      FUN = "-"
    )
    periodic_component <- sin(pi * differences / period)^2
    result <- result + periodic_component / (resolved[[dimension]]^2)
  }

  result
}
