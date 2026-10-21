# Internal validation helpers ------------------------------------------------

.default_tolerance <- function() {
  sqrt(.Machine$double.eps)
}


.assert_numeric_scalar <- function(
    x,
    name,
    lower = -Inf,
    lower_open = FALSE,
    finite = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    .gp_abort(
      sprintf("'%s' must be a single numeric value.", name),
      "gaussianprocesses_validation_error"
    )
  }

  if (finite && !is.finite(x)) {
    .gp_abort(
      sprintf("'%s' must be finite.", name),
      "gaussianprocesses_validation_error"
    )
  }

  violates_lower_bound <- if (lower_open) x <= lower else x < lower

  if (violates_lower_bound) {
    relation <- if (lower_open) "greater than" else "greater than or equal to"
    .gp_abort(
      sprintf("'%s' must be %s %s.", name, relation, format(lower)),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(x)
}


.assert_positive_integer <- function(x, name) {
  .assert_numeric_scalar(x, name, lower = 1)

  if (x != floor(x)) {
    .gp_abort(
      sprintf("'%s' must be an integer greater than or equal to 1.", name),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(as.integer(x))
}


.assert_numeric_vector <- function(
    x,
    name = deparse(substitute(x)),
    finite = TRUE,
    allow_empty = FALSE) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    .gp_abort(
      sprintf("'%s' must be a numeric vector.", name),
      "gaussianprocesses_validation_error"
    )
  }

  if (!allow_empty && length(x) == 0L) {
    .gp_abort(
      sprintf("'%s' must not be empty.", name),
      "gaussianprocesses_validation_error"
    )
  }

  if (finite && any(!is.finite(x))) {
    .gp_abort(
      sprintf("'%s' must contain only finite values.", name),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(x)
}


.assert_numeric_matrix <- function(
    x,
    name = deparse(substitute(x)),
    finite = TRUE,
    allow_empty = FALSE) {
  if (!is.matrix(x) || !is.numeric(x)) {
    .gp_abort(
      sprintf("'%s' must be a numeric matrix.", name),
      "gaussianprocesses_validation_error"
    )
  }

  if (!allow_empty && (nrow(x) == 0L || ncol(x) == 0L)) {
    .gp_abort(
      sprintf("'%s' must have at least one row and one column.", name),
      "gaussianprocesses_validation_error"
    )
  }

  if (finite && any(!is.finite(x))) {
    .gp_abort(
      sprintf("'%s' must contain only finite values.", name),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(x)
}


.assert_square_matrix <- function(
    x,
    name = deparse(substitute(x)),
    finite = TRUE) {
  .assert_numeric_matrix(x, name = name, finite = finite)

  if (nrow(x) != ncol(x)) {
    .gp_abort(
      sprintf("'%s' must be square.", name),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(x)
}


.assert_symmetric_matrix <- function(
    x,
    name = deparse(substitute(x)),
    tolerance = .default_tolerance()) {
  .assert_square_matrix(x, name = name)
  .assert_numeric_scalar(tolerance, "tolerance", lower = 0)

  magnitude <- max(1, max(abs(x)))
  asymmetry <- max(abs(x - t(x)))

  if (asymmetry > tolerance * magnitude) {
    .gp_abort(
      sprintf(
        "'%s' must be symmetric within relative tolerance %.3e.",
        name,
        tolerance
      ),
      "gaussianprocesses_validation_error"
    )
  }

  invisible(x)
}
