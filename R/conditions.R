# Internal condition helpers -------------------------------------------------

# `...` adds named data fields to the condition object.
.gp_abort <- function(message, class, ...) {
  condition <- structure(
    c(list(message = message, call = NULL), list(...)),
    class = c(class, "gaussianprocesses_error", "error", "condition")
  )

  stop(condition)
}


