# Gauss-Hermite quadrature -------------------------------------------------------
#
# Nodes t_j and weights w_j with sum_j w_j g(t_j) ~ int g(t) exp(-t^2) dt,
# exact for polynomials of degree up to 2n - 1. Golub and Welsch (1969): the
# nodes are the eigenvalues of the symmetric tridiagonal Jacobi matrix of the
# Hermite recurrence, with zero diagonal and off-diagonal sqrt(k / 2). The
# weights are Christoffel numbers (below). For a Gaussian expectation,
#   int g(f) N(f | mu, sigma^2) df ~ pi^(-1/2) sum_j w_j g(mu + sqrt(2) sigma t_j).
# Rules are cached by order.

.quadrature_cache <- new.env(parent = emptyenv())


.gauss_hermite <- function(order) {
  order <- .assert_positive_integer(order, "order")
  key <- as.character(order)

  if (!is.null(.quadrature_cache[[key]])) {
    return(.quadrature_cache[[key]])
  }

  nodes <- 0

  if (order > 1L) {
    jacobi <- matrix(0, order, order)
    off_diagonal <- sqrt(seq_len(order - 1L) / 2)
    jacobi[cbind(seq_len(order - 1L), 2:order)] <- off_diagonal
    jacobi[cbind(2:order, seq_len(order - 1L))] <- off_diagonal
    nodes <- sort(eigen(jacobi, symmetric = TRUE, only.values = TRUE)$values)

    # The rule is symmetric; symmetrize away rounding.
    nodes <- (nodes - rev(nodes)) / 2
  }

  # Christoffel numbers w_j = 1 / sum_k p_k(t_j)^2 over the orthonormal
  # Hermite polynomials p_0 = pi^(-1/4), p_1 = sqrt(2) t p_0,
  #   p_(k+1) = sqrt(2 / (k + 1)) t p_k - sqrt(k / (k + 1)) p_(k-1),
  # for k < order. A sum of positive terms keeps full relative accuracy for
  # the smallest weights, which sqrt(pi) v_1j^2 from the eigenvectors loses.
  previous <- rep(0, order)
  current <- rep(pi^(-1 / 4), order)
  total <- current^2

  for (k in seq_len(order - 1L) - 1L) {
    following <- sqrt(2 / (k + 1)) * nodes * current -
      sqrt(k / (k + 1)) * previous
    previous <- current
    current <- following
    total <- total + current^2
  }

  rule <- list(nodes = nodes, weights = 1 / total)
  assign(key, rule, envir = .quadrature_cache)
  rule
}


# Expectations E[g(f)] under N(mean, variance), one per element of `mean`.
# `g(f, index)` evaluates the integrand at a matrix of latent values f with
# one row per element and one column per node; `index` gives the element of
# each row.
.gaussian_expectation <- function(g, mean, variance, order) {
  rule <- .gauss_hermite(order)
  latent <- outer(mean, rep(1, order)) +
    outer(sqrt(2 * variance), rule$nodes)
  values <- g(latent, seq_along(mean))

  drop(values %*% rule$weights) / sqrt(pi)
}


# log E[exp(h(f))] under N(mean, variance), computed with log-sum-exp so
# that small densities do not underflow. With `center` and `scale`, the rule
# is placed at f_j = center + sqrt(2) scale t_j instead (adaptive
# Gauss-Hermite; Liu and Pierce, 1994):
#   int e^h(f) N(f) df = sqrt(2) scale sum_j w_j e^(t_j^2) e^h(f_j) N(f_j),
# which is exact for every order when h is quadratic and center and scale are
# the mode and curvature scale of h + log N. When e^h is peaked away from the
# mean this needs far fewer nodes.
.gaussian_log_expectation <- function(
    h,
    mean,
    variance,
    order,
    center = NULL,
    scale = NULL) {
  rule <- .gauss_hermite(order)
  ones <- rep(1, length(mean))

  if (is.null(center)) {
    latent <- outer(mean, rep(1, order)) +
      outer(sqrt(2 * variance), rule$nodes)
    terms <- h(latent, seq_along(mean)) + outer(ones, log(rule$weights))
    offset <- -0.5 * log(pi)
  } else {
    latent <- outer(center, rep(1, order)) + outer(sqrt(2) * scale, rule$nodes)
    terms <- h(latent, seq_along(mean)) +
      stats::dnorm(latent, mean, sqrt(variance), log = TRUE) +
      outer(ones, log(rule$weights) + rule$nodes^2)
    offset <- log(sqrt(2) * scale)
  }

  largest <- apply(terms, 1L, max)

  largest + log(rowSums(exp(terms - largest))) + offset
}
