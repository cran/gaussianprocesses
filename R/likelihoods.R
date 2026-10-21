# Likelihoods ------------------------------------------------------------------
#
# A likelihood specification is data only: list(family, link, parameters)
# with class "gaussianprocesses_likelihood", like a kernel specification, so
# it serializes safely. Its functions live in the internal table
# .likelihood_entry(), one entry per family and link, with
#   log_density(y, f, parameters, exposure)   log p(y_i | f_i),
#   derivatives(y, f, parameters, exposure)   d1, d2, d3: derivatives of the
#                                             log density in f_i,
#   parameter_derivatives(...)                derivatives of the log density,
#                                             d1, and d2 in each parameter's
#                                             optimizer coordinate,
#   predictive(mean, variance, parameters, exposure, order)
#                                             mean and variance of y_* when
#                                             f_* ~ N(mean, variance),
#   predictive_log_density(y, mean, variance, parameters, exposure, order)
#                                             log int p(y | f) N(f) df in
#                                             closed form, or NULL for
#                                             adaptive Gauss-Hermite
#                                             quadrature,
#   simulate(f, parameters, exposure)         responses drawn from
#                                             p(y | f), with the shape of f,
#   conditional_mean(f, parameters, exposure) E[y | f], increasing in f: the
#                                             class probability or the rate,
# a log-concavity flag, and response validation. Exposure multiplies the
# Poisson rate and is data, not a parameter.


# Validates a likelihood specification, including one restored with
# readRDS(), and returns its entry.
.likelihood_entry_of <- function(likelihood, name = "likelihood") {
  is_scalar_string <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value)
  }

  if (!inherits(likelihood, "gaussianprocesses_likelihood") ||
      !is.list(likelihood) ||
      !is_scalar_string(likelihood$family) ||
      !is_scalar_string(likelihood$link) ||
      !is.list(likelihood$parameters)) {
    .gp_abort(
      sprintf("'%s' must be a likelihood specification.", name),
      "gaussianprocesses_validation_error"
    )
  }

  entry <- .likelihood_entry(likelihood)
  definitions <- entry$parameters

  if (!setequal(names(likelihood$parameters), names(definitions))) {
    .gp_abort(
      sprintf(
        "'%s' must have the parameters %s.",
        name,
        if (length(definitions) == 0L) {
          "none"
        } else {
          paste(names(definitions), collapse = ", ")
        }
      ),
      "gaussianprocesses_validation_error"
    )
  }

  for (parameter in names(definitions)) {
    .validate_kernel_parameter(
      likelihood$parameters[[parameter]],
      parameter,
      definitions[[parameter]]
    )
  }

  entry
}


# Responses ---------------------------------------------------------------------

# Binary responses as 0/1: numeric 0/1, logical (TRUE is 1), or a factor with
# two levels whose second level is 1, as in glm().
.binary_response <- function(y, name = "y") {
  if (is.factor(y)) {
    if (nlevels(y) != 2L || anyNA(y)) {
      .gp_abort(
        sprintf(
          "Factor '%s' must have exactly two levels and no missing values.",
          name
        ),
        "gaussianprocesses_validation_error"
      )
    }

    return(as.numeric(as.integer(y) == 2L))
  }

  if (is.logical(y)) {
    if (anyNA(y) || length(y) == 0L) {
      .gp_abort(
        sprintf("Logical '%s' must be non-empty without missing values.", name),
        "gaussianprocesses_validation_error"
      )
    }

    return(as.numeric(y))
  }

  if (!is.numeric(y) || !is.null(dim(y)) || length(y) == 0L ||
      any(!(y %in% c(0, 1)))) {
    .gp_abort(
      sprintf(
        "'%s' must be 0/1 numbers, logical values, or a two-level factor.",
        name
      ),
      "gaussianprocesses_validation_error"
    )
  }

  as.numeric(y)
}


.count_response <- function(y, name = "y") {
  if (!is.numeric(y) || !is.null(dim(y)) || length(y) == 0L ||
      any(!is.finite(y)) || any(y < 0) || any(y != round(y))) {
    .gp_abort(
      sprintf("'%s' must be non-negative integer counts.", name),
      "gaussianprocesses_validation_error"
    )
  }

  as.numeric(y)
}


.continuous_response <- function(y, name = "y") {
  .assert_numeric_vector(y, name)
  as.numeric(y)
}


# Exposure per observation: NULL (1), one positive value, or one per
# observation.
.resolve_exposure <- function(exposure, n_observations) {
  if (is.null(exposure)) {
    return(rep(1, n_observations))
  }

  if (!is.numeric(exposure) || !is.null(dim(exposure)) ||
      !(length(exposure) %in% c(1L, n_observations)) ||
      any(!is.finite(exposure)) || any(exposure <= 0)) {
    .gp_abort(
      sprintf(
        "'exposure' must be one positive value or one per observation (%d).",
        n_observations
      ),
      "gaussianprocesses_validation_error"
    )
  }

  rep_len(as.numeric(exposure), n_observations)
}


# Families ------------------------------------------------------------------------

.gaussian_likelihood_entry <- function() {
  list(
    label = "GaussianLikelihood",
    log_concave = TRUE,
    response = .continuous_response,
    parameters = list(variance = .kernel_parameter(1)),
    log_density = function(y, f, parameters, exposure) {
      variance <- parameters$variance
      -0.5 * log(2 * pi * variance) - (y - f)^2 / (2 * variance)
    },
    derivatives = function(y, f, parameters, exposure) {
      list(
        first = (y - f) / parameters$variance,
        second = rep(-1 / parameters$variance, length(f)),
        third = numeric(length(f))
      )
    },
    # In the log coordinate of the variance s: d log p = -1/2 + r^2 / (2 s),
    # d d1 = -r / s, d d2 = 1 / s, with r = y - f.
    parameter_derivatives = function(y, f, parameters, exposure) {
      residual <- y - f
      variance <- parameters$variance

      list(
        variance = list(
          log_density = -0.5 + residual^2 / (2 * variance),
          first = -residual / variance,
          second = rep(1 / variance, length(f))
        )
      )
    },
    predictive = function(mean, variance, parameters, exposure, order) {
      list(mean = mean, variance = variance + parameters$variance)
    },
    predictive_log_density = function(y, mean, variance, parameters,
                                      exposure, order, ...) {
      stats::dnorm(y, mean, sqrt(variance + parameters$variance), log = TRUE)
    },
    simulate = function(f, parameters, exposure) {
      f + sqrt(parameters$variance) * stats::rnorm(length(f))
    },
    conditional_mean = function(f, parameters, exposure) f
  )
}


# Bernoulli with the logistic link, pi = 1 / (1 + exp(-f)):
#   log p = y f - log(1 + e^f), d1 = y - pi, d2 = -pi (1 - pi),
#   d3 = -pi (1 - pi) (1 - 2 pi).
# pi and 1 - pi are computed separately, without cancellation.
.logit_likelihood_entry <- function() {
  list(
    label = "BernoulliLikelihood",
    log_concave = TRUE,
    response = .binary_response,
    parameters = list(),
    log_density = function(y, f, parameters, exposure) {
      stats::plogis((2 * y - 1) * f, log.p = TRUE)
    },
    derivatives = function(y, f, parameters, exposure) {
      success <- stats::plogis(f)
      failure <- stats::plogis(-f)
      spread <- success * failure

      list(
        first = y * failure - (1 - y) * success,
        second = -spread,
        third = -spread * (failure - success)
      )
    },
    parameter_derivatives = NULL,
    predictive = function(mean, variance, parameters, exposure, order) {
      probability <- exp(.predictive_log_density(
        rep(1, length(mean)),
        mean,
        variance,
        parameters,
        exposure,
        order,
        entry = .logit_likelihood_entry()
      ))
      list(
        mean = probability,
        variance = probability * (1 - probability),
        probability = probability
      )
    },
    predictive_log_density = NULL,
    simulate = function(f, parameters, exposure) {
      .shaped_like(stats::rbinom(length(f), 1L, stats::plogis(f)), f)
    },
    conditional_mean = function(f, parameters, exposure) stats::plogis(f)
  )
}


# The inverse Mills ratio r = phi(z) / Phi(z), the shift z + r, and the
# curvature term (z + r)(z + 2 r) - 1. For z >= -5, r = exp(log phi -
# log Phi) directly. Below, r grows like -z and both z + r and the curvature
# term cancel, so with u = -z they come from the continued fraction of the
# Mills ratio, r = u + 1 / K_1 with K_m = u + (m + 1) / K_(m + 1):
#   z + r = 1 / K_1,
#   (z + r)(z + 2 r) - 1 = 2 (3 / K_3 - 2 / K_2) / (K_1^2 K_2),
# which follow from u K_1 - K_1^2 = -2 K_1 / K_2 and have no cancellation.
.probit_ratio_terms <- function(z) {
  ratio <- exp(stats::dnorm(z, log = TRUE) - stats::pnorm(z, log.p = TRUE))
  shift <- z + ratio
  curvature <- shift * (z + 2 * ratio) - 1
  tail <- z < -5

  if (any(tail)) {
    u <- -z[tail]
    k3 <- u

    for (m in 60:3) {
      k3 <- u + (m + 1) / k3
    }

    k2 <- u + 3 / k3
    k1 <- u + 2 / k2
    ratio[tail] <- u + 1 / k1
    shift[tail] <- 1 / k1
    curvature[tail] <- 2 * (3 / k3 - 2 / k2) / (k1^2 * k2)
  }

  list(ratio = ratio, shift = shift, curvature = curvature)
}


# Bernoulli with the probit link, p(y | f) = Phi(t f) with t = 2 y - 1. With
# z = t f and the inverse Mills ratio r(z) = phi(z) / Phi(z),
#   d1 = t r,  d2 = -r (z + r),  d3 = t r ((z + r)(z + 2 r) - 1).
# The predictive probability is Phi(mu / sqrt(1 + sigma^2)).
.probit_likelihood_entry <- function() {
  list(
    label = "BernoulliLikelihood",
    log_concave = TRUE,
    response = .binary_response,
    parameters = list(),
    log_density = function(y, f, parameters, exposure) {
      stats::pnorm((2 * y - 1) * f, log.p = TRUE)
    },
    derivatives = function(y, f, parameters, exposure) {
      sign <- 2 * y - 1
      terms <- .probit_ratio_terms(sign * f)

      list(
        first = sign * terms$ratio,
        second = -terms$ratio * terms$shift,
        third = sign * terms$ratio * terms$curvature
      )
    },
    parameter_derivatives = NULL,
    predictive = function(mean, variance, parameters, exposure, order) {
      probability <- stats::pnorm(mean / sqrt(1 + variance))
      list(
        mean = probability,
        variance = probability * (1 - probability),
        probability = probability
      )
    },
    predictive_log_density = function(y, mean, variance, parameters,
                                      exposure, order, ...) {
      stats::pnorm((2 * y - 1) * mean / sqrt(1 + variance), log.p = TRUE)
    },
    simulate = function(f, parameters, exposure) {
      .shaped_like(stats::rbinom(length(f), 1L, stats::pnorm(f)), f)
    },
    conditional_mean = function(f, parameters, exposure) stats::pnorm(f)
  )
}


# Poisson with the log link and exposure E, rate lambda = E e^f:
#   log p = y (f + log E) - E e^f - log y!,  d1 = y - E e^f,
#   d2 = d3 = -E e^f.
# For f ~ N(mu, sigma^2), lambda is lognormal, so
#   E[y] = E e^(mu + sigma^2 / 2),
#   Var[y] = E[y] + E^2 e^(2 mu + sigma^2) (e^(sigma^2) - 1).
.poisson_likelihood_entry <- function() {
  list(
    label = "PoissonLikelihood",
    log_concave = TRUE,
    response = .count_response,
    parameters = list(),
    log_density = function(y, f, parameters, exposure) {
      y * (f + log(exposure)) - exposure * exp(f) - lgamma(y + 1)
    },
    derivatives = function(y, f, parameters, exposure) {
      rate <- exposure * exp(f)
      list(first = y - rate, second = -rate, third = -rate)
    },
    parameter_derivatives = NULL,
    predictive = function(mean, variance, parameters, exposure, order) {
      expected <- exposure * exp(mean + variance / 2)
      list(
        mean = expected,
        variance = expected + expected^2 * expm1(variance)
      )
    },
    predictive_log_density = NULL,
    simulate = function(f, parameters, exposure) {
      .shaped_like(stats::rpois(length(f), exposure * exp(f)), f)
    },
    conditional_mean = function(f, parameters, exposure) exposure * exp(f)
  )
}


# Draws as doubles with the dimensions of the latent values.
.shaped_like <- function(values, f) {
  values <- as.numeric(values)
  dim(values) <- dim(f)
  values
}


.likelihood_entry <- function(likelihood) {
  switch(
    paste(likelihood$family, likelihood$link),
    "gaussian identity" = .gaussian_likelihood_entry(),
    "bernoulli logit" = .logit_likelihood_entry(),
    "bernoulli probit" = .probit_likelihood_entry(),
    "poisson log" = .poisson_likelihood_entry(),
    .gp_abort(
      sprintf(
        "Unknown likelihood '%s' with link '%s'.",
        likelihood$family,
        likelihood$link
      ),
      "gaussianprocesses_validation_error"
    )
  )
}


.new_likelihood <- function(family, link, parameters = list()) {
  likelihood <- structure(
    list(family = family, link = link, parameters = parameters),
    class = "gaussianprocesses_likelihood"
  )
  .likelihood_entry_of(likelihood)

  likelihood
}


# Constructors ------------------------------------------------------------------

#' Likelihoods for Gaussian-process models
#'
#' Observation models \eqn{p(y_i \mid f_i)} that link the latent Gaussian
#' process \eqn{f} to the responses. Like kernel specifications, they are
#' plain data with a class, a print method, and constrained parameters.
#'
#' @param variance Positive noise variance of the Gaussian likelihood.
#' @param link Link function: `"logit"` or `"probit"` for Bernoulli
#'   responses, and `"log"` for Poisson counts.
#'
#' @details
#' \describe{
#'   \item{`gaussian_likelihood()`}{\eqn{y_i \sim N(f_i, \sigma^2)}, the
#'   observation model of [fit_gp()] with `noise_variance` \eqn{\sigma^2}.}
#'   \item{`bernoulli_likelihood()`}{\eqn{\Pr(y_i = 1) = \sigma(f_i)} with
#'   the logistic function, \eqn{\log p = y f - \log(1 + e^f)}, or
#'   \eqn{\Phi(f_i)} with the probit link. Responses may be 0/1 numbers,
#'   logical values (`TRUE` is 1), or a two-level factor whose second level
#'   is 1, as in [stats::glm()]. Multi-class (softmax) likelihoods are not
#'   supported.}
#'   \item{`poisson_likelihood()`}{\eqn{y_i \sim \mathrm{Poisson}(E_i
#'   e^{f_i})} with exposure \eqn{E_i}, so
#'   \eqn{\log p = y (f + \log E) - E e^f - \log y!}. The exposure is data,
#'   supplied where the likelihood is used, not part of the specification.
#'   Counts more variable than a Poisson process with a smooth rate allows
#'   need negative-binomial or zero-inflated likelihoods, which are not yet
#'   supported; [gp_count_scores()] checks for such overdispersion.}
#' }
#'
#' Every likelihood provides, vectorized over observations, the log density
#' and its first three derivatives in \eqn{f_i} ([evaluate_likelihood()]),
#' derivatives in its parameters, the predictive distribution of a new
#' response when \eqn{f_*} is Gaussian ([likelihood_predictive()]), and
#' whether the log density is concave in \eqn{f}. All three are
#' log-concave, which makes Laplace approximations well behaved. Gaussian
#' expectations without a closed form use Gauss-Hermite quadrature: nodes
#' from the Golub-Welsch eigenvalue method, weights from Christoffel sums,
#' cached by order.
#'
#' The probit derivatives use the inverse Mills ratio
#' \eqn{r = \varphi(z) / \Phi(z)} computed in log space, with
#' \eqn{z = (2y - 1) f}. For \eqn{z < -5}, where \eqn{r} grows like
#' \eqn{-z} and \eqn{z + r} cancels, they use the continued fraction of the
#' Mills ratio instead. The logistic derivatives compute \eqn{\pi} and
#' \eqn{1 - \pi} separately. All derivatives stay finite and accurate for
#' large \eqn{|f|}.
#'
#' Latent models with any of these likelihoods are fitted with
#' [fit_latent_gp()] and [optimize_latent_gp()]; [fit_gp()] remains the
#' exact Gaussian model.
#'
#' @return An object of class `gaussianprocesses_likelihood`.
#' @references Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss
#'   quadrature rules. *Mathematics of Computation*, 23(106), 221--230.
#'
#'   Rasmussen, C. E. and Williams, C. K. I. (2006). *Gaussian Processes for
#'   Machine Learning*. MIT Press. Chapter 3.
#' @seealso [evaluate_likelihood()], [likelihood_predictive()],
#'   [fit_latent_gp()].
#' @examples
#' gaussian_likelihood(variance = 0.1)
#' bernoulli_likelihood("probit")
#' poisson_likelihood()
#'
#' evaluate_likelihood(bernoulli_likelihood(), c(TRUE, FALSE), f = c(1, 1))
#' @template experimental
#' @name gp_likelihoods
NULL


#' @rdname gp_likelihoods
#' @export
gaussian_likelihood <- function(variance = 1) {
  .new_likelihood("gaussian", "identity", list(variance = variance))
}


#' @rdname gp_likelihoods
#' @export
bernoulli_likelihood <- function(link = c("logit", "probit")) {
  .new_likelihood("bernoulli", match.arg(link))
}


#' @rdname gp_likelihoods
#' @export
poisson_likelihood <- function(link = "log") {
  .new_likelihood("poisson", match.arg(link))
}


#' @export
print.gaussianprocesses_likelihood <- function(x, ...) {
  entry <- .likelihood_entry_of(x, "x")
  details <- c(
    if (!identical(x$family, "gaussian")) paste0("link=", x$link),
    if (length(x$parameters) > 0L) .format_kernel_parameters(x$parameters)
  )
  cat(entry$label, "(", paste(details, collapse = ", "), ")\n", sep = "")
  invisible(x)
}


# Evaluation ------------------------------------------------------------------------

#' Evaluate a likelihood and its derivatives
#'
#' Computes \eqn{\log p(y_i \mid f_i)} and its first three derivatives with
#' respect to the latent values \eqn{f_i}, the quantities that Laplace and
#' related approximations need.
#'
#' @param likelihood A likelihood specification ([gp_likelihoods]).
#' @param y Responses, as the likelihood accepts them.
#' @param f Numeric vector of latent values, one per response.
#' @param exposure Optional positive exposure for Poisson likelihoods: one
#'   value or one per response.
#'
#' @return A data frame with one row per response: `y` in the internal coding
#'   (0/1 for binary responses), `f`, `log_density`, and the derivatives
#'   `first`, `second`, and `third`.
#' @examples
#' f <- c(-40, -2, 0, 2, 40)
#' evaluate_likelihood(bernoulli_likelihood("probit"), rep(1, 5), f)
#' evaluate_likelihood(poisson_likelihood(), c(0, 1, 2, 3, 4), f = log(1:5),
#'   exposure = 2)
#' @template experimental
#' @export
evaluate_likelihood <- function(likelihood, y, f, exposure = NULL) {
  entry <- .likelihood_entry_of(likelihood)
  y <- entry$response(y)
  .assert_numeric_vector(f, "f")

  if (length(f) != length(y)) {
    .gp_abort(
      "'y' and 'f' must have the same length.",
      "gaussianprocesses_validation_error"
    )
  }

  exposure <- .likelihood_exposure(likelihood, exposure, length(y))
  derivatives <- entry$derivatives(y, f, likelihood$parameters, exposure)

  data.frame(
    y = y,
    f = f,
    log_density = entry$log_density(y, f, likelihood$parameters, exposure),
    first = derivatives$first,
    second = derivatives$second,
    third = derivatives$third
  )
}


# Exposure applies to Poisson likelihoods only.
.likelihood_exposure <- function(likelihood, exposure, n_observations) {
  if (!is.null(exposure) && !identical(likelihood$family, "poisson")) {
    .gp_abort(
      "'exposure' applies to Poisson likelihoods only.",
      "gaussianprocesses_validation_error"
    )
  }

  .resolve_exposure(exposure, n_observations)
}


#' Predictive distribution of new responses
#'
#' Integrates a likelihood against a Gaussian distribution of the latent
#' value, \eqn{\int p(y_* \mid f_*) N(f_* \mid \mu, \sigma^2) df_*}: the
#' predictive distribution of a new response given the posterior of the
#' latent function at its input.
#'
#' @inheritParams evaluate_likelihood
#' @param mean,variance Mean and variance of the latent value, one per
#'   prediction.
#' @param observed Optional responses whose log predictive density is
#'   returned.
#' @param quadrature_order Number of Gauss-Hermite nodes for integrals
#'   without a closed form.
#'
#' @details
#' Closed forms are used where they exist: the Gaussian likelihood adds its
#' variance, the probit link gives \eqn{\Pr(y_* = 1) =
#' \Phi(\mu / \sqrt{1 + \sigma^2})}, and the Poisson mean and variance are
#' those of a Poisson with a lognormal rate,
#' \eqn{E e^{\mu + \sigma^2 / 2}} and
#' \eqn{\mathrm{E}[y_*] + E^2 e^{2\mu + \sigma^2}(e^{\sigma^2} - 1)}.
#' The logistic probability and the Poisson and logistic predictive
#' densities use adaptive Gauss-Hermite quadrature (Liu and Pierce, 1994).
#' With \eqn{\hat f} the mode of \eqn{\psi(f) = \log p(y \mid f) +
#' \log N(f \mid \mu, \sigma^2)}, found by Newton's method, and
#' \eqn{s = (-\psi''(\hat f))^{-1/2}},
#' \deqn{\int p(y \mid f) N(f \mid \mu, \sigma^2) df \approx \sqrt{2} s
#' \sum_j w_j e^{t_j^2} e^{\psi(\hat f + \sqrt{2} s t_j)},}
#' computed in logs. Centring the nodes on the integrand rather than on
#' \eqn{\mu} matters when the observation lies in the tail of the latent
#' distribution, as for large counts.
#'
#' Quadrature converges more slowly as \eqn{\sigma^2} grows, because the
#' logistic function has poles at \eqn{\pm i\pi}. With the default 40
#' nodes, the logistic probability is accurate to \eqn{10^{-10}} for
#' \eqn{\sigma^2} up to 3; increase `quadrature_order` for wider latent
#' distributions.
#'
#' @references Liu, Q. and Pierce, D. A. (1994). A note on Gauss-Hermite
#'   quadrature. *Biometrika*, 81(3), 624--629.
#' @return A list with the predictive `mean` and `variance` of \eqn{y_*},
#'   `probability`, \eqn{\Pr(y_* = 1)}, for Bernoulli likelihoods, and
#'   `log_density` of `observed` when it is given.
#' @examples
#' likelihood_predictive(bernoulli_likelihood("probit"), mean = c(-1, 0, 2),
#'   variance = c(0.5, 1, 4))$probability
#'
#' likelihood_predictive(poisson_likelihood(), mean = log(3), variance = 0.2,
#'   observed = 4, exposure = 2)
#' @template experimental
#' @export
likelihood_predictive <- function(
    likelihood,
    mean,
    variance,
    observed = NULL,
    exposure = NULL,
    quadrature_order = 40L) {
  entry <- .likelihood_entry_of(likelihood)
  .assert_numeric_vector(mean, "mean")
  .assert_numeric_vector(variance, "variance")
  quadrature_order <- .assert_positive_integer(
    quadrature_order,
    "quadrature_order"
  )

  if (length(variance) != length(mean) || any(variance < 0)) {
    .gp_abort(
      "'variance' must be non-negative, one value per mean.",
      "gaussianprocesses_validation_error"
    )
  }

  exposure <- .likelihood_exposure(likelihood, exposure, length(mean))
  result <- entry$predictive(
    mean,
    variance,
    likelihood$parameters,
    exposure,
    quadrature_order
  )

  if (!is.null(observed)) {
    observed <- entry$response(observed, "observed")

    if (length(observed) != length(mean)) {
      .gp_abort(
        "'observed' must have one value per mean.",
        "gaussianprocesses_validation_error"
      )
    }

    log_density <- if (is.null(entry$predictive_log_density)) {
      .predictive_log_density
    } else {
      entry$predictive_log_density
    }
    result$log_density <- log_density(
      observed,
      mean,
      variance,
      likelihood$parameters,
      exposure,
      quadrature_order,
      entry = entry
    )
  }

  result
}


# log int p(y | f) N(f | mean, variance) df by adaptive Gauss-Hermite
# quadrature, centred at the mode of the integrand and scaled by its
# curvature there. Newton's method with step halving finds the mode; for a
# log-concave likelihood the log integrand is strictly concave. Negative
# curvature of the likelihood is clipped at zero, so that the Newton
# direction stays an ascent direction for any likelihood.
.predictive_log_density <- function(
    y,
    mean,
    variance,
    parameters,
    exposure,
    order,
    entry) {
  result <- entry$log_density(y, mean, parameters, exposure)
  random <- variance > 0

  if (!any(random)) {
    return(result)
  }

  y <- y[random]
  mean <- mean[random]
  variance <- variance[random]
  exposure <- exposure[random]
  log_integrand <- function(f) {
    entry$log_density(y, f, parameters, exposure) -
      (f - mean)^2 / (2 * variance)
  }
  curvature_at <- function(f) {
    pmin(entry$derivatives(y, f, parameters, exposure)$second, 0) -
      1 / variance
  }

  mode <- mean
  value <- log_integrand(mode)

  for (iteration in seq_len(100L)) {
    derivatives <- entry$derivatives(y, mode, parameters, exposure)
    curvature <- pmin(derivatives$second, 0) - 1 / variance
    step <- -(derivatives$first - (mode - mean) / variance) / curvature

    if (all(abs(step) * sqrt(-curvature) <= 1e-10)) {
      break
    }

    fraction <- rep(1, length(mode))
    candidate <- mode + step
    candidate_value <- log_integrand(candidate)

    for (halving in seq_len(60L)) {
      worse <- !(candidate_value >= value - 1e-12 * (1 + abs(value)))

      if (!any(worse)) {
        break
      }

      fraction[worse] <- fraction[worse] / 2
      candidate[worse] <- mode[worse] + fraction[worse] * step[worse]
      candidate_value <- log_integrand(candidate)
    }

    mode <- candidate
    value <- candidate_value
  }

  result[random] <- .gaussian_log_expectation(
    function(f, index) {
      entry$log_density(
        y[index],
        f,
        parameters,
        exposure[index]
      )
    },
    mean,
    variance,
    order,
    center = mode,
    scale = 1 / sqrt(-curvature_at(mode))
  )
  result
}
