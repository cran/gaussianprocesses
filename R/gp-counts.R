# Predictive distributions and scores of counts -------------------------------
#
# With f ~ N(mu, sigma^2) and y | f ~ Poisson(E e^f), the predictive
# distribution of a count is a Poisson-lognormal mixture. Its quantiles are
# found by bisection over the integers, and its probabilities
# p(k) = F(k) - F(k - 1) come from the adaptive quadrature of
# likelihood_predictive().


# The distribution function F(k) = int ppois(k, E e^f) N(f | mu, sigma^2) df.
# As a function of f, ppois(k, E e^f) falls from 1 to 0 over a width of about
# 1 / sqrt(k + 1), too sharply for Gauss-Hermite quadrature in f when sigma
# is larger. Because ppois(k, lambda) = Pr(G > lambda) with
# G ~ Gamma(k + 1, 1) independent of f,
#   F(k) = Pr(E e^f < G) = E_G[Phi((log(G / E) - mu) / sigma)],
# whose integrand is smooth in z when G = qgamma(Phi(z), k + 1) and sigma is
# large. Each count uses the form in which its integrand is smooth: f when
# sigma sqrt(k + 1) <= 1, and z otherwise.
.count_distribution <- function(k, mean, variance, exposure, order) {
  sd <- sqrt(variance)
  shape <- k + 1
  result <- numeric(length(k))
  latent_form <- which(k >= 0 & sd * sqrt(shape) <= 1)
  gamma_form <- which(k >= 0 & sd * sqrt(shape) > 1)

  if (length(latent_form) > 0L) {
    i <- latent_form
    result[i] <- .gaussian_expectation(
      function(f, index) {
        matrix(
          stats::pgamma(
            exposure[i][index] * exp(f),
            shape[i][index],
            lower.tail = FALSE
          ),
          nrow = nrow(f)
        )
      },
      mean[i],
      variance[i],
      order
    )
  }

  if (length(gamma_form) > 0L) {
    i <- gamma_form
    result[i] <- .gaussian_expectation(
      function(z, index) {
        # Quantiles of G from log probabilities of the nearer tail, accurate
        # in both tails.
        shapes <- matrix(shape[i][index], nrow(z), ncol(z))
        lower <- z <= 0
        gamma <- z
        gamma[lower] <- stats::qgamma(
          stats::pnorm(z[lower], log.p = TRUE),
          shapes[lower],
          log.p = TRUE
        )
        gamma[!lower] <- stats::qgamma(
          stats::pnorm(-z[!lower], log.p = TRUE),
          shapes[!lower],
          lower.tail = FALSE,
          log.p = TRUE
        )

        matrix(
          stats::pnorm(
            (log(gamma / exposure[i][index]) - mean[i][index]) /
              sd[i][index]
          ),
          nrow = nrow(z)
        )
      },
      rep(0, length(i)),
      rep(1, length(i)),
      order
    )
  }

  result
}


# Smallest k >= 0 with F(k) >= probability, for each prediction. The search
# keeps a lower bound (F(lower - 1) < probability) and an upper bound
# (F(upper) >= probability). It starts from the quantile of the negative
# binomial distribution with the same mean and variance, which is usually
# within a count of the answer, widens the bracket in doubling steps, and
# finishes by bisection.
.count_quantile <- function(probability, mean, variance, exposure, order) {
  n <- length(mean)
  probability <- rep_len(probability, n)
  distribution <- function(k, rows) {
    .count_distribution(k, mean[rows], variance[rows], exposure[rows], order)
  }

  # Poisson-lognormal mean m = E e^(mu + sigma^2 / 2) and variance
  # m + m^2 (e^(sigma^2) - 1): a negative binomial of size 1 / expm1(sigma^2).
  rate_mean <- exposure * exp(mean + variance / 2)
  size <- 1 / expm1(variance)
  guess <- ifelse(
    is.finite(size) & size < 1e12,
    stats::qnbinom(probability, size = pmin(size, 1e12), mu = rate_mean),
    stats::qpois(probability, rate_mean)
  )

  lower <- rep(0, n)
  upper <- rep(NA_real_, n)
  reached <- distribution(guess, seq_len(n)) >= probability
  upper[reached] <- guess[reached]
  lower[!reached] <- guess[!reached] + 1

  # Upwards from the guess until F(upper) >= probability.
  step <- rep(1, n)

  repeat {
    open <- which(is.na(upper))

    if (length(open) == 0L) {
      break
    }

    candidate <- lower[open] + step[open] - 1
    reached <- distribution(candidate, open) >= probability[open]
    upper[open[reached]] <- candidate[reached]
    lower[open[!reached]] <- candidate[!reached] + 1
    step[open] <- 2 * step[open]
  }

  # Downwards from the upper bound until F(candidate) < probability.
  step <- rep(1, n)
  open <- which(lower < upper)

  while (length(open) > 0L) {
    candidate <- upper[open] - step[open]
    below <- candidate < lower[open]
    open <- open[!below]
    candidate <- candidate[!below]

    if (length(open) == 0L) {
      break
    }

    reached <- distribution(candidate, open) >= probability[open]
    upper[open[reached]] <- candidate[reached]
    lower[open[!reached]] <- candidate[!reached] + 1
    step[open] <- 2 * step[open]
    open <- open[reached]
  }

  # Bisection within the remaining brackets.
  while (any(lower < upper)) {
    open <- which(lower < upper)
    middle <- (lower[open] + upper[open]) %/% 2
    reached <- distribution(middle, open) >= probability[open]
    upper[open[reached]] <- middle[reached]
    lower[open[!reached]] <- middle[!reached] + 1
  }

  upper
}


# Equal-tailed prediction intervals of counts at `level`: the quantiles at
# (1 - level) / 2 and (1 + level) / 2. Because counts are discrete, they
# cover at least `level`.
.count_prediction_interval <- function(level, mean, variance, exposure,
                                       order) {
  cbind(
    lower = .count_quantile((1 - level) / 2, mean, variance, exposure, order),
    upper = .count_quantile((1 + level) / 2, mean, variance, exposure, order)
  )
}


#' Scores of count forecasts
#'
#' Scores Poisson-lognormal predictive distributions of counts, those of
#' latent Gaussian-process models with a [poisson_likelihood()], with proper
#' scores, a dispersion check, interval coverage, and randomized probability
#' integral transform (PIT) values. [gp_holdout_scores()] returns these
#' scores for latent Poisson models.
#'
#' @param observed Non-negative integer counts.
#' @param latent_mean,latent_variance Mean and variance of the latent log
#'   rate \eqn{f} at each observation, as returned by [predict_latent_gp()]
#'   (`mean` and `latent_variance`).
#' @param exposure Optional positive exposure: one value or one per count.
#' @param interval_levels Levels of the equal-tailed prediction intervals
#'   whose coverage is reported.
#' @param seed Optional seed for the randomization of the PIT values. The
#'   caller's random state is restored afterwards.
#' @param quadrature_order Number of Gauss-Hermite nodes.
#'
#' @details
#' Given \eqn{f \sim N(\mu, \sigma^2)}, a count \eqn{y \sim
#' \mathrm{Poisson}(E e^f)} has the Poisson-lognormal predictive
#' distribution, with mean \eqn{E e^{\mu + \sigma^2 / 2}} and probabilities,
#' distribution function, and quantiles computed by Gauss-Hermite
#' quadrature.
#' \describe{
#'   \item{log score}{The mean log predictive probability of the counts.
#'   Proper; higher is better.}
#'   \item{Dawid-Sebastiani score}{\eqn{(y - m)^2 / s^2 + 2 \log s} with the
#'   predictive mean \eqn{m} and standard deviation \eqn{s}; proper for the
#'   first two moments and recommended for counts by Czado, Gneiting, and
#'   Held (2009). Lower is better.}
#'   \item{dispersion}{The mean squared Pearson residual
#'   \eqn{(y - m)^2 / s^2}. It is about 1 when the predictive variance is
#'   right; clearly larger values mean that the counts are overdispersed
#'   relative to the model.}
#'   \item{coverage}{Of the equal-tailed prediction intervals at each level.
#'   Counts are discrete, so the intervals cover at least their level: their
#'   coverage is conservative, by up to the probability of their end
#'   points.}
#'   \item{randomized PIT}{\eqn{u = F(y - 1) + v\, p(y)} with
#'   \eqn{v \sim U(0, 1)}: uniform on \eqn{[0, 1]} when the forecasts are
#'   calibrated, unlike the plain PIT \eqn{F(y)} of a discrete
#'   distribution. The summary reports their mean (1/2), variance (1/12),
#'   and the p-value of a Kolmogorov-Smirnov test of uniformity.}
#' }
#'
#' @return An object of class `gaussianprocesses_count_scores` with
#' \describe{
#'   \item{`summary`}{Named numeric vector: `n`, `log_score`,
#'   `dawid_sebastiani`, `dispersion`, `pit_mean`, `pit_variance`, and
#'   `pit_ks_p_value`.}
#'   \item{`calibration`}{Data frame with one row per interval level: the
#'   `level`, empirical `coverage`, and mean `width` (upper minus lower
#'   count).}
#'   \item{`pointwise`}{Data frame with `observed`, `predictive_mean`,
#'   `predictive_variance`, `log_score`, `dawid_sebastiani`,
#'   `pearson_residual`, and `pit`.}
#' }
#' @references Czado, C., Gneiting, T., and Held, L. (2009). Predictive
#'   model assessment for count data. *Biometrics*, 65(4), 1254--1261.
#' @seealso [gp_scores()], [gp_classification_scores()],
#'   [predict_latent_gp()].
#' @examples
#' set.seed(1)
#' x <- seq(0, 10, length.out = 60)
#' counts <- rpois(60, exp(1 + sin(x)))
#' training <- seq(1, 60, by = 2)
#'
#' model <- fit_latent_gp(x[training], counts[training], rbf_kernel(0.5, 1.5),
#'   poisson_likelihood(), mean = constant_mean(1))
#' gp_holdout_scores(model, x[-training], counts[-training], seed = 1)
#' @template experimental
#' @export
gp_count_scores <- function(
    observed,
    latent_mean,
    latent_variance,
    exposure = NULL,
    interval_levels = c(0.5, 0.8, 0.9, 0.95),
    seed = NULL,
    quadrature_order = 40L) {
  observed <- .count_response(observed, "observed")
  .assert_numeric_vector(latent_mean, "latent_mean")
  .assert_numeric_vector(latent_variance, "latent_variance")
  quadrature_order <- .assert_positive_integer(
    quadrature_order,
    "quadrature_order"
  )
  .validate_seed(seed)

  if (length(latent_mean) != length(observed) ||
      length(latent_variance) != length(observed)) {
    .gp_abort(
      "'observed', 'latent_mean', and 'latent_variance' must have the same length.",
      "gaussianprocesses_validation_error"
    )
  }

  if (any(latent_variance < 0)) {
    .gp_abort(
      "'latent_variance' must be non-negative.",
      "gaussianprocesses_validation_error"
    )
  }

  interval_levels <- .validate_interval_levels(interval_levels)
  exposure <- .resolve_exposure(exposure, length(observed))
  likelihood <- poisson_likelihood()
  predictive <- likelihood_predictive(
    likelihood,
    latent_mean,
    latent_variance,
    observed = observed,
    exposure = exposure,
    quadrature_order = quadrature_order
  )
  log_score <- predictive$log_density
  predictive_sd <- sqrt(predictive$variance)
  pearson_residual <- (observed - predictive$mean) / predictive_sd
  dawid_sebastiani <- pearson_residual^2 + 2 * log(predictive_sd)

  below <- .count_distribution(
    observed - 1,
    latent_mean,
    latent_variance,
    exposure,
    quadrature_order
  )
  uniform <- .with_local_seed(
    seed,
    function() stats::runif(length(observed))
  )
  pit <- pmin(pmax(below + uniform * exp(log_score), 0), 1)

  calibration <- do.call(
    rbind,
    lapply(
      interval_levels,
      function(level) {
        interval <- .count_prediction_interval(
          level,
          latent_mean,
          latent_variance,
          exposure,
          quadrature_order
        )

        data.frame(
          level = level,
          coverage = mean(
            observed >= interval[, "lower"] & observed <= interval[, "upper"]
          ),
          width = mean(interval[, "upper"] - interval[, "lower"])
        )
      }
    )
  )

  structure(
    list(
      summary = c(
        n = length(observed),
        log_score = mean(log_score),
        dawid_sebastiani = mean(dawid_sebastiani),
        dispersion = mean(pearson_residual^2),
        pit_mean = mean(pit),
        pit_variance = if (length(pit) > 1L) stats::var(pit) else NA_real_,
        # Randomized PIT values are continuous, but extreme counts can give
        # values of exactly 0 or 1; ks.test() then warns about ties.
        pit_ks_p_value = if (length(pit) > 1L) {
          suppressWarnings(stats::ks.test(pit, "punif"))$p.value
        } else {
          NA_real_
        }
      ),
      calibration = calibration,
      pointwise = data.frame(
        observed = observed,
        predictive_mean = predictive$mean,
        predictive_variance = predictive$variance,
        log_score = log_score,
        dawid_sebastiani = dawid_sebastiani,
        pearson_residual = pearson_residual,
        pit = pit
      )
    ),
    class = "gaussianprocesses_count_scores"
  )
}


#' @export
print.gaussianprocesses_count_scores <- function(x, ...) {
  summary <- x$summary
  number <- function(value) format(signif(value, 4), trim = TRUE)

  cat("Count scores for ", summary[["n"]], " observations\n", sep = "")
  cat("  log score (mean): ", number(summary[["log_score"]]), "\n", sep = "")
  cat(
    "  Dawid-Sebastiani score (mean): ",
    number(summary[["dawid_sebastiani"]]),
    "\n",
    sep = ""
  )
  cat(
    "  dispersion: ", number(summary[["dispersion"]]),
    " (about 1 if the predictive variance is right)\n",
    sep = ""
  )
  cat(
    "  randomized PIT mean and variance: ", number(summary[["pit_mean"]]),
    ", ", number(summary[["pit_variance"]]),
    " (0.5 and 0.0833 if calibrated)\n",
    sep = ""
  )
  cat(
    "  coverage: ",
    paste0(
      format(100 * x$calibration$level, trim = TRUE), "% ",
      format(round(100 * x$calibration$coverage, 1), trim = TRUE), "%",
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
  invisible(x)
}
