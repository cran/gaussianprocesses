# Scores of binary probability forecasts --------------------------------------
#
# For responses y_i in {0, 1} and forecast probabilities p_i = Pr(y_i = 1):
#   log loss     -mean(log p(y_i)), with p(1) = p_i and p(0) = 1 - p_i;
#   Brier score  mean((p_i - y_i)^2);
#   ROC AUC      Pr(p_i > p_j) + Pr(p_i = p_j) / 2 for a random event i and
#                non-event j, from mid-ranks (the Mann-Whitney statistic);
#   reliability  observed frequency against mean probability in equal-width
#                bins of p, and the expected calibration error
#                ECE = sum_b (n_b / n) |mean(p_b) - mean(y_b)|.


# Mann-Whitney estimate of the area under the ROC curve, NA without both
# classes.
.roc_auc <- function(observed, probability) {
  n_events <- sum(observed == 1)
  n_non_events <- length(observed) - n_events

  if (n_events == 0L || n_non_events == 0L) {
    return(NA_real_)
  }

  ranks <- rank(probability, ties.method = "average")

  (sum(ranks[observed == 1]) - n_events * (n_events + 1) / 2) /
    (n_events * n_non_events)
}


.reliability_bins <- function(observed, probability, n_bins) {
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  bin <- findInterval(probability, breaks, rightmost.closed = TRUE,
                      all.inside = TRUE)
  occupied <- sort(unique(bin))

  data.frame(
    bin = occupied,
    lower = breaks[occupied],
    upper = breaks[occupied + 1L],
    n = vapply(occupied, function(b) sum(bin == b), integer(1)),
    mean_probability = vapply(
      occupied,
      function(b) mean(probability[bin == b]),
      numeric(1)
    ),
    observed_frequency = vapply(
      occupied,
      function(b) mean(observed[bin == b]),
      numeric(1)
    )
  )
}


#' Scores of binary probability forecasts
#'
#' Scores predicted class probabilities against binary responses with proper
#' scoring rules (log loss and Brier score), discrimination (ROC AUC), and
#' calibration (reliability diagram data and the expected calibration
#' error). [gp_holdout_scores()] returns these scores for latent models with
#' a Bernoulli likelihood.
#'
#' @param observed Binary responses: 0/1 numbers, logical values, or a
#'   two-level factor whose second level is the event coded 1, as in
#'   [bernoulli_likelihood()].
#' @param probability Forecast probabilities \eqn{p_i = \Pr(y_i = 1)}, one
#'   per response.
#' @param n_bins Number of equal-width probability bins of the reliability
#'   diagram and the expected calibration error.
#' @param log_density Optional log predictive probabilities of the observed
#'   responses, for the log loss. They default to \eqn{\log p_i} or
#'   \eqn{\log(1 - p_i)}; [gp_holdout_scores()] supplies them from the
#'   predictive distribution, which keeps their accuracy when \eqn{p_i} is
#'   close to 0 or 1.
#'
#' @details
#' \describe{
#'   \item{log loss}{\eqn{-\frac1n \sum_i \log p(y_i)}, the negative mean log
#'   predictive probability. Proper; lower is better; \eqn{\log 2 \approx
#'   0.693} for \eqn{p_i = 1/2}.}
#'   \item{Brier score}{\eqn{\frac1n \sum_i (p_i - y_i)^2}. Proper; lower is
#'   better; 0.25 for \eqn{p_i = 1/2}.}
#'   \item{ROC AUC}{The probability that a random event receives a higher
#'   forecast than a random non-event, with ties counted as one half: the
#'   Mann-Whitney statistic from mid-ranks. 0.5 for uninformative forecasts;
#'   `NA` without both classes. It measures discrimination only and ignores
#'   calibration.}
#'   \item{reliability}{The forecasts are grouped into `n_bins` equal-width
#'   bins of \eqn{[0, 1]}. For calibrated forecasts the observed frequency
#'   of each bin is close to its mean probability, within about
#'   `standard_error` \eqn{= \sqrt{\bar p (1 - \bar p) / n_b}}.}
#'   \item{ECE}{The expected calibration error
#'   \eqn{\sum_b \frac{n_b}{n} |\bar p_b - \bar y_b|}. It depends on the
#'   binning, and even calibrated forecasts have a positive ECE of the order
#'   of the bin standard errors.}
#' }
#'
#' @return An object of class `gaussianprocesses_classification_scores` with
#' \describe{
#'   \item{`summary`}{Named numeric vector: `n`, `log_loss`, `brier_score`,
#'   `auc`, `ece`, `base_rate` (the mean response), and `mean_probability`.}
#'   \item{`reliability`}{Data frame with one row per non-empty bin: `bin`,
#'   `lower`, `upper`, `n`, `mean_probability`, `observed_frequency`, and
#'   `standard_error`.}
#'   \item{`pointwise`}{Data frame with `observed` (0/1), `probability`,
#'   `log_density`, and `squared_error`.}
#' }
#' @references Gneiting, T. and Raftery, A. E. (2007). Strictly proper
#'   scoring rules, prediction, and estimation. *Journal of the American
#'   Statistical Association*, 102(477), 359--378.
#' @seealso [gp_scores()] for Gaussian predictions, [predict_latent_gp()].
#' @examples
#' set.seed(1)
#' x <- runif(80, 0, 6)
#' y <- rbinom(80, 1, plogis(2 * sin(x)))
#' training <- 1:50
#'
#' model <- fit_latent_gp(x[training], y[training], rbf_kernel(4, 1),
#'   bernoulli_likelihood())
#' scores <- gp_holdout_scores(model, x[-training], y[-training], n_bins = 5)
#' scores
#' scores$reliability
#' @template experimental
#' @export
gp_classification_scores <- function(
    observed,
    probability,
    n_bins = 10L,
    log_density = NULL) {
  observed <- .binary_response(observed, "observed")
  .assert_numeric_vector(probability, "probability")
  n_bins <- .assert_positive_integer(n_bins, "n_bins")

  if (length(probability) != length(observed)) {
    .gp_abort(
      "'observed' and 'probability' must have the same length.",
      "gaussianprocesses_validation_error"
    )
  }

  if (any(probability < 0 | probability > 1)) {
    .gp_abort(
      "'probability' must lie in [0, 1].",
      "gaussianprocesses_validation_error"
    )
  }

  if (is.null(log_density)) {
    log_density <- ifelse(observed == 1, log(probability), log1p(-probability))
  } else {
    .assert_numeric_vector(log_density, "log_density", finite = FALSE)

    if (length(log_density) != length(observed) || anyNA(log_density) ||
        any(log_density > 0)) {
      .gp_abort(
        paste0(
          "'log_density' must contain one non-positive log probability per ",
          "response."
        ),
        "gaussianprocesses_validation_error"
      )
    }
  }

  reliability <- .reliability_bins(observed, probability, n_bins)
  reliability$standard_error <- sqrt(
    reliability$mean_probability * (1 - reliability$mean_probability) /
      reliability$n
  )
  squared_error <- (probability - observed)^2

  structure(
    list(
      summary = c(
        n = length(observed),
        log_loss = -mean(log_density),
        brier_score = mean(squared_error),
        auc = .roc_auc(observed, probability),
        ece = sum(
          reliability$n *
            abs(reliability$mean_probability - reliability$observed_frequency)
        ) / length(observed),
        base_rate = mean(observed),
        mean_probability = mean(probability)
      ),
      reliability = reliability,
      pointwise = data.frame(
        observed = observed,
        probability = probability,
        log_density = log_density,
        squared_error = squared_error
      )
    ),
    class = "gaussianprocesses_classification_scores"
  )
}


#' @export
print.gaussianprocesses_classification_scores <- function(x, ...) {
  summary <- x$summary
  number <- function(value) format(signif(value, 4), trim = TRUE)

  cat(
    "Classification scores for ", summary[["n"]], " responses (base rate ",
    number(summary[["base_rate"]]), ")\n",
    sep = ""
  )
  cat("  log loss: ", number(summary[["log_loss"]]), "\n", sep = "")
  cat("  Brier score: ", number(summary[["brier_score"]]), "\n", sep = "")
  cat("  ROC AUC: ", number(summary[["auc"]]), "\n", sep = "")
  cat(
    "  expected calibration error: ", number(summary[["ece"]]), " (",
    nrow(x$reliability), " occupied bins)\n",
    sep = ""
  )
  invisible(x)
}
