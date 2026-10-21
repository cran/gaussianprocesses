# gaussianprocesses 1.0.0

This is the first release on CRAN and the first since 0.1.0.

* Breaking changes. The names renamed after 0.1.0 are removed, without the
  release of deprecation warnings that NEWS promised: version 1.0.0 is the
  first release after 0.1.0, and the stable API starts without deprecated
  names (#40).
  * `fit_gp_time()` is removed; use `fit_time_series_gp()`.
  * The `tolerance` argument of `fit_gp()`, `optimize_gp()`,
    `fit_sparse_gp()`, `sample_gp_prior()`, `sample_gp_posterior()`, and
    `simulate_gp_data()` is removed; use `symmetry_tolerance`. In
    `sample_gp_posterior()`, which passes `...` to its methods, an old
    `tolerance` argument is now ignored like any other unknown argument.
  * In `fit_heteroscedastic_gp()`, `tolerance` and `numerical_tolerance` are
    removed; use `convergence_tolerance` and `symmetry_tolerance`.
  * `predict_heteroscedastic_gp()` no longer returns `noise_variance`; use
    `observation_noise_variance`.
* CRAN preparation (#40).
  * `citation("gaussianprocesses")` gives a citation from `inst/CITATION`.
  * `DESCRIPTION` has no `BugReports` field: GitLab serves the issue list
    only to logged-in users, and R's incoming check flags the public work
    items page. The README links the tracker.
  * CI runs R 4.6.1, adds a manual `cran-check` job that installs every
    suggested package and checks the built tarball with its tests, and a
    manual `text-checks` job for spelling and URLs. RELEASING.md describes
    CRAN submission and resubmission.

* Stability tiers and the versioning policy for 1.0.
  * Every exported function is stable, experimental, or deprecated. Its help
    page says which in a new "Stability" section, and the reference index is
    grouped by tier. `?gaussianprocesses` sets out the policy. From 1.0.0,
    stable interfaces change incompatibly only in a major release, after a
    deprecation period, and experimental interfaces may change in a minor
    release, with a NEWS entry.
  * Experimental:
    * the sparse approximations;
    * latent models with their likelihoods and scores;
    * the heteroscedastic method;
    * state-space inference (`optimize_time_series_gp()` and
      `fit_time_series_gp(method = "state_space")`);
    * the spectral-mixture and changepoint kernels;
    * derivative observations (`fit_gp(derivative = )`) and
      `predict_gradient_gp()`;
    * hyperparameter uncertainty;
    * the numerical and calibration diagnostics;
    * the benchmark and simulation-scenario helpers.

    Everything else is stable.
  * Schema 4 is the model format of version 1.0, and every 1.x release will
    read models saved by any earlier version. The tests now read models
    saved by the 0.1.0 release and a set of schema-4 models, and check that
    they give the results they were saved with.
  * CONTRIBUTING.md records the versioning and deprecation policy and two
    decisions for 1.0: the kernel registry stays internal, and the package
    stays pure R, with the measurements behind both.

* Validation against independent references.
  * `inst/validation/references.R` defines 35 references, at least one for
    every inference method: other implementations (DiceKriging, kernlab,
    nlme, gplite, and LAPACK's pivoted Cholesky decomposition), brute-force
    computations (dense linear algebra, refits, numerical integration and
    differentiation, `optim()`), and closed forms. Each has a tolerance and,
    for other packages, the version it was set with. A test runs every
    reference, skipping those whose package is not installed, and checks
    that every exported function belongs to a validated method or is listed
    as performing no inference.
  * The DiceKriging and kernlab comparison scripts in `inst/benchmarks` and
    the gplite, nlme, and DiceKriging comparison tests are replaced by
    references. CI now installs DiceKriging and kernlab, so every reference
    runs there.
  * A new validation article tabulates the references as it is built and
    reproduces three analyses on data shipped with R: Rasmussen and
    Williams's Mauna Loa CO2 kernel on `datasets::co2`, a classifier for
    versicolor and virginica irises against logistic regression, and a
    Poisson model of `datasets::discoveries` against a quadratic Poisson
    regression.

* State-space (Kalman) inference for time series.
  * `fit_time_series_gp(method = "state_space")` computes the exact GP
    posterior and log marginal likelihood in O(n) time by Kalman filtering
    and Rauch-Tung-Striebel smoothing, for Matérn-1/2, 3/2, and 5/2
    kernels, their sums, and scaled versions of them. The kernel registry
    supplies their stochastic differential equations, and the transition
    matrices are in closed form. Results agree with exact inference to a
    relative 1e-8 on irregular and tied inputs. Other kernels raise an
    error of class `gaussianprocesses_state_space_error`.
  * `forecast_gp()`, `predict()`, `logLik()`, `fitted()`, `residuals()`,
    `summary()`, `log_marginal_likelihood()`, and
    `rolling_origin_gp(method = "state_space")` accept the new models, and
    `sample_gp_posterior()` draws by forward filtering, backward sampling.
  * New `optimize_time_series_gp()` estimates the hyperparameters of a
    time-series GP with either method; the state-space likelihood is
    optimized with finite-difference gradients.
  * `inst/benchmarks/state-space-scaling.R` times the method from 1000 to
    100 000 observations, and the time-series vignette compares state-space
    and exact forecasts.

* Inducing-point selection and optimization for sparse Gaussian processes.
  * `select_inducing_points(method = "variance", kernel = ...)` performs
    greedy conditional-variance selection, the pivot order of a pivoted
    Cholesky decomposition of K(X, X), in O(n m^2) time without forming
    K(X, X). The sparse fitters accept `selection = "variance"`.
  * `optimize_sparse_gp(optimize_inducing = TRUE)` optimizes the inducing
    inputs jointly with the hyperparameters, with analytical gradients from
    the kernel input derivatives. It is recommended for VFE; with FITC it
    overfits.
  * Sparse fits detect inducing points that are nearly identical under the
    kernel, record them in `$inducing_collapse`, and warn (class
    `gaussianprocesses_inducing_collapse_warning`).
  * The research note on sparse approximations reports the gap to the exact
    marginal likelihood that optimizing the inducing points closes on the
    simulation scenarios, reproduces FITC's noise underestimation with free
    inducing points, and times the greedy selection.

* Variational sparse Gaussian processes, gradients, and sparse sampling.
  * `fit_sparse_gp()` gains `method = "vfe"`, the collapsed variational
    bound of Titsias (2009), which lies below the exact log marginal
    likelihood and equals it when every input is an inducing point. VFE
    predicts with the deterministic training conditional. FITC stays the
    default.
  * Sparse models report the trace term tr(K_ff - Q_ff) / sigma^2, a
    diagnostic of how well the inducing points cover the data, and record
    their `method`; models now use schema 4, and older sparse models are read
    as FITC.
  * New `optimize_sparse_gp()` optimizes kernel parameters and noise for VFE
    (the default) or FITC with fixed inducing points, using the optimizer of
    `optimize_gp()` and analytical gradients computed in O(n m^2) time
    without n x n matrices.
  * `sample_gp_posterior()` draws from FITC and VFE posteriors, and
    `predict_sparse_gp()` predicts marginal variances in blocks
    (`block_size`), with memory linear in the number of prediction inputs.
  * A research note compares VFE and FITC on the failure modes reported by
    Bauer, van der Wilk, and Rasmussen (2016), and reports time and memory
    for up to 50 000 observations and 500 inducing points.

* Poisson count models with latent Gaussian processes.
  * Latent models accept mean coefficients estimated by maximum likelihood
    (`estimate_coefficients("ml")`, the default of `linear_mean()`):
    `fit_latent_gp()` maximizes the approximate marginal likelihood over
    them with its analytical gradient, and `optimize_latent_gp()` estimates
    them jointly with the hyperparameters as `mean.<name>`. A constant mean
    estimated this way is the log level of a count model. `coef()` returns
    the coefficients of latent models.
  * `predict_latent_gp()` returns, for Poisson likelihoods, `rate_interval`,
    the exact interval of the lognormal rate, and `prediction_interval`,
    equal-tailed intervals for new counts from the Poisson-lognormal
    predictive distribution, which are conservative because counts are
    discrete; for Gaussian likelihoods it returns the `prediction_interval`
    of the exact model. `sample_gp_posterior()` returns `rate` draws.
  * The Poisson-lognormal distribution function is computed by quadrature in
    the latent value or, where the Poisson probabilities change too sharply
    for that, through the gamma representation of the Poisson distribution
    function; it agrees with numerical integration to 1e-9 for rates from
    0.3 to 5000.
  * New `gp_count_scores()` scores count forecasts with the log score,
    Dawid-Sebastiani score, Pearson dispersion (an overdispersion check),
    coverage of the count intervals, and randomized PIT values;
    `gp_holdout_scores()` returns them for latent Poisson models.
  * A new article analyses `datasets::discoveries` and reports a
    simulation: rate intervals cover the true rate at their nominal level
    with known hyperparameters, count intervals are conservative, randomized
    PIT values are uniform, and the dispersion check detects
    overdispersion.
  * The Poisson fits agree with `gplite`.

* Binary classification with latent Gaussian processes.
  * New `gp_classification_scores()` scores probability forecasts with the
    log loss, Brier score, ROC AUC, reliability bins, and expected
    calibration error. `gp_holdout_scores()` returns these scores for latent
    models with a Bernoulli likelihood, with the log loss from predictive
    log probabilities, and the scores of the Gaussian predictive for a
    Gaussian likelihood.
  * `predict_latent_gp()` returns `probability_interval`, the interval of
    the class probability itself, and `sample_gp_posterior()` returns
    `probability` draws for Bernoulli likelihoods.
  * `optimize_latent_gp()` warns (class `gaussianprocesses_bound_warning`)
    when a parameter ends at a bound of the search and records it in
    `$optimization$at_bound`. With separable classes the signal variance
    grows to a bound or to a very large value.
  * A new research note simulates 300 data sets from a known GP classifier.
    With the true hyperparameters, the Laplace predictive probabilities are
    slightly underconfident (pooled reliability deviations of up to 7
    standard errors), while their Brier score exceeds the oracle's by only
    0.015; the note explains why and documents the plug-in and MacKay
    shortcuts.
  * The Laplace fits agree with the `gplite` package (now suggested) to
    1e-6 in the mode, approximate marginal likelihood, and latent
    predictions.

* Latent Gaussian-process models with the Laplace approximation, for
  responses observed through a likelihood such as binary labels or counts.
  * New `fit_latent_gp(x, y, kernel, likelihood, ...)` finds the posterior
    mode by Newton's method on the well-conditioned matrix
    B = I + W^(1/2) K W^(1/2) (Rasmussen and Williams, Algorithm 3.1), with
    the package's jitter policy, a step-halving line search, and a
    convergence test on the Newton decrement. `$laplace` reports the
    iterations, step halvings, decrement, and gradient norm, and a warning
    is raised if the iterations do not converge. Likelihoods that are not
    log-concave are refused. The Newton update is written so that it stays
    accurate for counts in the thousands, where the usual form loses all
    digits.
  * `log_marginal_likelihood()` and `log_marginal_likelihood_gradient()`
    return the Laplace approximation and its analytical gradient, including
    the implicit dependence of the mode on the hyperparameters, for kernel
    and likelihood parameters.
  * New `optimize_latent_gp()` optimizes the hyperparameters with the
    optimizer of `optimize_gp()`, which now shares its setup, multi-start,
    and diagnostics code, and warm-starts each Newton iteration from the
    previous mode.
  * New `predict_latent_gp()` gives the latent posterior (Algorithm 3.2) and
    the predictive distribution of new responses through the likelihood,
    including class probabilities. `sample_gp_posterior()` draws latent
    values from the Gaussian approximation and observations from the
    likelihood. `predict()`, `logLik()`, `nobs()`, `fitted()`, and
    `summary()` have methods for latent models.
  * Means may have fixed coefficients or a Gaussian coefficient prior.
  * With `gaussian_likelihood()` the approximation is exact, and the model
    equals `fit_gp()` to 1e-10. `?fit_latent_gp` reports the error of the
    approximation against tensor-product quadrature for up to three
    observations.

* Likelihoods for non-Gaussian responses, the foundation for latent
  Gaussian-process models.
  * New `gaussian_likelihood()`, `bernoulli_likelihood()` (logit and probit
    links), and `poisson_likelihood()` (log link, with exposure) build
    likelihood specifications: plain data with a print method, like kernel
    specifications. Binary responses may be 0/1, logical, or a two-level
    factor whose second level is the event, as in `glm()`.
  * New `evaluate_likelihood()` returns the log density and its first three
    derivatives in the latent values. They stay finite and accurate for
    latent values up to 40 in magnitude: the probit derivatives use the
    inverse Mills ratio in log space, and its continued fraction in the far
    tail.
  * New `likelihood_predictive()` gives the predictive mean, variance, class
    probability, and log density of new responses when the latent value is
    Gaussian, in closed form where one exists and otherwise by adaptive
    Gauss-Hermite quadrature centred at the mode of the integrand.
  * Gauss-Hermite rules come from the Golub-Welsch eigenvalue method, with
    weights from Christoffel sums that keep the smallest weights accurate.
  * Non-Gaussian models will be fitted through one entry point that takes a
    likelihood; `fit_gp()` remains the exact Gaussian model (see the API
    audit).

* Proper scoring rules, hyperparameter uncertainty, and multi-start
  diagnostics.
  * New `gp_scores()` scores Gaussian predictions with the log score, the
    CRPS in closed form, interval scores and coverage across a grid of
    levels, PIT values, SMSE, and MSLL. `gp_holdout_scores()` scores any
    model on held-out data, and `gp_loo_scores()` on its exact
    leave-one-out predictions.
  * New `gp_hyperparameter_uncertainty()` computes the observed information
    at the estimate from differences of the analytical gradient on the
    optimizer scale, with standard errors, correlations, intervals, and a
    flag for poorly identified hyperparameters. It detects the
    variance-length-scale ridge of a Matérn process observed within one
    length scale, and stays silent when the data span twenty length scales.
    New `gp_profile_likelihood()` profiles one hyperparameter.
  * `optimize_gp()` records the distinct optima its starts reached in
    `$optimization$optima`, and `summary()` reports them.
  * A new research note checks the standard errors in a 200-replicate
    simulation: for the noise variance they match the spread of the
    estimates, and for the variance and length scale they are 13% and 20%
    too small at 60 observations, with 95% intervals covering about 90% of
    the time.

* New `spectral_mixture_kernel()` (Wilson and Adams, 2013): a sum of
  components w exp(-2 pi^2 tau^2 v) cos(2 pi tau mu), one weight, frequency,
  and spectral variance per component and input dimension. A component with
  frequency 0 is an RBF kernel.
  * Each component is a kernel registry entry, so hyperparameter gradients,
    input derivatives, and derivative observations work for it.
    Frequencies are real-valued (identity coordinate).
  * New `initialize_spectral_mixture()` places the components at the peaks
    of the data's spectrum: the periodogram for regular inputs and a
    Lomb-Scargle periodogram for irregular ones, with a mixture of Gaussians
    and a uniform background, for the noise floor, fitted by EM. On noisy
    sums of two sinusoids it recovers both frequencies within a fraction of
    a frequency-resolution bin, for regular and irregular sampling.
  * `optimize_gp()` bounds every frequency by 0 and the Nyquist frequency of
    its input column, and its later starts displace frequencies by one
    frequency-resolution bin rather than one unit.
  * A new worked example compares a spectral mixture with a composite
    kernel on `datasets::co2`, by held-out log predictive density.

* New `changepoint_kernel(before, after, location, steepness, column)`
  switches between two kernels around an estimable location along one input
  column: (1 - s(x)) k1 (1 - s(x')) + s(x) k2 s(x') with a sigmoid s.
  * It is the first composite with a real-valued hyperparameter: the
    location is optimized on the identity scale and the steepness on the
    log scale. Parameter paths are `location`, `steepness`, `before.`, and
    `after.`, and nesting gives several changes.
  * Hyperparameter gradients, input derivatives, and derivative
    observations support it, through the chain and product rules.
  * `optimize_gp()` bounds every changepoint location by the training range
    of its column and starts it at evenly spaced quantiles of that column,
    because the likelihood is often multimodal in the location.
  * `time_series_kernel()` gains `changepoint`, which wraps the
    trend-periodic-local kernel in a changepoint with independent copies
    before and after the break.
  * The time-series vignette gains a structural-break example: on a
    simulated series, the estimated break is one day from the true one, and
    the changepoint model keeps forecasting the new level where a single
    kernel falls back towards zero.

* Derivative observations. `fit_gp()` and `optimize_gp()` take
  `derivative`, the type of each observation: 0 for a function value and d
  for the derivative with respect to input d.
  * The joint covariance of values and derivatives is built from the
    kernel's input derivatives, and every observation has its own noise
    variance. Prediction of values (`predict_gp()`) and derivatives
    (`predict_gradient_gp()`), posterior draws, leave-one-out diagnostics,
    and the likelihood and its gradient use it.
  * The likelihood gradient uses closed-form hyperparameter derivatives of
    the input derivatives: a new kernel registry field,
    `input_parameter_gradient`, for every differentiable kernel, combined
    for sums, products, scaled kernels, and selections of input columns.
  * Parametric means contribute their input derivatives. Kernels that are
    not mean-square differentiable raise
    `gaussianprocesses_smoothness_error`.
  * `gp_numerical_diagnostics()` reports conditioning per observation type
    in `observation_types`.
  * Models now use schema 3, which adds `derivative`; older models are read
    as before.
  * The kernel-derivative vignette derives the joint model and shows the
    variance reduction from gradient observations.

* Parametric mean functions — universal kriging.
  * `linear_mean()`, `polynomial_mean()`, and `basis_mean()` define trends
    m(x) = h(x)'β, and `constant_mean()` accepts an unknown constant.
    Polynomial bases standardize their inputs to [-1, 1] for conditioning.
  * The coefficients are fixed (a numeric vector), estimated by generalized
    least squares (`estimate_coefficients()`, with the profile likelihood
    or REML), or marginalized under a Gaussian prior or its vague limit
    (`coefficient_prior()`). Marginalized coefficients add their
    uncertainty to predictions, posterior draws, and leave-one-out
    diagnostics, and a Gaussian prior is included in prior draws.
  * `log_marginal_likelihood()`, its gradient, and `optimize_gp()` use the
    profile, restricted, or marginal likelihood that matches the mean.
    Coefficients are profiled out in closed form, and
    `$optimization$likelihood` and `$optimization$mean_coefficients`
    record the result. `coef()` and `vcov()` return the coefficients and
    their covariance. Rank-deficient bases raise
    `gaussianprocesses_rank_deficiency_error`.
  * Tests check the estimates and likelihoods against dense formulas,
    `nlme::gls()` (ML and REML), and `DiceKriging::km()` (estimates, profile
    likelihood, and universal kriging), and the Gaussian-prior model against
    the Gaussian process with the augmented kernel.
  * Sparse and heteroscedastic models accept parametric means with fixed
    coefficients only.
  * Models now use schema 2, which adds `mean_fit`; older models are read
    as before.
  * The marginal-likelihood vignette gains a universal-kriging section.
* `inst/benchmarks/external-dicekriging.R` used DiceKriging length scales
  √2 times too long. DiceKriging's Gaussian covariance has the same form
  as the RBF kernel, so the length scales are equal.

* New `select_dimensions()` restricts a kernel to some input columns,
  k(x_S, x'_S). Sums of restricted kernels give additive models such as
  f(x) = f1(x1) + f2(x2, x3).
  * The selection is not a hyperparameter and adds no level to parameter
    paths. ARD parameters have one value per selected column.
  * Evaluation, diagonals, hyperparameter and input derivatives, fitting,
    optimization, prediction, sampling, sparse FITC and heteroscedastic
    models all accept restricted kernels. Derivatives along unselected
    columns are zero.
  * Fitting and sampling functions now check, before computing any
    covariance, that selected columns exist and that every ARD parameter
    has one value or one per input dimension of its kernel. The error
    names the parameter path, such as `kernel2.length_scale`.
  * The ARD vignette gains an additive model on simulated data.

* New `kernel_input_gradient()` and `predict_gradient_gp()`.
  * `kernel_input_gradient()` differentiates a kernel with respect to its
    inputs: the covariance between the gradient of the process and the
    process.
  * `predict_gradient_gp()` returns the posterior mean, variance, and
    credible intervals of each partial derivative of the latent function
    of an exact model, computed in blocks like `predict_gp()`.
  * Every differentiable built-in kernel has closed-form input
    derivatives that need no limit at coincident inputs, and sums,
    products, and scaled kernels follow the sum and product rules.
    Kernels whose process is not mean-square differentiable, Matérn 1/2
    and white noise and any composite that contains them, raise an error
    of class `gaussianprocesses_smoothness_error`.
  * The kernel-derivative vignette derives the formulas and tabulates the
    smoothness of every kernel. Kernel registry entries declare a
    `smoothness` and their input derivatives (see CONTRIBUTING.md).

* `gp_diagnostics()` computes leave-one-out once instead of twice, which
  removes about a third of its run time (1.11 s to 0.78 s at n = 1500).
* Test coverage is 96% of lines, up from 92%, and every file in `R/` is
  above 89%. A `coverage` CI job, which runs in scheduled pipelines or when
  started by hand, measures it with covr and fails if a file falls below
  85%. The built-in kernels are now registered on first use, so coverage
  tools see their code; their results are unchanged.

* **Behaviour change:** `fit_heteroscedastic_gp()` corrects its noise
  targets for the leave-one-out uncertainty of the latent function. Squared
  leave-one-out residuals estimate the noise variance plus that
  uncertainty, so the noise was overestimated. In a 50-replicate simulation
  study the median log noise error was +0.24, +0.14, and +0.08 for n = 40,
  100, and 200. It is now -0.06, +0.02, and +0.02, each within one Monte
  Carlo standard error of zero, with a lower root-mean-square error at every
  n. `inst/benchmarks/heteroscedastic-noise.R` holds the study, and the
  heteroscedastic research note derives the correction.
* **Behaviour change:** `fit_heteroscedastic_gp()` defaults to
  `max_iterations = 50` instead of 6. In the study no fit converged within 6
  iterations, and every fit that converged did so within 50. It now warns,
  with class `gaussianprocesses_convergence_warning`, when it stops before
  converging; about 2 to 10% of the simulated fits cycle without
  converging.

* Observation noise is handled consistently.
  * `optimize_gp()` accepts one noise variance per observation. With
    `optimize_noise = FALSE` the vector is fixed. With
    `optimize_noise = TRUE` it is a known relative noise pattern whose
    overall scale is estimated with an analytical gradient, and
    `$optimization$noise_parameterization` records which case applied.
  * `rolling_origin_gp()` accepts one noise variance per observation and
    gives each training window and each forecast the values of its rows.
  * `benchmark_sparse_gp()` gives both models the same scalar or
    per-observation noise and takes prediction-time noise through the new
    `observation_noise_variance` argument.
  * Noise of the wrong length raises a validation error.

* API audit (`inst/notes/api-audit.md` records the decision for every
  exported function). Renamed functions and arguments keep working for at
  least one minor release, with a deprecation warning:
  * `fit_gp_time()` is now `fit_time_series_gp()`.
  * `tolerance`, which meant different things in different functions, is
    now `symmetry_tolerance` in `fit_gp()`, `optimize_gp()`,
    `fit_sparse_gp()`, `sample_gp_prior()`, `sample_gp_posterior()`, and
    `simulate_gp_data()`.
  * In `fit_heteroscedastic_gp()`, `tolerance` is now
    `convergence_tolerance` and `numerical_tolerance` is now
    `symmetry_tolerance`.
* Predictions share one set of fields with the same meaning, documented in
  `?predict_gp`. `predict_heteroscedastic_gp()` now returns
  `observation_noise_variance` and, with `include_covariance = TRUE`,
  `observation_covariance`; its `noise_variance` field is a deprecated copy
  of `observation_noise_variance`.
* Exact, sparse, and heteroscedastic models have methods for `predict()`,
  `logLik()`, `nobs()`, `fitted()`, `residuals()`, and `summary()`, so
  `AIC()` and `BIC()` work too (see `?gp_model_methods`). `logLik()` is the
  type-II marginal likelihood; comparing kernels by AIC or BIC is a
  heuristic. Heteroscedastic models have no `logLik()`, because their noise
  function is a plug-in estimate. The package now declares `stats` in
  `Imports`.
* The documentation explains the two kernel families: specifications such as
  `rbf_kernel()` for models, and functions such as `kernel_rbf()` that return
  one covariance matrix.
* The project is public: the README no longer asks for an access token, and
  `BugReports` points to the issue tracker page that GitLab serves to
  visitors who are not signed in.

* **Behaviour change:** numerical jitter is now relative to the scale of the
  covariance being factorized. When a Cholesky factorization fails, the
  diagonal jitter is the jitter ladder (`initial_jitter`, `fallback_jitter`,
  `jitter_multiplier`) times the mean of the covariance's diagonal, instead
  of the ladder itself. Results therefore no longer depend on the units of
  `y`: rescaling `y` by 1e-4 or 1e4, which previously changed the log
  marginal likelihood, posterior mean, and function draws by up to 170%,
  now changes them by at most about 1e-6 (relative, after undoing the
  scaling). Fits that need no jitter, which is most of them, are unchanged;
  fits that need jitter on covariances whose mean diagonal is not 1 now use
  a different jitter. This applies to every factorization: exact and
  sparse fits, optimization, sampling, and simulation.
  `inst/benchmarks/jitter-scaling.R` compares the two policies.
* Jitter is reported consistently. `gp_numerical_diagnostics()` returns
  `relative_jitter` and `covariance_scale` (the mean diagonal of the
  training covariance) next to `numerical_jitter`; they replace
  `jitter_relative_to_covariance_scale`, which divided by the larger of 1
  and the largest diagonal entry. Sample objects record `relative_jitter`
  and `covariance_scale`, and models, diagnostics, and samples print the
  jitter with its size relative to the covariance scale.
* The log marginal likelihood gradient keeps forming the inverse covariance
  from its Cholesky factor with `chol2inv()`. Against an exact reference it was as
  accurate as two triangular solves per parameter for condition numbers up
  to 1e12, and 4 to 18 times faster. `inst/benchmarks/gradient-numerics.R`
  holds the study, and `vignette("v04-numerical-stability")` summarizes both
  studies.
* Fitted exact models no longer store the kernel matrix and the
  observation covariance, only the Cholesky factor of the covariance, which
  is all that prediction, likelihoods, gradients, leave-one-out diagnostics,
  and sampling need. A model with 2000 observations shrinks from 91.6 MB to
  30.6 MB. The matrices are recomputed exactly when needed. Code that read
  `model$kernel_covariance` or `model$observation_covariance` should use
  `evaluate_kernel(model$kernel, model$x)`, plus
  `diag(model$training_noise_variance)` for the observation covariance.
* Every fitted model now records `schema_version` and `package_version`,
  and exact models record the numerical settings used to fit them. Every
  function that accepts a model converts models saved by earlier versions,
  which give the same results, and refuses models with an unknown schema
  with a `gaussianprocesses_model_format_error`. The new help page
  `?gaussianprocesses_model` documents the fields of fitted models, their
  schema, and how to save them.

* `predict_gp()` no longer forms the covariance between prediction inputs
  when only marginal predictions are requested (the default,
  `include_covariance = FALSE`). It evaluates the prior variances with
  `kernel_diagonal()` and the training-by-prediction cross-covariance in
  blocks, so memory no longer grows with the square of the number of
  prediction inputs. With 500 training inputs, 8000 prediction inputs
  previously needed about 2 GB and now need about 65 MB; 50000 now need
  about 100 MB. Results agree with the full-covariance path to about 1e-14,
  and `include_covariance = TRUE` is unchanged. The new `block_size`
  argument sets the block size; by default it is chosen from a memory
  budget. `predict_heteroscedastic_gp()` and `forecast_gp()` use the same
  path. `inst/benchmarks/prediction-memory.R` measures time and memory.

* Internal: each leaf kernel is now defined once, in a kernel registry that
  supplies its parameters, constraints, covariance, diagonal, and
  derivatives. Constructors, `evaluate_kernel()`, `kernel_diagonal()`,
  `kernel_gradient()`, parameter paths and updates, printing, sparse FITC
  models, and `optimize_gp()` all read the registry instead of listing kernel
  types separately. Results are unchanged: golden tests generated before the
  change reproduce every covariance, diagonal, gradient, parameter path, and
  printout exactly, and optimizer results and evaluation counts are
  identical. `optimize_gp()` now records each parameter's optimizer
  coordinate in `$optimization$coordinates`; for the built-in kernels every
  coordinate is the log. CONTRIBUTING.md describes how to add a kernel.

* New `sample_gp_prior()` and `sample_gp_posterior()` draw reproducible
  functions from a GP prior and from the latent posterior, and noisy
  observations from the posterior predictive distribution. Draws use the
  Cholesky factor of the covariance with the package's jitter policy, never
  an inverse. Models with known observation-specific noise take
  prediction-time noise variances explicitly, and heteroscedastic models are
  sampled conditionally on their plug-in noise estimate. The regression
  vignette now shows function draws.
* `optimize_gp()` diagnostics now report `function_evaluations`,
  `gradient_evaluations`, and `model_fits`, per start and summed over starts,
  and the printout shows them. Unlike `optim()`'s own `counts`, the function
  evaluations include the objective calls made for finite-difference
  gradients. `inst/benchmarks/optimizer-gradients.R` compares analytical and
  numerical gradients with these counts.
* Every exported function's help page now has a runnable example.
  `R CMD check` and the documentation site run them all.
* `optimize_gp()` now uses the exact gradient of the log marginal likelihood
  instead of finite differences, which needs one GP fit per optimizer step
  rather than two per parameter. `gradient = "numerical"` restores the
  previous behaviour, and the optimization diagnostics and printout record
  which gradient was used.
* `optimize_gp()` now allows 500 L-BFGS-B iterations per start by default
  instead of 100, so composite kernels with periodic components, which can
  need a few hundred iterations, converge.
* New `log_marginal_likelihood_gradient()` returns that gradient with respect
  to log kernel parameters and, for homoscedastic models, the log noise
  variance.
* New `kernel_gradient()` returns analytical derivatives of a covariance
  matrix with respect to the log of every kernel parameter, named like
  `kernel_parameters(kernel, flatten = TRUE)`. It supports all kernels, ARD
  length scales, and sum, product, and scaled kernels. A new vignette
  derives each formula.
* The documentation site now includes the research notes, a worked
  time-series example, and installation instructions.

# gaussianprocesses 0.1.0

First release.

* Exact Gaussian-process regression with zero and constant mean functions,
  known observation-specific noise, and posterior prediction that reports
  latent and observation uncertainty separately.
* Covariance kernels: RBF, Matérn 1/2, 3/2 and 5/2, rational quadratic,
  periodic, linear, and white noise, with ARD length scales and sum, product,
  and scale composition.
* Exact log marginal likelihood and multi-start hyperparameter optimization
  in log-parameter space.
* Exact leave-one-out diagnostics, calibration summaries, and
  numerical-conditioning diagnostics.
* Time-series workflows: explicit time indexing, interpolation and
  forecasting, and rolling-origin evaluation.
* An approximate heteroscedastic GP and sparse FITC regression.
* Reproducible simulation scenarios, benchmarks, and mathematical vignettes.
