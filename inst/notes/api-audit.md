# API audit

This audit records every exported function of gaussianprocesses, the
conventions they follow, and the decision taken for each (issue #20). The
audit left 56 exported functions (55 before, plus `fit_time_series_gp()`,
the new name of `fit_gp_time()`) and the S3 methods listed at the end.
Functions added since are listed with the issue that added them:
`kernel_input_gradient()` and `predict_gradient_gp()` (#24),
`select_dimensions()` (#25), the five parametric-mean functions (#26),
`changepoint_kernel()` (#28), `spectral_mixture_kernel()` and
`initialize_spectral_mixture()` (#29), five diagnostics functions (#30),
five likelihood functions (#31), three latent-model functions (#32),
`gp_classification_scores()` (#33), `gp_count_scores()` (#34),
`optimize_sparse_gp()` (#35), and `optimize_time_series_gp()` (#37) brought
the total to 84, and removing `fit_gp_time()` (#40) brings it to 83. Since #39 every export also has a stability tier; see
"Stability tiers" below.

## Conventions

- **Function names.**
  - Verbs that act on Gaussian processes end in `_gp`: `fit_gp()`, `predict_gp()`, `optimize_gp()`, `loo_gp()`, `forecast_gp()`.
  - Variants put their qualifier before `_gp`: `fit_sparse_gp()`, `fit_heteroscedastic_gp()`, `fit_time_series_gp()`, `predict_sparse_gp()`.
  - Functions that return a noun start with `gp_`: `gp_diagnostics()`, `gp_time_index()`, `gp_simulation_scenario()`.
- **Two kernel families.**
  - `name_kernel()` builds a reusable kernel specification.
  - `kernel_name()` computes one covariance matrix directly from parameter values.
  - Both families use the same formulas, from the kernel registry.
- **Argument order.**
  - Fitting: `x, y, kernel, noise_variance, mean`, then method options, then numerical options (`initial_jitter`, `fallback_jitter`, `jitter_multiplier`, `max_attempts`, `symmetry_tolerance`).
  - Prediction: `model, x, interval_level, include_covariance, variance_tolerance`, then model-specific options.
  - Issue #20 proposed `x, y, kernel, mean, noise_variance`. That order was not adopted: positional calls such as `fit_gp(x, y, kernel, 0.02)` are common, and reordering would silently change their meaning.
- **Tolerances.** Every tolerance names what it controls:
  - `symmetry_tolerance` for the matrix symmetry check;
  - `convergence_tolerance` for an iterative estimator;
  - `variance_tolerance` for clipping negative posterior variances;
  - `condition_tolerance` for flagging ill-conditioning.

  No exported function takes a plain `tolerance` any more.
- **Prediction results.** Every prediction shares these fields, with the same meaning:
  - `x` and `mean`;
  - `latent_variance` and `latent_sd`;
  - `observation_noise_variance`, `observation_variance`, and `observation_sd`;
  - `interval_level`, `latent_interval`, and `prediction_interval`;
  - with `include_covariance = TRUE`, `latent_covariance` and `observation_covariance`.

  Specialized predictions may add fields; `?predict_gp` documents the shared ones.

  `predict_gradient_gp()` is the exception. It predicts derivatives of the latent function, which have no observation noise and one value per input dimension, so it returns `x`, `interval_level`, and the matrices `mean`, `variance`, `sd`, `lower`, and `upper`, with one row per input and one column per input dimension.

  `predict_latent_gp()` (#32) returns the shared latent fields (`x`, `mean`, `latent_variance`, `latent_sd`, `interval_level`, `latent_interval`, and `latent_covariance`) but replaces the observation fields by `response_mean`, `response_variance`, `probability` and `probability_interval` (Bernoulli likelihoods), and `rate_interval` and `exposure` (Poisson likelihoods). `prediction_interval` (#34) keeps its shared meaning: the interval for a new count, discrete and conservative, or for a new Gaussian observation. A non-Gaussian response has no observation-noise variance, and its prediction intervals are discrete (#34).
- **Deprecation.** A renamed function or argument keeps working for at least one minor release.
  - A renamed function warns through base `.Deprecated()`.
  - A renamed argument warns with a `deprecatedWarning`, and supplying both the old and the new name is an error.
  - NEWS lists every rename, and `?"gaussianprocesses-deprecated"` lists the deprecated functions.

## Decisions

| Decision | Before | After | Rationale |
|---|---|---|---|
| Rename | `fit_gp_time()` | `fit_time_series_gp()` | Qualifier before `_gp`, like the other fitters, and matching `time_series_kernel()`. `fit_gp_time()` was deprecated, and removed in 1.0.0 (#40). |
| Rename argument | `tolerance` (symmetry) in `fit_gp()`, `optimize_gp()`, `fit_sparse_gp()`, `sample_gp_prior()`, both `sample_gp_posterior()` methods, `simulate_gp_data()` | `symmetry_tolerance` | `tolerance` meant different things in different functions. |
| Rename argument | `tolerance` (convergence) in `fit_heteroscedastic_gp()` | `convergence_tolerance` | Same. |
| Rename argument | `numerical_tolerance` in `fit_heteroscedastic_gp()` | `symmetry_tolerance` | Same name as everywhere else. |
| Prediction field | `predict_heteroscedastic_gp()` returned `noise_variance` and no `observation_covariance` | adds `observation_noise_variance` and, with `include_covariance = TRUE`, `observation_covariance`; `noise_variance` kept as a deprecated copy | Shared prediction fields. A list field cannot warn when read, so the copy is documented as deprecated. It was removed in 1.0.0 (#40). |
| Kernel families | `kernel_name()` and `name_kernel()`, undocumented distinction | both kept, documented and cross-linked | 0.1.0 was released with both. Specifications are needed for models; direct matrices are convenient and cost nothing to keep. |
| Base R methods | none | `predict()`, `logLik()`, `nobs()`, `fitted()`, `residuals()`, `summary()` for the exact, sparse, and heteroscedastic classes | Standard R workflows (`predict(model, newdata)`, `AIC()`, `BIC()`). Heteroscedastic models refuse `logLik()`, because the plug-in noise estimate does not maximize a likelihood. |
| Imports | `stats` used through `::`, undeclared | `Imports: stats` | Required for the method registrations and good practice for CRAN. |
| Fitting entry point for non-Gaussian likelihoods (#31) | none | one entry point, `fit_latent_gp(x, y, kernel, likelihood, mean, exposure, ...)`, for every likelihood, with the approximation chosen by an argument (Laplace first, #32); `fit_gp()` remains the exact Gaussian model and gets no `likelihood` argument; no per-likelihood wrappers such as `fit_classification_gp()` before the API freeze (#39) | The likelihood argument already names the model, so a wrapper would save one argument while duplicating every other one and adding an export per likelihood. A `likelihood` argument to `fit_gp()` would give it two result classes and options that apply to only one of them. Likelihoods are specifications like kernels and means, so the entry point follows the fitting argument order with `likelihood` after `kernel`. |
| Latent-model functions (#32) | none | `fit_latent_gp()`, `optimize_latent_gp()`, and `predict_latent_gp()`; `log_marginal_likelihood()`, `log_marginal_likelihood_gradient()`, `sample_gp_posterior()`, and the base methods accept latent models | The entry point decided in #31, with optimization and prediction named like their exact counterparts. Functions that compute something every model has accept the new class instead of gaining latent twins. `optimize_latent_gp()` shares the optimizer of `optimize_gp()` (coordinates, bounds, starts, evaluation counts, and the `$optimization` record), with `optimize_likelihood` in place of `optimize_noise`. |
| Latent-model means (#32, #34) | none | fixed coefficients, coefficients with a Gaussian `coefficient_prior()`, or coefficients estimated by maximum likelihood (#34) | A Gaussian prior folds into the prior covariance exactly. Estimated coefficients maximize the approximate marginal likelihood: `fit_latent_gp()` optimizes them at fixed hyperparameters and `optimize_latent_gp()` jointly as `mean.<name>`, so that both give the same estimate at the optimum. REML has no Laplace counterpart and vague priors are improper, so both are refused. |
| Sparse approximations (#35) | FITC only, fixed hyperparameters | `fit_sparse_gp(method = c("fitc", "vfe"))`; new `optimize_sparse_gp()` with `method = c("vfe", "fitc")`; `predict_sparse_gp()` gains `block_size`; `sample_gp_posterior()` accepts sparse models; sparse models record `method` and `trace_term` (schema 4) | `method` comes after `selection`, so no positional call changes meaning, and FITC stays the default of `fit_sparse_gp()` so that existing fits are unchanged. `optimize_sparse_gp()` is new and defaults to VFE, whose objective is a lower bound and the safer one for learning hyperparameters; the defaults differ deliberately, and both are documented. The optimizer is that of `optimize_gp()`, including noise patterns; inducing points stay fixed until #36. |
| Inducing points (#36) | farthest-point and quantile selection, fixed `Z` | `select_inducing_points(method = "variance", kernel)`, greedy conditional-variance selection; `selection = "variance"` in the sparse fitters; `optimize_sparse_gp(optimize_inducing = FALSE)` optimizes `Z` as `inducing.<j>.<d>`; sparse models record `inducing_collapse` | No new exports. `kernel` comes last in `select_inducing_points()` and is needed only by the variance method. `optimize_inducing` defaults to `FALSE` so that `optimize_sparse_gp()` keeps the meaning it had in #35, and its documentation recommends it for VFE and warns against it for FITC. |
| State-space inference (#37) | exact time-series fits only | `fit_time_series_gp(method = c("exact", "state_space"))`; new `optimize_time_series_gp()` with the same `method`; class `gaussianprocesses_state_space_model` (schema 4), accepted by `forecast_gp()`, `predict_gp()` (at time-index inputs), `log_marginal_likelihood()`, `sample_gp_posterior()`, `rolling_origin_gp()` (through `...`), and the base methods | Kalman smoothing gives the exact posterior, so it is a `method` of the existing fitter, as the issue proposed for `fit_gp_time()`, rather than a new fitter; `method` comes after `unit`, before `...`, so no named call changes. The state-space method needs a time index, so it is offered only through the time-series functions, which also guarantee one sorted input; the engine accepts tied and unsorted inputs, which the tests use. `optimize_time_series_gp()` is new because `optimize_gp()` has no time metadata and no `method`; for `method = "exact"` it calls `optimize_gp()` and records the time index. Arguments that apply only to exact inference (jitter options, `derivative`, `gradient`) are refused for `"state_space"` with a validation error rather than ignored. Unsupported kernels raise `gaussianprocesses_state_space_error`. |
| Counts (#34) | none | `gp_count_scores(observed, latent_mean, latent_variance, exposure, interval_levels, seed)`; `predict_latent_gp()` adds `rate_interval` and `prediction_interval` for Poisson likelihoods, and the samples `rate` | The count predictive is defined by the latent mean and variance and the exposure, so the scores take those rather than a mean and standard deviation. `prediction_interval` has its meaning from `predict_gp()`, an interval for a new observation; it is also filled for Gaussian likelihoods. `gp_holdout_scores()` gains `seed` for the randomized PIT. |
| Classification (#33) | none | no `fit_classification_gp()` wrapper; `gp_classification_scores()` for binary probability forecasts, which `gp_holdout_scores()` returns for latent Bernoulli models; `predict_latent_gp()` adds `probability_interval` and the samples `probability` | Binary classification is `fit_latent_gp()` with `bernoulli_likelihood()`, as decided in #31. Classification scores differ from Gaussian ones (log loss, Brier score, ROC AUC, reliability), so they get their own function with the `gp_*_scores` name; `gp_holdout_scores()` dispatches to it so that held-out scoring has one entry point. |
| Binary responses (#31) | none | 0/1 numbers, logical values (`TRUE` is 1), or a two-level factor whose second level is 1 | The factor coding is that of `glm()`, so a factor response means the same thing in both. Other values, missing values, and other factors are errors. |
| Unchanged | `loo_gp()` against the `gp_*_diagnostics()` names | kept | `loo_gp()` is a verb (compute leave-one-out); the diagnostics functions return named summaries. |
| Unchanged | `simulate_gp_data(x, kernel, mean, noise_variance)` | kept | It has no `y`. Its leading arguments match `sample_gp_prior(x, kernel, mean)`, and reordering would break positional calls. |
| Unchanged | heteroscedastic model field `tolerance` | kept | The model field is internal (`?gaussianprocesses_model`), and renaming it would change the model schema for no user-visible gain. |

## Exported functions

The tables give each function's family, role, arguments in order (`...` where
the function forwards them), what it returns, its main result fields, and
the decision.

### Kernel specifications

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `rbf_kernel` | squared-exponential specification | `variance, length_scale` | `gaussianprocesses_kernel` | `type, parameters, children` | keep | Registry-backed constructor. |
| `matern12_kernel` | Matérn 1/2 specification | `variance, length_scale` | `gaussianprocesses_kernel` | as above | keep | |
| `matern32_kernel` | Matérn 3/2 specification | `variance, length_scale` | `gaussianprocesses_kernel` | as above | keep | |
| `matern52_kernel` | Matérn 5/2 specification | `variance, length_scale` | `gaussianprocesses_kernel` | as above | keep | |
| `rational_quadratic_kernel` | rational-quadratic specification | `variance, length_scale, alpha` | `gaussianprocesses_kernel` | as above | keep | |
| `periodic_kernel` | periodic specification | `variance, length_scale, period` | `gaussianprocesses_kernel` | as above | keep | |
| `linear_kernel` | linear specification | `variance` | `gaussianprocesses_kernel` | as above | keep | |
| `white_noise_kernel` | white-noise specification | `variance` | `gaussianprocesses_kernel` | as above | keep | |
| `sum_kernel` | sum of kernels | `...` | `gaussianprocesses_kernel` | `children` | keep | |
| `product_kernel` | product of kernels | `...` | `gaussianprocesses_kernel` | `children` | keep | |
| `scale_kernel` | scaled kernel | `kernel, scale` | `gaussianprocesses_kernel` | `parameters$scale` | keep | |
| `select_dimensions` | kernel on selected input columns | `kernel, columns` | `gaussianprocesses_kernel` | `columns`, `children` | new (#25) | Named for what it does rather than `*_kernel()`: it wraps a kernel, like `scale_kernel()`, but has no hyperparameter and adds no level to parameter paths. |
| `time_series_kernel` | trend + periodic + local composite | `period, trend_variance, periodic_variance, periodic_length_scale, local_variance, local_length_scale, changepoint, changepoint_steepness` | `gaussianprocesses_kernel` | composite | keep; extended (#28) | `changepoint` wraps the kernel in a `changepoint_kernel()` with independent copies before and after the break. |
| `spectral_mixture_kernel` | spectral mixture of Q components | `n_components, weights, frequencies, spectral_variances` | `gaussianprocesses_kernel` | a `sum_kernel()` of spectral components | new (#29) | Each component is a registry kernel, so gradients, input derivatives, and derivative observations follow the registry contract. Frequencies are real-valued (identity coordinate), which allows the RBF case at 0. |
| `initialize_spectral_mixture` | spectral-mixture kernel from the data's spectrum | `x, y, n_components` | `gaussianprocesses_kernel` | | new (#29) | A verb because it builds its result from data, unlike the `*_kernel()` constructors. |
| `changepoint_kernel` | change between two kernels along one column | `before, after, location, steepness, column` | `gaussianprocesses_kernel` | `parameters$location, parameters$steepness, column, children` | new (#28) | The first composite with a real-valued parameter (`location`, identity coordinate). Child paths are `before.` and `after.`, which read better than `kernel1.` and `kernel2.` for two roles that are not interchangeable. |

### Kernel operations

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `evaluate_kernel` | covariance of a specification | `kernel, x, y` | numeric matrix | | keep | Now explains the two kernel families. |
| `kernel_diagonal` | diagonal without the matrix | `kernel, x` | numeric vector | | keep | |
| `kernel_gradient` | derivatives in log parameters | `kernel, x, y` | named list of matrices | one per parameter path | keep | |
| `kernel_input_gradient` | derivatives in the inputs | `kernel, x, y` | list of matrices | one per input dimension | new (#24) | Same arguments as `kernel_gradient()`. Differentiable kernels only. |
| `kernel_parameters` | nested or flattened parameters | `kernel, flatten` | list or named numeric | parameter paths | keep | |
| `update_kernel_parameters` | update by parameter path | `kernel, parameters` | `gaussianprocesses_kernel` | | keep | |

### Covariance-matrix functions

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `kernel_rbf` | RBF covariance matrix | `x, y, variance, length_scale` | numeric matrix | | keep | Cross-linked to `rbf_kernel()`. |
| `kernel_matern12` | Matérn 1/2 matrix | `x, y, variance, length_scale` | numeric matrix | | keep | |
| `kernel_matern32` | Matérn 3/2 matrix | `x, y, variance, length_scale` | numeric matrix | | keep | |
| `kernel_matern52` | Matérn 5/2 matrix | `x, y, variance, length_scale` | numeric matrix | | keep | |
| `kernel_rational_quadratic` | rational-quadratic matrix | `x, y, variance, length_scale, alpha` | numeric matrix | | keep | |
| `kernel_periodic` | periodic matrix | `x, y, variance, length_scale, period` | numeric matrix | | keep | |
| `kernel_linear` | linear matrix | `x, y, variance` | numeric matrix | | keep | |
| `kernel_white_noise` | white-noise matrix | `x, y, variance` | numeric matrix | | keep | |

### Mean functions

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `zero_mean` | zero mean | none | `gaussianprocesses_mean` | `type, parameters` | keep | |
| `constant_mean` | constant mean | `value` | `gaussianprocesses_mean` | `type, parameters`, or `coefficients` | keep; extended (#26) | `value` also accepts a coefficient treatment, for an estimated or marginalized constant. A number keeps the earlier representation, so saved models read unchanged. |
| `linear_mean` | trend on (1, x1, ..., xD) | `coefficients` | `gaussianprocesses_mean` | `type, coefficients` | new (#26) | Raw inputs, so coefficients are an intercept and slopes. |
| `polynomial_mean` | polynomial trend without interactions | `degree, coefficients, center, scale` | `gaussianprocesses_mean` | `type, degree, center, scale, coefficients` | new (#26) | Inputs are standardized to [-1, 1] by default for conditioning; fitted models store the standardization. |
| `basis_mean` | trend on a user basis | `basis, coefficients` | `gaussianprocesses_mean` | `type, basis, coefficients` | new (#26) | |
| `estimate_coefficients` | GLS coefficients with profile ML or REML | `method` | `gaussianprocesses_coefficients` | `treatment, method` | new (#26) | Coefficients are profiled out in closed form rather than optimized numerically, which is exact and adds no optimizer dimensions. |
| `coefficient_prior` | Gaussian or vague coefficient prior | `mean, covariance` | `gaussianprocesses_coefficients` | `treatment, mean, covariance` | new (#26) | `covariance = NULL` is the vague limit. The spelling *marginalized* follows the package's American English (`optimize_gp()`). |
| `evaluate_mean` | mean values at inputs | `mean, x` | numeric vector | | keep | |

### Likelihoods

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `gaussian_likelihood` | Gaussian observation model | `variance` | `gaussianprocesses_likelihood` | `family, link, parameters` | new (#31) | Data only, like kernel specifications, so models that store it serialize safely; the functions live in an internal table. |
| `bernoulli_likelihood` | binary observation model | `link` (`"logit"`, `"probit"`) | `gaussianprocesses_likelihood` | as above | new (#31) | Named after the distribution, with the link as an argument as in `stats::binomial()`. |
| `poisson_likelihood` | count observation model | `link` (`"log"`) | `gaussianprocesses_likelihood` | as above | new (#31) | Exposure is data, so it is an argument where the likelihood is used, not a parameter. |
| `evaluate_likelihood` | log density and derivatives in the latent values | `likelihood, y, f, exposure` | data frame | `y, f, log_density, first, second, third` | new (#31) | Exposes the quantities that Laplace approximations use, so users can check a likelihood. `y` is returned in the internal coding. |
| `likelihood_predictive` | predictive distribution for a Gaussian latent value | `likelihood, mean, variance, observed, exposure, quadrature_order` | list | `mean, variance, probability, log_density` | new (#31) | Closed forms where they exist, adaptive Gauss-Hermite quadrature otherwise. |

### Latent Gaussian-process models

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `fit_latent_gp` | Laplace approximation at given hyperparameters | `x, y, kernel, likelihood, mean, exposure, method, max_iterations, convergence_tolerance, initial_jitter, fallback_jitter, jitter_multiplier, max_attempts, symmetry_tolerance` | `gaussianprocesses_latent_model` | `mode, alpha, curvature, b_cholesky, jitter, log_marginal_likelihood, laplace` | new (#32) | Fitting argument order with `likelihood` after `kernel`. `method` reserves room for other approximations; only `"laplace"` exists. |
| `optimize_latent_gp` | approximate type-II maximum likelihood | `x, y, kernel, likelihood, mean, exposure, optimize_likelihood, n_starts, start_spread, lower, upper, control, method, max_iterations, convergence_tolerance`, numerical options, `gradient` | `gaussianprocesses_latent_model` | `optimization` | new (#32) | Arguments of `optimize_gp()` where they apply. Likelihood parameters are named `likelihood.<name>`. |
| `predict_latent_gp` | latent posterior and predictive distribution | `model, x, interval_level, include_covariance, variance_tolerance, exposure, quadrature_order, block_size` | `gaussianprocesses_latent_prediction` | shared latent fields, `response_mean, response_variance, probability, probability_interval, rate_interval, prediction_interval, exposure` | new (#32); extended (#33, #34) | See Conventions. `probability_interval` and `rate_interval` are the link applied to the latent interval; `prediction_interval` is the interval for a new count or, with a Gaussian likelihood, a new observation. |

### Exact Gaussian-process regression

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `fit_gp` | exact conditioning | `x, y, kernel, noise_variance, mean, initial_jitter, fallback_jitter, jitter_multiplier, max_attempts, symmetry_tolerance, derivative` | `gaussianprocesses_model` | see `?gaussianprocesses_model` | keep; rename argument; extended (#27) | `tolerance` becomes `symmetry_tolerance`. `derivative` (#27) gives the observation type of each row; it comes after the numerical options so that no positional call changes meaning. |
| `optimize_gp` | type-II maximum likelihood | `x, y, kernel, noise_variance, mean, optimize_noise, n_starts, start_spread, lower, upper, control, initial_jitter, fallback_jitter, jitter_multiplier, max_attempts, symmetry_tolerance, gradient, derivative` | `gaussianprocesses_model` | `optimization` | keep; rename argument; extended (#27) | Same. |
| `predict_gp` | posterior prediction | `model, x, interval_level, include_covariance, variance_tolerance, observation_noise_variance, block_size` | `gaussianprocesses_prediction` | the shared prediction fields | keep | Documents the shared fields. |
| `predict_gradient_gp` | posterior gradient of the latent function | `model, x, interval_level, variance_tolerance, block_size` | `gaussianprocesses_gradient_prediction` | `mean, variance, sd, lower, upper` | new (#24) | Prediction argument order; no `include_covariance`, since only marginal variances are returned. |
| `log_marginal_likelihood` | type-II log likelihood | `model` | numeric | | keep | Also available as `logLik()`. |
| `log_marginal_likelihood_gradient` | its gradient in log parameters | `model` | named numeric | parameter paths, `noise_variance` | keep | |

### Sampling

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `sample_gp_prior` | prior function draws | `x, kernel, mean, n_draws, seed, initial_jitter, fallback_jitter, jitter_multiplier, max_attempts, symmetry_tolerance` | `gaussianprocesses_samples` | `latent, mean, latent_covariance, jitter, relative_jitter, covariance_scale` | keep; rename argument | |
| `sample_gp_posterior` | posterior draws (generic) | `model, x, n_draws, ...`; methods add `type, observation_noise_variance, seed`, numerical options, `symmetry_tolerance` | `gaussianprocesses_samples` | adds `observation` | keep; rename argument | |

### Diagnostics

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `loo_gp` | exact leave-one-out | `model, interval_level` | data frame | `mean, variance, residual, standardized_residual, pit, covered, log_predictive_density` | keep | |
| `gp_calibration_diagnostics` | calibration summaries | `model, interval_level` | list | `empirical_coverage, coverage_error, pit_mean, mean_log_predictive_density`, … | keep | |
| `gp_numerical_diagnostics` | numerical conditioning | `model, condition_tolerance` | list | `numerical_jitter, relative_jitter, covariance_scale, cholesky_attempts, reciprocal_condition_number`, `observation_types`, … | keep; extended (#27) | `observation_types` reports conditioning per observation type. |
| `gp_scores` | proper scores of Gaussian predictions | `observed, mean, sd, interval_levels, reference_mean, reference_variance` | `gaussianprocesses_scores` | `summary, calibration, pointwise` | new (#30) | Log score, CRPS, interval scores and coverage per level, PIT, SMSE, MSLL. Takes plain vectors, so it scores predictions from any source. |
| `gp_holdout_scores` | scores of a model on held-out data | `model, x, y, interval_levels, n_bins, seed, ...` | `gaussianprocesses_scores`, `gaussianprocesses_classification_scores` for latent Bernoulli models, or `gaussianprocesses_count_scores` for latent Poisson models | as above | new (#30); extended (#33, #34) | Predicts with `predict()`, so every model class works. `n_bins` (#33) and `seed` (#34) come before `...`, so no call changes meaning. |
| `gp_count_scores` | scores of count forecasts | `observed, latent_mean, latent_variance, exposure, interval_levels, seed, quadrature_order` | `gaussianprocesses_count_scores` | `summary, calibration, pointwise` | new (#34) | Log score, Dawid-Sebastiani score, Pearson dispersion (the overdispersion check), discrete interval coverage, and randomized PIT values. |
| `gp_classification_scores` | scores of binary probability forecasts | `observed, probability, n_bins, log_density` | `gaussianprocesses_classification_scores` | `summary, reliability, pointwise` | new (#33) | Log loss, Brier score, ROC AUC from mid-ranks, reliability bins, and the expected calibration error. Takes plain vectors, like `gp_scores()`. |
| `gp_loo_scores` | scores of exact leave-one-out predictions | `model, interval_levels` | `gaussianprocesses_scores` | as above | new (#30) | Reuses `loo_gp()`. |
| `gp_hyperparameter_uncertainty` | Laplace approximation and identifiability | `model, step, interval_level, condition_threshold, log_scale_threshold` | `gaussianprocesses_hyperparameter_uncertainty` | `parameters, covariance, correlation, information, poorly_identified, weakest_direction` | new (#30) | Observed information from differences of the analytical gradient on the optimizer scale. |
| `gp_profile_likelihood` | profile likelihood of one hyperparameter | `model, parameter, values, n_values, step` | `gaussianprocesses_profile_likelihood` data frame | `value, log_marginal_likelihood, deviance, converged` | new (#30) | |
| `gp_diagnostics` | all three combined | `model, interval_level, condition_tolerance` | `gaussianprocesses_diagnostics` | `loo, calibration, numerical` | keep | |

### Time series

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `gp_time_index` | time to numeric index | `time, unit` | `gaussianprocesses_time_index` | attributes `origin, unit, input_type` | keep | |
| `fit_time_series_gp` | GP on a time index | `time, y, kernel, unit, method, ...` | `gaussianprocesses_model` or `gaussianprocesses_state_space_model` with `time` | `time` | new name; extended (#37) | Replaces `fit_gp_time()`. `method = "state_space"` (#37) uses Kalman smoothing. |
| `optimize_time_series_gp` | type-II maximum likelihood on a time index | `time, y, kernel, unit, method, ...` | as `fit_time_series_gp()` | `optimization` | new (#37) | `optimize_gp()` with time metadata, or the Kalman-filter likelihood with finite-difference gradients. |
| `forecast_gp` | prediction at time values | `model, future_time, interval_level, include_covariance, variance_tolerance, observation_noise_variance` | `gaussianprocesses_time_forecast` | shared prediction fields, `time, time_index, region` | keep | |
| `rolling_origin_gp` | rolling-origin backtest | `time, y, kernel, initial_window, horizon, step, noise_variance, mean, unit, interval_level, ...` | `gaussianprocesses_backtest` (data frame) | `observed, mean, lower, upper, covered, log_predictive_density` | keep | Scalar noise only until #21. |
| `forecast_metrics` | backtest summaries | `backtest` | list | `rmse, mae, empirical_coverage, mean_log_predictive_density` | keep | |

### Heteroscedastic models

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `fit_heteroscedastic_gp` | plug-in heteroscedastic GP | `x, y, kernel, noise_kernel, mean, initial_noise_variance, noise_floor, noise_ceiling, damping, max_iterations, convergence_tolerance, initial_jitter, fallback_jitter, jitter_multiplier, max_attempts, symmetry_tolerance` | `gaussianprocesses_heteroscedastic_model` | `mean_model, noise_model, training_noise_variance, converged, iterations, history` | keep; rename arguments | `tolerance` and `numerical_tolerance` become specific. |
| `predict_heteroscedastic_gp` | prediction with plug-in noise | `model, x, interval_level, include_covariance, variance_tolerance` | `gaussianprocesses_heteroscedastic_prediction` | shared fields, `noise_log_variance, noise_log_variance_uncertainty, noise_variance_interval` | keep; add shared fields | See Decisions. |

### Sparse approximation

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `fit_sparse_gp` | FITC or VFE approximation | `x, y, kernel, noise_variance, mean, n_inducing, inducing_points, selection, method, fitc_diagonal_floor, initial_jitter, fallback_jitter, jitter_multiplier, max_attempts, symmetry_tolerance` | `gaussianprocesses_sparse_model` | `inducing_points, n_inducing, method, approximation, log_marginal_likelihood, trace_term` | keep; rename argument; extended (#35) | `method` (#35) after `selection`. |
| `optimize_sparse_gp` | sparse type-II maximum likelihood | `x, y, kernel, noise_variance, mean, n_inducing, inducing_points, selection, method, optimize_inducing, optimize_noise, n_starts, start_spread, lower, upper, control, fitc_diagonal_floor`, numerical options, `gradient` | `gaussianprocesses_sparse_model` | `optimization`, `inducing_collapse` | new (#35); extended (#36) | Arguments of `fit_sparse_gp()` and `optimize_gp()`. Analytical gradients in $`O(nm^2)`$ without $`n \times n`$ matrices, including the inducing inputs (#36). |
| `predict_sparse_gp` | FITC or VFE prediction | `model, x, interval_level, include_covariance, variance_tolerance, observation_noise_variance, block_size` | `gaussianprocesses_sparse_prediction` | shared prediction fields, `approximation` | keep; extended (#35) | `block_size` as in `predict_gp()`. |
| `select_inducing_points` | inducing-point selection | `x, n_inducing, method, kernel` | numeric matrix | attributes `indices`, `method` | keep; extended (#36) | `method = "variance"` (#36) with a `kernel`. |
| `benchmark_sparse_gp` | exact against FITC | `x, y, kernel, n_inducing, prediction_x, noise_variance, mean, selection` | list | `comparison, exact_prediction, sparse_prediction` | keep | Scalar noise only until #21. |

### Simulation and benchmarks

| Name | Role | Arguments | Returns | Main fields | Decision | Rationale |
|---|---|---|---|---|---|---|
| `simulate_gp_data` | draw data from a GP | `x, kernel, mean, noise_variance, seed, initial_jitter, fallback_jitter, jitter_multiplier, max_attempts, symmetry_tolerance` | `gaussianprocesses_simulation` | `latent, observed, observation_error, simulation_jitter` | keep; rename argument | |
| `gp_simulation_scenario` | named scenarios | `scenario, n, seed` | `gaussianprocesses_simulation` | adds `scenario` | keep | |
| `run_gp_benchmark_suite` | benchmark suite | `scenarios, n, n_inducing, repeats, seed` | `gaussianprocesses_benchmark_suite` | `accuracy, performance, stability, settings` | keep | |
| `benchmark_exact_gp_scaling` | exact-GP timing | `sizes, repeats, seed` | data frame | `n, median_fit_elapsed_seconds, retained_model_bytes`, … | keep | |
| `compare_gp_predictions` | prediction differences | `reference, candidate` | list | `posterior_mean_rmse, latent_variance_mae`, … | keep | |

## Stability tiers (#39)

`R/stability.R` is the source of truth for the tiers, and the tests check that the help pages and the reference index agree with it. The policy is on `?gaussianprocesses` and in CONTRIBUTING.md.

| Tier | Exports |
|---|---|
| Experimental: approximate inference | `fit_sparse_gp`, `optimize_sparse_gp`, `predict_sparse_gp`, `select_inducing_points`, `fit_latent_gp`, `optimize_latent_gp`, `predict_latent_gp`, `gaussian_likelihood`, `bernoulli_likelihood`, `poisson_likelihood`, `evaluate_likelihood`, `likelihood_predictive`, `gp_classification_scores`, `gp_count_scores`, `fit_heteroscedastic_gp`, `predict_heteroscedastic_gp` |
| Experimental: benchmark helpers | `benchmark_exact_gp_scaling`, `benchmark_sparse_gp`, `run_gp_benchmark_suite`, `compare_gp_predictions`, `gp_simulation_scenario` |
| Experimental: recent kernels and methods | `optimize_time_series_gp`, `spectral_mixture_kernel`, `initialize_spectral_mixture`, `changepoint_kernel`, `predict_gradient_gp`, `kernel_input_gradient`; also `fit_gp(derivative = )`, `fit_time_series_gp(method = "state_space")`, and `time_series_kernel(changepoint = )` |
| Experimental: uncertainty and diagnostics | `gp_hyperparameter_uncertainty`, `gp_profile_likelihood`, `gp_numerical_diagnostics`, `gp_calibration_diagnostics` |
| Deprecated | none. The names deprecated after 0.1.0 (`fit_gp_time`, the `tolerance` and `numerical_tolerance` arguments, and the heteroscedastic `noise_variance` field) were removed in 1.0.0 without a release of warnings (#40). |
| Stable | every other export (52) |

A stable function's results and options that concern an experimental model class, kernel, or argument follow the experimental tier: for example, `predict_gp()` for a state-space model.

## S3 methods

| Generic | Classes | Notes |
|---|---|---|
| `predict` | exact, sparse, heteroscedastic, latent, state-space models | Delegates to `predict_gp()`, `forecast_gp()` (time-series models), `predict_sparse_gp()`, `predict_heteroscedastic_gp()`. Without `newdata` it predicts at the training inputs. |
| `logLik`, `nobs`, `fitted`, `residuals`, `summary`, `print` | state-space models | New (#37). As for exact models: the log marginal likelihood from the Kalman filter, and the smoother means as fitted values. |
| `logLik` | exact, sparse models | Type-II log marginal likelihood (the FITC objective or the VFE bound for sparse models). `df` is the number of hyperparameters estimated by `optimize_gp()` or `optimize_sparse_gp()`, otherwise 0. The heteroscedastic method raises an error that explains why there is none. |
| `logLik`, `nobs`, `fitted`, `summary` | latent models | New (#32). `logLik()` is the Laplace approximation, `fitted()` the posterior mode. There is no `residuals()` method, because the responses are not on the scale of the latent function. |
| `coef` | latent models | New (#34). Fixed, estimated, or approximate posterior-mean coefficients; no `vcov()`. |
| `nobs`, `fitted`, `residuals` | exact, sparse, heteroscedastic models | `fitted()` is the posterior latent mean at the training inputs. |
| `summary`, `print.summary.gaussianprocesses_model` | exact, sparse, heteroscedastic models | Compact summary of the model. |
| `print` | models, predictions, samples, diagnostics, simulations, backtests, kernels, benchmark suites | Unchanged. Gradient predictions (#24), mean specifications (#26), likelihoods (#31), latent models and their predictions (#32), classification scores (#33), and count scores (#34) have one too. |
| `coef`, `vcov` | exact models | New (#26). Mean coefficients and their covariance given the hyperparameters; empty for the zero mean, zero covariance for fixed coefficients. |
| `sample_gp_posterior` | exact, heteroscedastic, latent, sparse (#35), and state-space (#37) models, default | Unchanged apart from `symmetry_tolerance`. The latent method (#32) draws observations from the likelihood and takes `exposure`; it returns `probability` (#33) or `rate` (#34) draws. |
