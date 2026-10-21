# Heteroscedastic Gaussian Processes

## Scope

This project uses a transparent two-process approximation for input-dependent
Gaussian observation noise.

The observation model is

```text
y_i = f(x_i) + epsilon_i
epsilon_i ~ Normal(0, sigma^2(x_i))
```

with one GP for the latent regression function `f` and a second GP for
`g(x) = log sigma^2(x)`.

This follows the general two-GP idea introduced by Goldberg, Williams and
Bishop (1998) and the practical most-likely-noise direction of Kersting,
Plagemann, Pfaff and Burgard (2007). The implementation here is deliberately
simpler and is not claimed to reproduce either paper exactly.

References:

- Goldberg, P. W., Williams, C. K. I., & Bishop, C. M. (1998).
  *Regression with Input-Dependent Noise: A Gaussian Process Treatment*.
  Advances in Neural Information Processing Systems 10.
- Kersting, K., Plagemann, C., Pfaff, P., & Burgard, W. (2007).
  *Most Likely Heteroscedastic Gaussian Process Regression*.
  ICML 2007, 393-400. DOI: 10.1145/1273496.1273546.

## Approximation used here

1. Fit a mean GP with an initial scalar noise estimate.
2. Compute exact leave-one-out predictive residuals from that fit.
3. Treat corrected log squared residuals as noisy observations of
   `log sigma^2(x)`; the correction is described below.
4. Fit a second GP to those log-noise targets.
5. Predict a noise variance at each training input.
6. Refit the mean GP using those values as a diagonal observation-noise
   covariance.
7. Repeat until the maximum log-scale change in training noise variance is
   small or the iteration limit is reached.

For a Gaussian residual

```text
e ~ Normal(0, sigma^2)
```

we have

```text
log(e^2) = log(sigma^2) + log(chi^2_1)
```

so the transformed residual has known additive-noise moments. The
implementation subtracts

```text
digamma(1/2) + log(2)
```

from the log squared residual and uses

```text
trigamma(1/2)
```

as the observation variance of the transformed target.

## Leave-one-out residuals include latent uncertainty

A leave-one-out residual is not a pure noise draw. With `v_i` the
leave-one-out posterior variance of the latent function at `x_i`, under the
model

```text
e_i = y_i - mu_{-i}(x_i) ~ Normal(0, v_i + sigma^2(x_i))
E[log(e_i^2)] - digamma(1/2) - log(2) = log(sigma^2(x_i)) + log(1 + v_i / sigma^2(x_i))
```

so the log-chi-squared correction alone leaves an upward bias of
`log(1 + v_i / sigma^2(x_i))`. The bias is largest where the latent function
is poorly determined: small samples, sparse inputs, and low noise.

`loo_gp()` returns `v_i + s_i^2`, where `s_i^2` is the current noise
estimate. Three ways of forming the target were compared:

- **residual**, the original rule: `log(e_i^2) - digamma(1/2) - log(2)`.
- **standardized**: the same, minus `log((v_i + s_i^2) / s_i^2)`. This
  amounts to `log(s_i^2) + log(z_i^2) - digamma(1/2) - log(2)` with the
  standardized residual `z_i = e_i / sqrt(v_i + s_i^2)`. At a fixed point
  `s_i^2 = sigma^2(x_i)`, the target has expectation `log(sigma^2(x_i))`.
- **subtracted**: `log(max(e_i^2 - v_i, floor)) - digamma(1/2) - log(2)`.
  It subtracts the latent variance before taking logs.

## Simulation study

`inst/benchmarks/heteroscedastic-noise.R` fits `heteroscedastic_1d` data,
whose noise variance `0.015 + 0.25 x^2` is known. Each fit uses the true mean
kernel and an RBF noise kernel, with 50 replicates for each sample size and
up to 200 iterations. For each fit, the error is
`log(estimated / true noise variance)` at the training inputs. The bias is
the mean over replicates of the median error, and RMSE is the mean
root-mean-square error. Monte Carlo standard errors are in parentheses.

| Rule | n | Bias | RMSE | Converged within 6 / 30 / 50 / 200 iterations |
|---|---|---|---|---|
| residual | 40 | 0.235 (0.056) | 0.629 (0.030) | 0 / 0.92 / 0.92 / 0.92 |
| residual | 100 | 0.139 (0.035) | 0.387 (0.026) | 0 / 1.00 / 1.00 / 1.00 |
| residual | 200 | 0.082 (0.022) | 0.271 (0.014) | 0 / 0.96 / 0.98 / 0.98 |
| standardized | 40 | -0.060 (0.063) | 0.576 (0.032) | 0 / 0.88 / 0.90 / 0.90 |
| standardized | 100 | 0.017 (0.035) | 0.361 (0.022) | 0 / 0.98 / 0.98 / 0.98 |
| standardized | 200 | 0.015 (0.022) | 0.252 (0.012) | 0 / 0.98 / 0.98 / 0.98 |
| subtracted | 40 | -1.372 (0.067) | 1.456 (0.066) | 0 / 0.46 / 0.46 / 0.46 |
| subtracted | 100 | -1.032 (0.043) | 1.076 (0.038) | 0 / 0.60 / 0.62 / 0.62 |
| subtracted | 200 | -0.839 (0.026) | 0.884 (0.025) | 0 / 0.80 / 0.80 / 0.80 |

What the study shows:

- **Residual.** The original rule overestimates the noise at every sample
  size, by about four Monte Carlo standard errors. The bias shrinks with `n`,
  as the derivation predicts.
- **Standardized.** Its bias is within about one Monte Carlo standard error
  of zero at every sample size, and its RMSE is lower than the original
  rule's at every sample size. It is the default.
- **Subtracted.** Truncating `e_i^2 - v_i` biases this rule strongly
  downwards, and it often fails to converge.
- **Convergence.** No fit converged within the former default of six
  iterations. Most converge in 11 to 15. With 50 iterations, every fit that
  converges at all has done so, so the default `max_iterations` is now 50.
  The remaining fits, 10% at `n = 40` and 2% at larger `n`, cycle without
  converging even in 200 iterations, and `fit_heteroscedastic_gp()` now
  warns whenever it stops before converging.

The rule can be changed for such studies with the internal option
`gaussianprocesses.heteroscedastic_target`; it is not part of the documented
API.

## What the model returns

The prediction API keeps three ideas separate:

- **latent variance**: posterior uncertainty about `f(x)`;
- **noise variance**: the plug-in estimate of `sigma^2(x)`;
- **log-noise uncertainty**: posterior uncertainty of the second GP.

The ordinary observation predictive variance is

```text
latent variance + estimated noise variance
```

The uncertainty of the noise GP itself is not added to this value. Doing so
would mix epistemic uncertainty about the noise function with aleatoric
measurement noise.

## Known limitations

- This is not full joint Bayesian inference.
- Hyperparameters of the mean and noise kernels are not jointly integrated.
- Squared residuals are intrinsically noisy estimates of local variance.
- Sparse data can make the noise process weakly identified.
- The noise floor and damping parameter affect numerical robustness.
- Very sharp changes in noise variance can be oversmoothed by the noise kernel.
- The lognormal transformation means posterior summaries of noise variance are
  asymmetric; the implementation uses `exp(E[g(x)])` as a plug-in estimate.
- Mean-model misspecification can leak into the estimated noise process.
- The bias correction is unbiased only at the fixed point of the iteration;
  the study shows no bias beyond Monte Carlo error for the scenario tested,
  not for every design.
- A minority of fits cycle instead of converging. They return the last
  iterate with `converged = FALSE` and a warning.

These limitations are intentional and documented rather than hidden behind a
single black-box fit.
