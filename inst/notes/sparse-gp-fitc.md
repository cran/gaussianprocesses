# Sparse Gaussian Processes: Initial Approximation

## Choice: FITC

The first sparse approximation implemented in this package is the Fully
Independent Training Conditional (FITC) approximation.

For training inputs `X` and inducing inputs `Z`, define

```text
K_uu = K(Z, Z)
K_uf = K(Z, X)
Q_ff = K_fu K_uu^{-1} K_uf
```

FITC replaces the exact training covariance by

```text
Q_ff + diag(K_ff - Q_ff) + Sigma_noise
```

where `Sigma_noise` is the observation-noise covariance.

The diagonal correction is important: it restores the exact prior marginal
variance at every training input while retaining the low-rank structure.

## Why FITC first?

The package needs a sparse method that is:

- probabilistically interpretable;
- implementable from the existing kernel and Cholesky layers;
- explicit enough to inspect mathematically;
- useful for prediction as well as training;
- substantially cheaper than the exact GP when `m << n`.

FITC has dominant training complexity

```text
O(n m^2 + m^3)
```

and stores an `m x n` cross-covariance matrix. The exact implementation uses
an `n x n` covariance and an `O(n^3)` Cholesky factorization.

## References

- Quiñonero-Candela, J., & Rasmussen, C. E. (2005).
  *A Unifying View of Sparse Approximate Gaussian Process Regression*.
  Journal of Machine Learning Research, 6, 1939-1959.
- Snelson, E., & Ghahramani, Z. (2006).
  *Sparse Gaussian Processes using Pseudo-inputs*.
  Advances in Neural Information Processing Systems 18.

`optimize_sparse_gp(optimize_inducing = TRUE)` optimizes the inducing
locations jointly with the hyperparameters; with FITC that is prone to
overfitting, as the research note on VFE and FITC shows.

## Inducing-point initialization

Two deterministic strategies are provided.

### Farthest-point selection

Inputs are standardized by column, the point nearest the standardized centroid
is selected first, and later points maximize distance from the selected set.
This gives a reproducible space-filling subset without an external clustering
dependency.

### Quantile selection

For one-dimensional inputs, points nearest equally spaced empirical quantiles
are selected.

Explicit inducing points can also be supplied. They do not need to coincide
with training inputs.

## Numerical formulation

With

```text
Lambda = diag(K_ff - Q_ff) + Sigma_noise
A = K_uu + K_uf Lambda^{-1} K_fu
```

the implementation works with the `m x m` matrices `K_uu` and `A`.
All solves use Cholesky factorizations. No explicit matrix inverse is formed.

A very small positive floor can be applied to the diagonal of `Lambda` before
inversion. The number of floor applications is recorded in the fitted model.

## Prediction

The latent predictive mean is

```text
m_* + K_*u A^{-1} K_uf Lambda^{-1} (y - m)
```

and the latent covariance is

```text
K_** - K_*u K_uu^{-1} K_u*
      + K_*u A^{-1} K_u*
```

Observation noise is added separately, preserving the package-wide distinction
between latent uncertainty and measurement uncertainty.

## Known limitations

- FITC can differ materially from the exact GP when the inducing set is too
  small or poorly located.
- Inducing locations stay fixed unless `optimize_sparse_gp()` optimizes them,
  which is safe with VFE but lets FITC underestimate the noise.
- FITC's approximate marginal likelihood should not be compared numerically to
  the exact marginal likelihood as though they were the same model.
- `benchmark_sparse_gp()` reports retained model-object memory rather than
  peak allocation; `inst/benchmarks/sparse-scaling.R` measures peak
  allocation.
- Full posterior covariance at many test points is still quadratic in the
  number of prediction points; the default sparse prediction path returns only
  marginal variances.
- The collapsed variational approximation of Titsias (2009) is available as
  `fit_sparse_gp(method = "vfe")`; the research note on VFE and FITC compares
  them.

Exact GP regression remains the reference implementation throughout the
package.
