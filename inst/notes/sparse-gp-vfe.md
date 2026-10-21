# Sparse Gaussian Processes: VFE, FITC, and Inducing Points

## The collapsed variational bound

For training inputs `X`, inducing inputs `Z`, and noise variances
`Sigma = diag(sigma_i^2)`, both sparse approximations use the low-rank
covariance

```text
Q_ff = K_fu K_uu^{-1} K_uf
```

FITC models the responses as Gaussian with covariance
`Q_ff + diag(K_ff - Q_ff) + Sigma`. The collapsed variational free energy
(VFE) of Titsias (2009) instead keeps `Q_ff + Sigma` and subtracts a trace
penalty:

```text
F = log N(y | m(X), Q_ff + Sigma) - 1/2 sum_i (k(x_i, x_i) - [Q_ff]_ii) / sigma_i^2
```

`F` is a lower bound on the exact log marginal likelihood `log p(y)`, with
equality when the inducing points explain the whole prior variance at the
training inputs, for example when `Z = X`. Its gap is the Kullback-Leibler
divergence between the variational and the exact posterior, so maximizing it
over the hyperparameters can only trade fit for coverage, never overfit the
way an approximate likelihood can. VFE predicts with the deterministic
training conditional, which has the form of the FITC prediction with `Sigma`
in place of FITC's diagonal.

`fit_sparse_gp(method = "vfe")` computes `F` with the factorization of the
FITC implementation: `V = R_uu^{-T} K_uf` with `R_uu' R_uu = K_uu` and the
Cholesky factor of `B = I + V Sigma^{-1} V'`, whose eigenvalues are at least
one. Both approximations report the trace term

```text
t = sum_i (k(x_i, x_i) - [Q_ff]_ii) / sigma_i^2
```

the prior variance that the inducing points leave unexplained, in units of
noise variance. VFE's bound is loose by at least `t / 2`. A trace term much
larger than one says that the inducing points do not cover the data.

## Gradients

`optimize_sparse_gp()` maximizes the FITC objective or the VFE bound with
the optimizer of `optimize_gp()`, with the inducing points fixed. With
`C = Q_ff + Lambda`, `alpha = C^{-1} r`, `M = alpha alpha' - C^{-1}`, and
`P = K_uu^{-1} K_uf`, the Gaussian term changes by `tr(M (dQ + dLambda)) / 2`
and

```text
dQ       = dK_fu P + P' dK_uf - P' dK_uu P
tr(M dQ) = 2 sum(P M * dK_uf) - sum(P M P' * dK_uu)
```

In whitened form `V C^{-1} = B^{-1} V Lambda^{-1}`, so `P M` is an `m x n`
matrix and `P M P'` an `m x m` matrix computed without any `n x n` matrix.
FITC adds `tr(M dLambda) / 2` through its diagonal and VFE the derivative of
its trace term. Both enter only as weighted sums over observations, so the
`m x m` matrix `P W P'`, with the weights in `W`, is accumulated once for all
parameters, and the diagonal terms cost `O(n m^2)` in total rather than
`O(n m^2)` for each of the `p` parameters. The observations are processed in
blocks. The gradients agree
with finite differences to a relative 1e-6 for ARD, composite, and
heteroscedastic models.

## Hyperparameter learning: VFE against FITC

Bauer, van der Wilk, and Rasmussen (2016) report that FITC's objective is
not a bound and can exceed the exact marginal likelihood, that FITC tends to
underestimate the noise variance by using its diagonal correction as
input-dependent noise, and that VFE overestimates the noise when the
inducing points cover the data poorly but recovers the exact
hyperparameters as inducing points are added. `inst/benchmarks/sparse-vfe-fitc.R`
checks these claims with fixed inducing points chosen by farthest-point
selection.

Each of 20 replicates draws 400 training inputs uniformly on `[-4, 4]`,
responses `sin(2x)` plus noise of variance 0.1, and 500 test points, and
estimates an RBF kernel and the noise variance by the exact likelihood, FITC,
and VFE with the same three starts. The exact optimum was the best exact
likelihood reached in every replicate. Means over replicates, with Monte
Carlo standard errors in parentheses; "relative to exact" columns are paired
with the exact fit of the same data.

| m | method | noise relative to exact | objective minus exact log likelihood at the estimate | replicates with objective above exact | trace term | test log density relative to exact |
|---|---|---|---|---|---|---|
| 6 | FITC | 1.56 (0.13) | -72.7 (6.2) | 0 of 20 | 119 (13) | -0.28 (0.035) |
| 6 | VFE | 1.96 (0.12) | -66.2 (3.0) | 0 of 20 | 38 (2.4) | -0.28 (0.034) |
| 10 | FITC | 0.956 (0.009) | +0.47 (0.36) | 9 of 20 (largest +5.05) | 24 (4.2) | -0.00003 (0.0007) |
| 10 | VFE | 1.013 (0.002) | -1.72 (0.09) | 0 of 20 | 3.0 (0.08) | +0.0007 (0.0007) |
| 20 | FITC | 1.000 (< 0.0001) | -0.0003 (0.0002) | 10 of 20, all below 0.0003 | 0.006 | 0 |
| 20 | VFE | 1.000 (< 0.0001) | -0.003 (0.002) | 0 of 20 | 0.005 | 0 |
| 40 | both | 1.000 | about 0 | FITC 8 of 20, below 1e-7 | < 1e-6 | 0 |

What the run shows:

- VFE's bound stayed below the exact log marginal likelihood in all 80
  sparse fits, as the theory requires. FITC's objective exceeded it in 27 of
  80; materially only at `m = 10`, by up to 5.05 and on average 0.47, and by
  less than 3e-4 at `m >= 20`, where FITC is nearly exact.
- With moderate coverage (`m = 10`), FITC underestimated the noise variance
  by 4.4% relative to the exact estimate, about five standard errors, and
  chose hyperparameters under which its trace term was 24, eight times
  VFE's: it lets the inducing points cover the data less well and uses the
  diagonal correction for part of the variance, the mechanism Bauer et al.
  describe. VFE overestimated the noise by 1.3%.
- With poor coverage (`m = 6`) both approximations inflated the noise
  variance, VFE more (by 96%) than FITC (by 56%), and predicted much worse
  than the exact model. From `m = 20` both recovered the exact estimates.
- Neither the FITC bias at `m = 10` nor VFE's predicted differently from the
  exact model on test data within the Monte Carlo error.

With fixed inducing points the run does not show the severe noise
underestimation that Bauer et al. report for FITC; that needs free inducing
points, and the section on optimizing them below reproduces it. Conclusions
are limited to this signal, kernel, and selection rule.

## Choosing inducing points: greedy variance selection

`select_inducing_points(method = "variance", kernel = ...)`, and
`selection = "variance"` in the sparse fitters, add inducing points one at a
time, each at the input whose prior variance is least explained by those
chosen so far: the largest residual `k(x, x) - Q(x, x)`. Burt, Rasmussen, and
van der Wilk (2019, 2020) give approximation guarantees for this choice. It
is the pivot order of a pivoted Cholesky decomposition of `K(X, X)`, and the
tests check it against `chol(pivot = TRUE)`. It is computed incrementally,
evaluating only the `m` chosen columns of `K(X, X)`:

```text
l_t = (K(X, x_i) - L_{t-1} L_{t-1}[i, ]') / sqrt(d_i),   d <- d - l_t^2
```

in `O(n m^2)` time and `O(n m)` memory. Because it greedily minimizes
`tr(K_ff - Q_ff)`, it suits VFE, whose bound loses half that trace over the
noise variance. If the variance is explained before `m` points are chosen,
it stops and warns.

`inst/benchmarks/sparse-inducing.R` times it on two-dimensional inputs with
an RBF kernel of length scale 0.05 (median of three runs, seconds):

| n | m = 50 | m = 100 | m = 200 | m = 400 | m = 800 |
|---|---|---|---|---|---|
| 5 000 | 0.05 | 0.11 | 0.37 | 1.50 | 8.32 |
| 10 000 | 0.19 | 0.47 | 1.97 | 6.42 | 8.91 |
| 20 000 | 0.14 | 0.44 | 2.36 | 6.32 | 23.3 |
| 40 000 | 0.33 | 0.78 | 2.33 | 7.98 | 51.6 |

The log-log slope of time is 0.93 in `n` at `m = 800` and 1.79 in `m` at
`n = 40 000`, consistent with `O(n m^2)` once the `m^2` term dominates; at
small `m` the `O(n m)` kernel evaluations weigh more. Individual timings on
this desktop varied by up to a factor of three, as the `n = 10 000`,
`m = 800` entry shows.

## Optimizing the inducing inputs

`optimize_sparse_gp(optimize_inducing = TRUE)` optimizes the `m x D`
inducing inputs `Z` jointly with the hyperparameters. The gradient reuses
the sensitivities of the objective to `K_uf` and `K_uu` from the
hyperparameter gradient,

```text
G_uf = P M - 2 P W,   G_uu = -P M P' / 2 + P W P'
```

and, since `Z_jd` changes only row `j` of `K_uf` and row and column `j` of
`K_uu`,

```text
dF / dZ_jd = sum_i G_uf[j, i] dk(z_j, x_i) / dz_jd + 2 sum_l G_uu[j, l] dk(z_j, z_l) / dz_jd
```

with the kernel input derivatives of `kernel_input_gradient()`, at the same
`O(n m^2)` cost. The gradients agree with finite differences to a relative
1e-6 for `m` from 3 to 10 and `D` from 1 to 3. Inducing points that become
nearly identical under the kernel, with correlation above `1 - 1e-6`, are
recorded in `$inducing_collapse` with a warning, and `K_uu` is factorized
with the package's jitter policy.

### Gap to the exact marginal likelihood

For each simulation scenario of `gp_simulation_scenario()` with 400
observations, `inst/benchmarks/sparse-inducing.R` maximizes the exact log
marginal likelihood and the VFE bound with inducing points fixed at
farthest-point selection, fixed at greedy variance selection, or optimized
from the greedy selection. The gap is the exact maximum minus the bound's
maximum; the fractions are one minus the optimized gap over the fixed gap.

| scenario | m | gap, farthest | gap, greedy variance | gap, optimized | closed from farthest | closed from greedy |
|---|---|---|---|---|---|---|
| smooth_1d | 5 | 121 | 120 | 72.9 | 40% | 39% |
| smooth_1d | 10 | 4.04 | 4.03 | 1.40 | 65% | 65% |
| smooth_1d | 20 | 4.5e-4 | 1.5e-4 | 1.4e-5 | 97% | 91% |
| near_singular | 5 | 1510 | 1530 | 1330 | 12% | 13% |
| near_singular | 10 | 4.20 | 3.91 | 0.32 | 92% | 92% |
| near_singular | 20 | 0.39 | 0.036 | 2.2e-6 | 100% | 100% |
| ard_2d | 5 | 353 | 352 | 297 | 16% | 16% |
| ard_2d | 10 | 190 | 253 | 105 | 45% | 59% |
| ard_2d | 20 | 59.5 | 47.4 | 31.5 | 47% | 34% |
| heteroscedastic_1d | 5 | 42.3 | 45.0 | 26.1 | 38% | 42% |
| heteroscedastic_1d | 10 | 7.05 | 4.29 | 2.68 | 62% | 38% |
| heteroscedastic_1d | 20 | 0.25 | 0.20 | 0.065 | 74% | 67% |

Optimizing the inducing inputs closed between 12% and 100% of the gap in
these runs: least with `m = 5` on `near_singular`, whose noise variance is
nearly zero so that the trace penalty dominates the bound, and on `ard_2d`,
whose short length scale needs more inducing points; most once the
inducing points are nearly enough. Greedy variance selection gave a smaller
gap than farthest-point selection in nine of the twelve cases and a larger
one in three; on `ard_2d` with `m = 10` it was worse by 63. These numbers
are for one data set per scenario and support no general claim.

### FITC and VFE with free inducing points

On the design of the earlier comparison (20 replicates, `m = 10`, inducing
points initialized by greedy variance selection):

| method | inducing points | noise relative to exact | objective minus exact log likelihood | objective above exact | median trace term | replicates with collapsed points |
|---|---|---|---|---|---|---|
| FITC | fixed | 0.966 (0.007) | +0.17 (0.34) | 9 of 20 | 13.9 | 0 |
| VFE | fixed | 1.013 (0.002) | -1.69 (0.15) | 0 of 20 | 2.8 | 0 |
| FITC | optimized | **0.714 (0.011)** | **+22.9 (1.8)** | **20 of 20** | 163 | **12** |
| VFE | optimized | 1.007 (0.002) | -0.62 (0.04) | 0 of 20 | 1.9 | 0 |

With free inducing points FITC reproduces the failure Bauer et al. describe:
its objective exceeded the exact log marginal likelihood at its own
estimates in every replicate, by 23 on average; it underestimated the noise
variance by 29%; its trace term rose to a median of 163, meaning that the
inducing points moved to where they leave most of the prior variance to the
diagonal correction; and in 12 of 20 replicates two or more inducing points
collapsed onto each other. VFE, by contrast, moved closer to the exact model
when its inducing points were freed: its gap shrank from 1.69 to 0.62 and
its noise estimate from 1.3% to 0.7% above the exact one. Optimizing
inducing points is therefore recommended with VFE and not with FITC.

## Time and memory

`inst/benchmarks/sparse-scaling.R` measures, for two-dimensional inputs, an
ARD RBF kernel, and farthest-point inducing points, the median of three
fresh-process runs on the development machine (a Windows desktop, R 4.5.1,
reference BLAS). "Peak" is R's allocation high-water mark during the step,
the larger of FITC and VFE, which differ by little; the model keeps the
`m x n` cross-covariance.

| n | m | fit FITC (s) | fit VFE (s) | gradient FITC (s) | gradient VFE (s) | predict 1000 VFE (s) | peak fit (MB) | peak gradient (MB) | model (MB) |
|---|---|---|---|---|---|---|---|---|---|
| 5 000 | 100 | 0.31 | 0.15 | 0.58 | 0.34 | 0.01 | 62 | 55 | 4.5 |
| 5 000 | 200 | 0.31 | 0.19 | 0.82 | 0.52 | 0.03 | 66 | 97 | 9.2 |
| 5 000 | 500 | 1.06 | 1.07 | 2.19 | 2.97 | 0.29 | 117 | 182 | 27 |
| 10 000 | 100 | 0.33 | 0.28 | 0.94 | 0.87 | 0.02 | 63 | 103 | 8.7 |
| 10 000 | 200 | 0.70 | 0.78 | 1.66 | 2.52 | 0.03 | 86 | 149 | 17 |
| 10 000 | 500 | 2.16 | 5.35 | 9.64 | 13.06 | 0.11 | 182 | 189 | 47 |
| 20 000 | 100 | 0.64 | 0.53 | 1.61 | 2.11 | 0.01 | 85 | 148 | 17 |
| 20 000 | 200 | 1.03 | 1.31 | 2.45 | 2.36 | 0.04 | 161 | 158 | 33 |
| 20 000 | 500 | 8.77 | 4.89 | 12.36 | 12.08 | 0.29 | 392 | 229 | 86 |
| 50 000 | 100 | 1.91 | 1.78 | 5.60 | 4.94 | 0.03 | 166 | 241 | 42 |
| 50 000 | 200 | 3.41 | 2.46 | 13.44 | 17.92 | 0.06 | 392 | 233 | 81 |
| 50 000 | 500 | 14.34 | 13.69 | 29.47 | 24.43 | 0.17 | 796 | 405 | 202 |

Fitting time grows roughly linearly in `n`, and from `m = 200` to
`m = 500` by three to eight times in these runs, around the six-fold `m^2`
growth of its dominant term: about 14 seconds and 0.8 GB at `n = 50 000`, `m = 500`. FITC and VFE
cost the same, as they share the factorization. A gradient typically costs
two to three fits (between 1.4 and 7 in these runs), so an optimizer step
costs a few fits. Prediction is linear in the number of prediction inputs
and needs about 50 MB at `m = 500`, independent of `n`. Timings varied
between runs by up to a factor of two at the largest sizes, so ratios
between neighbouring rows are rough.

## Sampling

`sample_gp_posterior()` accepts FITC and VFE models. It forms the joint
latent covariance at the `k` sampling inputs,
`K_** - Q_** + K_*u A^{-1} K_u*`, in `O(k m^2 + k^2 m)` time, and factorizes
it with the package's jitter policy in `O(k^3)` time and `k^2` memory, as
for exact models. The number of inducing points enters only through the
`m x k` cross-covariances, so the practical limit is the number of
sampling inputs `k`, a few thousand, not `m`. Draws match the predictive
mean and covariance within five Monte Carlo standard errors in the tests.

## Limitations

- Optimizing the inducing points adds `m D` coordinates to the optimizer,
  and with FITC it overfits, as shown above.
- The VFE bound needs positive noise variances.
- Fitting stores the `m x n` cross-covariance, so memory grows as `n m`.
- The timings come from one machine with the reference BLAS; an optimized
  BLAS would change them, more for the `m^3` and `n m^2` steps than for the
  rest.

## References

- Burt, D. R., Rasmussen, C. E., and van der Wilk, M. (2019). Rates of
  convergence for sparse variational Gaussian process regression.
  *Proceedings of the 36th International Conference on Machine Learning*,
  862-871.
- Burt, D. R., Rasmussen, C. E., and van der Wilk, M. (2020). Convergence of
  sparse variational inference in Gaussian processes regression. *Journal of
  Machine Learning Research*, 21(131), 1-63.
- Titsias, M. (2009). Variational learning of inducing variables in sparse
  Gaussian processes. *Proceedings of the 12th International Conference on
  Artificial Intelligence and Statistics*, 567-574.
- Bauer, M., van der Wilk, M., and Rasmussen, C. E. (2016). Understanding
  probabilistic sparse Gaussian process approximations. *Advances in Neural
  Information Processing Systems*, 29.
- Quiñonero-Candela, J., and Rasmussen, C. E. (2005). A unifying view of
  sparse approximate Gaussian process regression. *Journal of Machine
  Learning Research*, 6, 1939-1959.
