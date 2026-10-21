# Simulation and Benchmark Protocol

## Goals

The benchmark suite separates four different questions:

1. **Statistical recovery** — how close is the posterior mean to the latent
   simulated function?
2. **Approximation accuracy** — how close is FITC to the exact GP posterior?
3. **Numerical stability** — was jitter required and how well-conditioned was
   the stabilized covariance?
4. **Performance** — how much elapsed time and retained model memory were used?

These are deliberately reported separately. A faster model is not labelled
more accurate, and a numerically stable fit is not automatically a better
statistical model.

## Standard scenarios

### smooth_1d

A one-dimensional RBF process with moderate homoscedastic noise.

Purpose:
- baseline correctness;
- interpolation accuracy;
- exact-versus-FITC comparison.

### near_singular

A one-dimensional RBF process with two almost coincident inputs and extremely
small observation noise.

Purpose:
- Cholesky stabilization;
- condition-number diagnostics;
- sensitivity to near-duplicate observations.

### ard_2d

A two-dimensional Matérn-3/2 process with strongly unequal length scales.

Purpose:
- multidimensional input handling;
- ARD covariance geometry;
- inducing-point selection in more than one dimension.

### heteroscedastic_1d

A one-dimensional Matérn-5/2 process with observation noise increasing with
input.

Purpose:
- observation-specific noise handling;
- sparse and exact inference with known heteroscedastic variance.

## Reproducibility

Every public simulation function accepts an explicit integer seed. Seeded
simulation restores the caller's previous global RNG state after drawing, so
running a benchmark does not silently alter later random calculations.

The default suite uses deterministic scenario-specific offsets from the base
seed.

## Timing

Timing uses repeated `system.time()` calls and reports the median elapsed
time. This reduces sensitivity to one unusually slow run but is still not a
microbenchmarking framework.

The benchmark functions are intended for local scientific comparison. They
should not be used as strict CI performance gates because shared-runner load,
BLAS implementations, CPU frequency, and operating system scheduling can
change timings materially.

## Memory

The suite reports retained R object size through `object.size()`.

This is not peak resident memory. It answers a narrower question: how much R
object memory remains associated with the fitted model after fitting.

## External implementations

Comparisons with other implementations (DiceKriging, kernlab, nlme, and
gplite) are part of the validation references in
`inst/validation/references.R`, which the test suite runs and the validation
article tabulates. They fix the hyperparameters where practical, so that
they compare inference rather than optimizers, and compare optimizers
separately by the likelihood they reach.

## Regenerating results

Run:

```r
source("inst/benchmarks/run-suite.R")
```

No proprietary data are required.
