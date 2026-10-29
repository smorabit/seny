# seny: sanity checks for differential expression in noisy data

For every cluster-by-contrast comparison in a multi-sample single-cell
dataset, seny decides whether a pseudobulk differential expression test
is Adequate, Underpowered, Caution or Invalid, and explains why. The
diagnostics are analytic and run on the user's own data: design validity
and confounding (layer 1), pseudobulk quality (layer 2), and variance
structure plus sensitivity (layer 4), the last giving a per-gene minimum
detectable effect from a negative-binomial variance model. No
per-dataset simulation.

## Details

The name is a backronym – Sanity checks for differential gene Expression
in NoisY data – and also the Catalan word for common sense.

The main entry point is
[`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md),
which chains
[`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md),
[`assessDesign()`](https://smorabit.github.io/seny/reference/assessDesign.md),
[`assessPseudobulk()`](https://smorabit.github.io/seny/reference/assessPseudobulk.md),
[`estimateVariance()`](https://smorabit.github.io/seny/reference/estimateVariance.md),
[`estimateSensitivity()`](https://smorabit.github.io/seny/reference/estimateSensitivity.md)
and
[`callComparisons()`](https://smorabit.github.io/seny/reference/callComparisons.md).

## Author

**Maintainer**: Samuel Morabito <sam.morabito@cnag.eu>

Authors:

- Samuel Morabito <sam.morabito@cnag.eu>
