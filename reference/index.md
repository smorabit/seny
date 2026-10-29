# Package index

## Entry point

Run the whole assessment, or build the pseudobulk object the layers
share.

- [`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
  : Assess every pseudobulk comparison end to end
- [`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md) :
  Set up pseudobulk comparisons from a single-cell dataset

## Diagnostic layers

The per-comparison diagnostics, in order.

- [`assessDesign()`](https://smorabit.github.io/seny/reference/assessDesign.md)
  : Layer 1: design validity and confounding
- [`assessPseudobulk()`](https://smorabit.github.io/seny/reference/assessPseudobulk.md)
  : Layer 2: pseudobulk quality
- [`estimateVariance()`](https://smorabit.github.io/seny/reference/estimateVariance.md)
  : Layer 4 (variance): null-model dispersion, decomposition and
  per-gene SE
- [`estimateSensitivity()`](https://smorabit.github.io/seny/reference/estimateSensitivity.md)
  : Layer 4 (sensitivity): minimum detectable effect and power

## Calling and reporting

Assemble the comparison table, make the call, and show the matrix.

- [`callComparisons()`](https://smorabit.github.io/seny/reference/callComparisons.md)
  : Assemble the comparison table and make the adequacy call
- [`defaultThresholds()`](https://smorabit.github.io/seny/reference/defaultThresholds.md)
  : Default decision thresholds for the adequacy call
- [`plotAdequacy()`](https://smorabit.github.io/seny/reference/plotAdequacy.md)
  : Traffic-light matrix of adequacy calls
