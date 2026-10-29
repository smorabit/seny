# Default decision thresholds for the adequacy call

Provisional values used by
[`callComparisons()`](https://smorabit.github.io/seny/reference/callComparisons.md).
The final numbers are to be derived from the validation experiments
(E1-E8); until then these are sensible placeholders, documented so a
user can override them.

## Usage

``` r
defaultThresholds()
```

## Value

A named list: `cramerV` (condition-covariate association above which a
comparison is flagged Caution), `vif` (variance inflation factor
threshold), `sensitivity` (sensitivity fraction at or above which a
comparison is adequately powered), `minSamples` (fewer samples per group
is Invalid), `cellsCondCor` (absolute correlation of log cells per
sample with condition above which cell-count imbalance is flagged
Caution), and `minTestable` (fewer testable genes is Underpowered). The
last two gate the layer-2 metrics and, like the others, are provisional
placeholders until E1-E8.
