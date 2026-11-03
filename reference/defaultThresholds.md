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
threshold), `sensitivity` (expected recall at the target effect at or
above which a comparison is adequately powered; 0.5 matches replication
failure defined as recall below 0.5), `minSamples` (fewer samples per
group is Invalid), `cellsCondCor` (absolute correlation of log cells per
sample with condition above which cell-count imbalance is flagged
Caution), `cellsCondAlpha` (the correlation must also exceed the
critical value at this level for the comparison's sample size, so chance
correlations in small designs are not flagged: the cutoff is
`max(cellsCondCor, r_crit(n))`, about 0.81 at 6 samples and
`cellsCondCor` from ~16 samples), and `minTestable` (fewer testable
genes is Underpowered). The layer-2 values, like the others, are
provisional placeholders until E1-E8. The composition values (see
[`assessComposition()`](https://smorabit.github.io/seny/reference/assessComposition.md))
were calibrated on imposed subtype shifts in OneK1K (E6): `compFDR`
(propeller FDR below which a child's share within its parent is called
shifted), `compFrac` (expected DEGs from that shift alone, as a share of
testable genes, at or above which the comparison is flagged Caution;
0.01 maximised Youden's J for predicting artefacts in at least 1% of
genes, with no flags without a shift) and `compExplained` (share of a
parent DEG's observed log2 fold change the predicted shift-only change
must explain for the gene to be attributed to composition; 0.25, since
the observed change of a detected DEG is inflated by selection).
