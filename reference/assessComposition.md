# Composition: cluster shares, subtype shifts and composition-driven DE

A broad cluster can hold subtypes whose shares shift between the groups.
Pseudobulk DE on the broad cluster then reports the subtypes' marker
genes as differentially expressed although no gene changes within any
subtype. This step measures that risk for every comparison, in four
parts:

## Usage

``` r
assessComposition(
  x,
  variance,
  sensitivity,
  power = 0.8,
  fdr = 0.05,
  reference = c("reference", "pooled"),
  thresholds = defaultThresholds()
)
```

## Arguments

- x:

  A
  [S4Vectors::SimpleList](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html)
  from
  [`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md).

- variance:

  The list returned by
  [`estimateVariance()`](https://smorabit.github.io/seny/reference/estimateVariance.md).

- sensitivity:

  The table returned by
  [`estimateSensitivity()`](https://smorabit.github.io/seny/reference/estimateSensitivity.md)
  (its per-gene `alpha` sets the MDE).

- power:

  Target power for the MDE. Default 0.8.

- fdr:

  FDR at which parent and child genes count as DE for the attribution.
  Default 0.05.

- reference:

  Where the per-cell child expression `c_gk` comes from: `"reference"`
  (the contrast's reference group A only; default, which avoids
  attributing real within-child DE in B to composition) or `"pooled"`
  (all samples of the contrast's groups).

- thresholds:

  A list from
  [`defaultThresholds()`](https://smorabit.github.io/seny/reference/defaultThresholds.md);
  uses `compExplained`.

## Value

A list with `summary` (a `data.frame`, one row per comparison: the
cluster share `share_A`, `share_B`, `share_fdr`; the `composition`
status; the most shifted child `comp_child` with `comp_share_A`,
`comp_share_B` and `comp_min_fdr`; `exp_comp_deg`, `frac_comp_expected`;
`n_comp_sensitive`, `frac_comp_sensitive`; `n_parent_deg`,
`n_comp_driven`, `frac_comp_driven`), `genes` (per gene and comparison:
`lfc_comp`, `comp_power`, `mde`, `sensitive`, the observed `logfc`,
`padj`, `parent_de`, `child_de` and the attribution `label`) and
`proportions` (per sample and child, the child's share of the parent).

## Details

1.  **Cluster share** (every annotation level, reported only): the
    cluster's share of each sample's cells, A vs B, tested with
    [`speckle::propeller`](https://rdrr.io/pkg/speckle/man/propeller.html)
    (arcsin-square-root proportions, limma). A cluster whose share
    shifts is not a problem for DE within it, so this never changes the
    call.

2.  **Subtype shift within the cluster** (needs a finer annotation
    nested in this one, see
    [`pbSetup()`](https://smorabit.github.io/seny/reference/pbSetup.md)):
    propeller on the children's shares of the parent's cells, over the
    comparison's samples.

3.  **Composition-sensitive genes** (prospective): per gene, the log2
    fold change the observed subtype shift alone would produce in the
    parent's library-size-normalised pseudobulk,
    `log2[(sum_k pB_k c_gk / sum_k pB_k L_k) / (sum_k pA_k c_gk / sum_k pA_k L_k)]`,
    with `p` the donor-averaged child shares, `c_gk` the per-cell mean
    count of gene g in child k and `L_k` the child's per-cell library
    size. Each gene's power to be called DE at that shift-only change
    (layer 4, at the comparison's per-gene alpha) is `comp_power`; their
    sum, `exp_comp_deg`, is the number of DEGs the shift alone is
    expected to produce, and `frac_comp_expected` its share of the
    testable genes.
    [`callComparisons()`](https://smorabit.github.io/seny/reference/callComparisons.md)
    flags the comparison **Caution** when a child's shift is significant
    and `frac_comp_expected` reaches `thresholds$compFrac`. Genes whose
    shift-only change exceeds their own MDE (`sensitive`, counted in
    `n_comp_sensitive`) are reported too, but many genes each at modest
    power add up to many artefacts without any of them clearing the MDE,
    so that count understates the risk (E6).

4.  **Attribution of the observed parent DEGs** (retrospective, reported
    only): each gene the method's own test calls DE in the parent at
    `fdr` is labelled *composition-driven* (DE in no child comparison,
    and the shift-only change has the same sign and explains at least
    `thresholds$compExplained` of the observed change), *within-child
    DE* (DE in at least one child), *mixed* (both) or *unresolved*
    (neither, or no child comparison could be fitted).

Clusters without a finer annotation nested in them get
`composition = "not assessed"`; parts 2-4 need two annotation levels.
