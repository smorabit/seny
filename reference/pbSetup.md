# Set up pseudobulk comparisons from a single-cell dataset

Aggregates a single-cell dataset to pseudobulk once per annotation level
and returns the object the rest of the package operates on: a named
[S4Vectors::SimpleList](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html)
holding one pseudobulk
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
per annotation set. Each pseudobulk is a gene-by-(sample x cluster)
matrix with per-pseudobulk metadata in `colData`, and the comparison
definitions, column names and per-cluster cell-level overdispersion
(`theta`) in `metadata()`.

## Usage

``` r
pbSetup(
  sce,
  sample,
  group,
  annotations,
  covariates = character(0),
  contrasts = NULL,
  pairing = NULL,
  minCells = 10,
  aggregation = c("sum", "mean"),
  assayName = "counts",
  thetaCells = 50,
  thetaDonors = 50,
  thetaReps = 3,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- sce:

  A
  [SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  or
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with a raw-count assay (non-negative integers).

- sample:

  Character scalar; `colData` column giving the biological sample
  (donor) each cell belongs to.

- group:

  Character scalar; `colData` column giving the condition compared.

- annotations:

  Character vector of `colData` columns, each a clustering (annotation
  level) to pseudobulk, e.g. `c("cell_type_broad")`.

- covariates:

  Character vector of `colData` columns to carry as candidate
  confounders. Default none.

- contrasts:

  Optional list of length-2 character vectors `c(A, B)` naming `group`
  levels to compare (B vs A). Default `NULL`: all non-reference levels
  versus the first level.

- pairing:

  Optional `colData` column for paired designs (reserved; not used in
  the MVP).

- minCells:

  Integer; pseudobulks with fewer contributing cells are dropped.
  Default 10.

- aggregation:

  Either `"sum"` (counts, for the negative-binomial model) or `"mean"`
  (reserved for the limma path). Default `"sum"`.

- assayName:

  Name of the raw-count assay in `sce`. Default `"counts"`.

- thetaCells, thetaDonors, thetaReps:

  Subsample size (cells per donor, donors) and number of draws for the
  cell-level `theta` ensemble.

- BPPARAM:

  A `BiocParallelParam` object for the per-cluster theta estimation.
  Default
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html).

## Value

A named
[S4Vectors::SimpleList](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html)
of
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
one per annotation level.

## Details

Cell-level `theta` (within-donor cell-to-cell overdispersion, needed by
[`estimateVariance()`](https://smorabit.github.io/seny/reference/estimateVariance.md))
is estimated here, while the cells are in hand, so the later steps never
need the single-cell object again.
