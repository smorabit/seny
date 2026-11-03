# Assessing pseudobulk DE adequacy with seny

## What seny does

For every cluster-by-contrast comparison in a multi-sample single-cell
dataset, `seny` decides whether a pseudobulk differential-expression
test is **Adequate**, **Underpowered**, **Caution** or **Invalid**, and
explains why. The diagnostics are analytic and run on your own data;
there is no per-dataset simulation.

The call has a severity order — Invalid \> Caution \> Underpowered \>
Adequate — and a `reason` string that lists *every* flag, so a single
headline never hides a second problem:

| Call | Meaning |
|----|----|
| **Invalid** | Type I error cannot be trusted (design not estimable, too few samples per group). |
| **Caution** | Valid, but a specific risk applies (condition confounded with a covariate, cell-count imbalance tracking the condition, or a subtype shift inside the cluster that would make genes look DE on its own). |
| **Underpowered** | Valid but unable to detect the target effect (sensitivity too low, or too few testable genes). |
| **Adequate** | No flags. |

## The example dataset

`seny` ships a small simulated dataset: 24 donors (12 per condition),
three cell types, each planted in a different regime so the call varies.
A finer annotation, `cell_type_fine`, splits CD4 T into naive and memory
subtypes; we come back to it under *Composition*.

``` r

library(seny)
sce <- readRDS(system.file("extdata", "example_sce.rds", package = "seny"))
table(sce$cell_type, sce$condition)
#> Loading required namespace: SingleCellExperiment
#>       
#>           A    B
#>   CD4T 1080 1080
#>   CD8T  840  180
#>   NK    240  240
```

## Running the assessment

[`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
chains the layers — pseudobulk setup, design/confounding (layer 1),
pseudobulk quality (layer 2), variance structure and sensitivity (layer
4) — and returns one row per comparison.

``` r

res <- assessAll(sce, sample = "donor", group = "condition",
                 annotations = "cell_type", targetLFC = 1)
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
```

The comparison table carries the layer-2 cell-quality metrics alongside
the call. `cells_cond_cor` is the correlation of log cells-per-sample
with the condition; `cells_ratio` is the smaller group median over the
larger; `expected_recall` is the share of true effects of the target
size the test is expected to detect (mean per-gene power at the target),
and drives the Adequate / Underpowered call.

``` r

cols <- c("cluster", "nA", "nB", "median_cells", "gini_cells", "cells_ratio",
          "cells_cond_cor", "n_testable", "expected_recall", "call")
knitr::kable(as.data.frame(res)[, cols], digits = 2, row.names = FALSE)
```

| cluster | nA | nB | median_cells | gini_cells | cells_ratio | cells_cond_cor | n_testable | expected_recall | call |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| CD4T | 12 | 12 | 90.0 | 0.00 | 1.00 | NA | 400 | 1.00 | Adequate |
| CD8T | 12 | 12 | 42.5 | 0.32 | 0.21 | -1 | 400 | 1.00 | Caution |
| NK | 12 | 12 | 20.0 | 0.00 | 1.00 | NA | 400 | 0.45 | Underpowered |

Each non-adequate call explains itself:

``` r

d <- as.data.frame(res)
for (i in seq_len(nrow(d))) {
    cat(d$cluster[i], "->", as.character(d$call[i]), "::", d$reason[i], "\n")
}
#> CD4T -> Adequate :: adequate 
#> CD8T -> Caution :: cell-count imbalance correlated with condition (r -1.00, cutoff 0.50 at 24 samples; consider sample weights) 
#> NK -> Underpowered :: underpowered (expected recall 45% at log2FC 1)
```

Here CD4 T is clean (**Adequate**); CD8 T has far fewer cells in
condition B, so its cell-count imbalance tracks the condition and breaks
the equal-variance assumption (**Caution** — consider voom
`sampleWeights` or dreamlet precision weights); NK has high
between-donor variance at a modest cell count, so few genes are
detectable at the target effect (**Underpowered**).

## The traffic-light matrix

[`plotAdequacy()`](https://smorabit.github.io/seny/reference/plotAdequacy.md)
renders the table as a clusters-by-contrasts grid, one plot per
annotation set, coloured by the call — the at-a-glance view.

``` r

plots <- plotAdequacy(res)
plots[["cell_type"]]
```

![Traffic-light matrix of adequacy calls by cell
type](seny_files/figure-html/matrix-1.png)

## Reading the diagnostics

A family of plotting functions explains the calls.
[`plotVariance()`](https://smorabit.github.io/seny/reference/plotVariance.md)
shows where each comparison’s noise comes from — the depth, cell and
donor shares of the per-sample variance. A donor-dominated bar means
only more donors will help.

``` r

plotVariance(res)
```

![Variance attribution bars](seny_files/figure-html/variance-1.png)

[`plotRegime()`](https://smorabit.github.io/seny/reference/plotRegime.md)
places the same decomposition on one triangle, coloured by the call and
sized by the median detectable effect: position shows where the noise
is, size shows how hard detection is.

``` r

plotRegime(res)
```

![Variance regime ternary](seny_files/figure-html/regime-1.png)

[`plotComparison()`](https://smorabit.github.io/seny/reference/plotComparison.md)
assembles a single comparison — the call and its reason, the variance
donut, the sensitivity curve and the cells per sample — into one card.

``` r

plotComparison(res, "cell_type:NK:B_vs_A")
```

![Comparison card for one comparison](seny_files/figure-html/card-1.png)

Other views:
[`plotConfounding()`](https://smorabit.github.io/seny/reference/plotConfounding.md)
(Layer-1 association),
[`plotCells()`](https://smorabit.github.io/seny/reference/plotCells.md)
(cells per sample by condition),
[`plotSensitivity()`](https://smorabit.github.io/seny/reference/plotSensitivity.md)
(detectable effect vs expression) and
[`plotPower()`](https://smorabit.github.io/seny/reference/plotPower.md)
(power vs effect size).

## Composition: subtype shifts inside a cluster

A broad cluster often holds subtypes, and their shares can differ
between the conditions. Pseudobulk DE on the broad cluster then reports
the subtypes’ marker genes as DE although no gene changes within any
subtype. CD4 T looked **Adequate** above, but its naive and memory
subtypes shift with the condition:

``` r

table(sce$cell_type_fine, sce$condition)
#>             
#>                A   B
#>   CD4 memory 432 708
#>   CD4 naive  648 372
#>   CD8T       840 180
#>   NK         240 240
```

Give
[`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
both annotation levels. `seny` works out which level nests in which from
the cells, and for each broad cluster:

- tests the subtypes’ shares of the cluster with
  [`speckle::propeller`](https://rdrr.io/pkg/speckle/man/propeller.html);
- predicts, per gene, the log2 fold change the observed shift alone
  would cause, and the power to call it DE, which sums to the number of
  DEGs expected from composition alone;
- labels each DEG of the broad cluster’s own test as
  *composition-driven*, *within-child DE*, *mixed* or *unresolved*, by
  comparing the observed change with the prediction and checking whether
  any subtype shows the gene as DE.

``` r

res2 <- assessAll(sce, sample = "donor", group = "condition",
                  annotations = c("cell_type", "cell_type_fine"), targetLFC = 1)
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
d2 <- as.data.frame(res2)
cols <- c("annotation", "cluster", "call", "composition", "comp_child",
          "exp_comp_deg", "n_parent_deg", "n_comp_driven")
knitr::kable(d2[d2$annotation == "cell_type", cols], digits = 2, row.names = FALSE)
```

| annotation | cluster | call | composition | comp_child | exp_comp_deg | n_parent_deg | n_comp_driven |
|:---|:---|:---|:---|:---|---:|---:|---:|
| cell_type | CD4T | Caution | assessed | CD4 memory | 44.36 | 70 | 40 |
| cell_type | CD8T | Caution | single child | NA | NA | NA | NA |
| cell_type | NK | Underpowered | single child | NA | NA | NA | NA |

``` r

d2$reason[d2$comparison == "cell_type:CD4T:B_vs_A"]
#> [1] "composition shift within the cluster (CD4 memory 40% -> 66%, FDR 0): ~44 genes (11.1%) expected to look DE from the shift alone"
```

CD4 T is now **Caution**: the memory share rises from about 40% to 66%,
enough to make about 44 genes look DE at the broad level.
[`plotComposition()`](https://smorabit.github.io/seny/reference/plotComposition.md)
shows the shift and which of the cluster’s DEGs it explains:
composition-driven genes sit on the diagonal, while the genes that
really change inside the subtypes do not.

``` r

plotComposition(res2, "cell_type:CD4T:B_vs_A")
```

![Subtype shares and composition attribution for CD4
T](seny_files/figure-html/compplot-1.png)

What to do with a composition Caution: test the subtypes separately
(their own calls are in the same table), or interpret the broad
cluster’s DEGs with the attribution in hand. A cluster whose overall
share differs between conditions is reported (`share_A`, `share_B`,
`share_fdr`) but not flagged, since that alone does not bias DE within
it. With a single annotation level the composition is `"not assessed"`.

## Thresholds and honest reporting

The numeric cut-offs that drive the call live in
[`defaultThresholds()`](https://smorabit.github.io/seny/reference/defaultThresholds.md)
and can be overridden. They are **provisional placeholders** — the final
values are being derived from a validation study (OneK1K, AIDA and an
in-house atlas), not chosen up front.

``` r

str(defaultThresholds())
#> List of 10
#>  $ cramerV       : num 0.3
#>  $ vif           : num 5
#>  $ sensitivity   : num 0.5
#>  $ minSamples    : num 2
#>  $ cellsCondCor  : num 0.5
#>  $ cellsCondAlpha: num 0.05
#>  $ minTestable   : num 200
#>  $ compFDR       : num 0.05
#>  $ compFrac      : num 0.01
#>  $ compExplained : num 0.5
```

This run assessed 3 comparisons. When you scan many comparisons at once,
treat the matrix as a screen, not a ranking: picking the standout cell
after the fact inflates false positives. The call tells you which
comparisons are *worth* testing — it does not replace multiple-testing
control across them.

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] seny_0.0.1
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3          jsonlite_2.0.0             
#>   [3] magrittr_2.0.5              ggbeeswarm_0.7.3           
#>   [5] spatstat.utils_3.2-5        farver_2.1.2               
#>   [7] rmarkdown_2.32              fs_2.1.0                   
#>   [9] ragg_1.5.2                  vctrs_0.7.3                
#>  [11] ROCR_1.0-12                 spatstat.explore_3.8-3     
#>  [13] DelayedMatrixStats_1.34.0   htmltools_0.5.9            
#>  [15] S4Arrays_1.12.1             SparseArray_1.12.3         
#>  [17] sass_0.4.10                 sctransform_0.4.3          
#>  [19] parallelly_1.50.0           KernSmooth_2.23-27         
#>  [21] bslib_0.12.0                htmlwidgets_1.6.4          
#>  [23] desc_1.4.3                  ica_1.0-3                  
#>  [25] plyr_1.8.9                  speckle_1.12.0             
#>  [27] plotly_4.12.1               zoo_1.9-1                  
#>  [29] cachem_1.1.0                igraph_2.3.4               
#>  [31] mime_0.13                   lifecycle_1.0.5            
#>  [33] pkgconfig_2.0.3             Matrix_1.7-6               
#>  [35] R6_2.6.1                    fastmap_1.2.0              
#>  [37] MatrixGenerics_1.24.0       fitdistrplus_1.2-6         
#>  [39] future_1.82.0               shiny_1.14.0               
#>  [41] digest_0.6.39               patchwork_1.3.2            
#>  [43] S4Vectors_0.50.3            tensor_1.5.1               
#>  [45] Seurat_5.6.0                RSpectra_0.16-2            
#>  [47] irlba_2.4.1                 textshaping_1.0.5          
#>  [49] GenomicRanges_1.64.0        beachmat_2.28.0            
#>  [51] labeling_0.4.3              progressr_1.1.0            
#>  [53] spatstat.sparse_3.2-0       polyclip_1.10-7            
#>  [55] httr_1.4.9                  abind_1.4-8                
#>  [57] compiler_4.6.2              withr_3.0.3                
#>  [59] S7_0.2.2                    BiocParallel_1.46.0        
#>  [61] fastDummies_1.7.6           MASS_7.3-66                
#>  [63] DelayedArray_0.38.2         tools_4.6.2                
#>  [65] vipor_0.4.7                 lmtest_0.9-40              
#>  [67] otel_0.2.0                  beeswarm_0.4.0             
#>  [69] httpuv_1.6.17               future.apply_1.21.0        
#>  [71] goftest_1.2-3               glmGamPoi_1.24.0           
#>  [73] glue_1.8.1                  nlme_3.1-171               
#>  [75] promises_1.5.0              grid_4.6.2                 
#>  [77] Rtsne_0.17                  cluster_2.1.8.3            
#>  [79] reshape2_1.4.5              generics_0.1.4             
#>  [81] spatstat.data_3.1-9         gtable_0.3.6               
#>  [83] tidyr_1.3.2                 data.table_1.18.6.1        
#>  [85] sp_2.2-4                    XVector_0.52.0             
#>  [87] spatstat.geom_3.8-3         BiocGenerics_0.58.1        
#>  [89] RcppAnnoy_0.0.23            ggrepel_0.9.8              
#>  [91] RANN_2.6.3                  pillar_1.11.1              
#>  [93] stringr_1.6.0               spam_2.11-4                
#>  [95] RcppHNSW_0.7.0              limma_3.68.5               
#>  [97] later_1.4.8                 splines_4.6.2              
#>  [99] dplyr_1.2.1                 lattice_0.23-1             
#> [101] deldir_2.0-4                survival_3.8-12            
#> [103] tidyselect_1.2.1            SingleCellExperiment_1.34.0
#> [105] locfit_1.5-9.12             miniUI_0.1.2               
#> [107] pbapply_1.7-5               knitr_1.52                 
#> [109] gridExtra_2.3.1             IRanges_2.46.0             
#> [111] Seqinfo_1.2.0               edgeR_4.10.7               
#> [113] SummarizedExperiment_1.42.0 scattermore_1.2            
#> [115] stats4_4.6.2                xfun_0.62                  
#> [117] Biobase_2.72.0              statmod_1.5.2              
#> [119] matrixStats_1.5.0           stringi_1.8.9              
#> [121] yaml_2.3.12                 evaluate_1.0.5             
#> [123] codetools_0.2-20            tibble_3.3.1               
#> [125] cli_3.6.6                   uwot_0.2.5                 
#> [127] xtable_1.8-8                reticulate_1.47.0          
#> [129] systemfonts_1.3.2           jquerylib_0.1.4            
#> [131] Rcpp_1.1.2                  spatstat.random_3.5-2      
#> [133] globals_0.20.0              png_0.1-9                  
#> [135] spatstat.univar_3.2-0       parallel_4.6.2             
#> [137] pkgdown_2.2.1               ggplot2_4.0.3              
#> [139] dotCall64_1.2               sparseMatrixStats_1.24.0   
#> [141] listenv_1.1.0               viridisLite_0.4.3          
#> [143] scales_1.4.0                ggridges_0.5.7             
#> [145] SeuratObject_5.4.0          purrr_1.2.2                
#> [147] rlang_1.3.0                 cowplot_1.2.0
```
