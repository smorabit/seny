# seny 0.99.0

First Bioconductor submission.

* Per-comparison adequacy calls (Adequate / Underpowered / Caution / Invalid)
  for pseudobulk differential expression, each with its reason: design and
  confounding (layer 1), pseudobulk quality and cell imbalance (layer 2), and
  variance structure and power at a target effect (layer 4), with edgeR, DESeq2
  and voom.
* The composition layer: subtype shifts within a cluster, the expected number of
  DEGs they produce on their own, and attribution of the cluster's DEGs.
* Count matrices beyond one sparse matrix (a `DelayedMatrix`, HDF5-backed) are
  checked and summed a block of cells at a time.
* Plots: the traffic-light matrix (`plotAdequacy()`), the comparison card
  (`plotComparison()`) and the `plotX()` family.
