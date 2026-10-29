# seny

**Sanity checks for differential gene Expression in NoisY data** (and *seny*, the
Catalan word for common sense).

For every cluster-by-contrast comparison in a multi-sample single-cell dataset,
`seny` decides whether a pseudobulk differential-expression test is **Adequate**,
**Underpowered**, **Caution** or **Invalid**, and explains why. The diagnostics
are analytic and run on your own data — there is no per-dataset simulation.

| Call | Meaning |
| --- | --- |
| **Invalid** | Type I error cannot be trusted (design not estimable, too few samples per group). |
| **Caution** | Valid, but a risk applies (condition confounded with a covariate, or cell-count imbalance tracking the condition). |
| **Underpowered** | Valid but unable to detect the target effect (sensitivity too low, or too few testable genes). |
| **Adequate** | No flags. |

The call has a severity order — Invalid > Caution > Underpowered > Adequate — and a
`reason` string that lists every flag, so a single headline never hides a second
problem.

## Installation

`seny` is not yet on Bioconductor. Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("smorabit/seny")
```

## Usage

```r
library(seny)

sce <- readRDS(system.file("extdata", "example_sce.rds", package = "seny"))

res <- assessAll(sce, sample = "donor", group = "condition",
                 annotations = "cell_type", targetLFC = 1)
res                         # one row per comparison: metrics, call, reason

plots <- plotAdequacy(res)  # traffic-light matrix, one ggplot per annotation set
plots[["cell_type"]]
```

See `vignette("seny")` for the full walk-through.

## Status

Experimental (P1 MVP). The package covers diagnostic layers 1, 2 and 4; the call
thresholds in `defaultThresholds()` are **provisional placeholders** being
recalibrated from a validation study (OneK1K, AIDA and an in-house atlas). Outlier
detection, empirical checks, composition assessment, calibrated power and an HTML
report are planned.
