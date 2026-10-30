#---------------------------------------------------------#
# Comparison enumeration and per-comparison subsetting (internal)
#
# A comparison is one (annotation level, cluster, contrast) triple. The layer
# functions iterate these, subsetting the pseudobulk SummarizedExperiment to the
# cluster's samples in the two contrast groups.
#---------------------------------------------------------#

#' Enumerate the comparisons in one annotation-level SummarizedExperiment
#' @return data.frame(annotation, cluster, A, B, comparison)
#' @noRd
enumerate_comparisons <- function(se) {
    md <- S4Vectors::metadata(se)
    # include clusters that lost every pseudobulk to minCells, so they are called
    # Invalid rather than silently dropped
    clusters <- unique(c(as.character(SummarizedExperiment::colData(se)$cluster),
                         as.character(md$cell_counts$cluster)))
    rows <- list()
    for (cl in clusters) for (ct in md$contrasts) {
        rows[[length(rows) + 1]] <- data.frame(
            annotation = md$annotation, cluster = cl, A = ct[1], B = ct[2],
            comparison = paste0(md$annotation, ':', cl, ':', ct[2], '_vs_', ct[1]),
            stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
}

#' Pull the pseudobulk data for one comparison
#' @return a list with counts, group (factor A/B), covariates, nCells, theta and
#'   samples, or NULL if the cluster/contrast is absent
#' @noRd
comparison_data <- function(se, cluster, A, B) {
    md <- S4Vectors::metadata(se)
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    sel <- cd$cluster == cluster & cd$group %in% c(A, B)
    if (sum(sel) == 0) return(NULL)
    counts <- as.matrix(SummarizedExperiment::assay(se, 'counts')[, sel, drop = FALSE])
    grp <- factor(ifelse(cd$group[sel] == B, 'B', 'A'), levels = c('A', 'B'))
    covs <- if (length(md$covariates)) cd[sel, md$covariates, drop = FALSE] else NULL
    list(counts = counts, group = grp, covariates = covs, nCells = cd$nCells[sel],
         theta = md$theta[[cluster]], samples = rownames(cd)[sel])
}
