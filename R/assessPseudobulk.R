#' Layer 2: pseudobulk quality
#'
#' For every comparison, summarises the pseudobulks the test will run on: the
#' cells contributed per sample and their inequality, whether cell counts or
#' library sizes track the condition (which breaks the equal-variance
#' assumption), and how many genes are testable. Specifically, per comparison:
#' cells per sample (`min_cells`, `median_cells`, `gini_cells`, and
#' `frac_below_min`, the fraction of the cluster's samples dropped for having
#' fewer than `minCells` cells); group imbalance (`cells_ratio`, the smaller
#' group median over the larger, and `cells_cond_cor`, the correlation of
#' log cells with condition); `libsize_cond_cor` (correlation of log library
#' size with condition); and testable genes (`n_testable` passing
#' [edgeR::filterByExpr] under the comparison's design, and `frac_testable`).
#' The testable-gene count is the denominator for the sensitivity fraction in
#' [estimateSensitivity()]. Outlier-sample detection is part of layer 3, not
#' here.
#'
#' @param x A [S4Vectors::SimpleList] from [pbSetup()].
#' @return A `data.frame`, one row per comparison, with `min_cells`,
#'   `median_cells`, `gini_cells`, `frac_below_min`, `cells_ratio`,
#'   `cells_cond_cor`, `libsize_cond_cor`, `n_testable` and `frac_testable`.
#' @export
#' @importFrom edgeR DGEList filterByExpr
assessPseudobulk <- function(x) {
    do.call(rbind, lapply(x, function(se) {
        md <- S4Vectors::metadata(se)
        cc <- md$cell_counts
        cmp <- enumerate_comparisons(se)
        do.call(rbind, lapply(seq_len(nrow(cmp)), function(i) {
            d <- comparison_data(se, cmp$cluster[i], cmp$A[i], cmp$B[i])
            key <- data.frame(comparison = cmp$comparison[i], stringsAsFactors = FALSE)
            if (is.null(d) || length(d$group) < 2) {
                return(cbind(key, min_cells = NA_integer_, median_cells = NA_real_,
                             gini_cells = NA_real_, frac_below_min = NA_real_,
                             cells_ratio = NA_real_, cells_cond_cor = NA_real_,
                             libsize_cond_cor = NA_real_, n_testable = 0L,
                             frac_testable = NA_real_))
            }
            cond <- as.integer(d$group == 'B')
            med_a <- stats::median(d$nCells[cond == 0])
            med_b <- stats::median(d$nCells[cond == 1])
            cells_ratio <- if (max(med_a, med_b) > 0) min(med_a, med_b) / max(med_a, med_b) else 0

            # fraction of this cluster's samples (both groups, pre-filter) below minCells
            below <- cc[cc$cluster == cmp$cluster[i] &
                        cc$group %in% c(cmp$A[i], cmp$B[i]), ]
            frac_below_min <- if (nrow(below)) mean(below$nCells < md$minCells) else NA_real_

            bd <- build_design(d$group, d$covariates)
            n_testable <- if (!bd$estimable) {
                NA_integer_
            } else if (identical(md$aggregation, 'mean')) {
                sum(edgeR::filterByExpr(d$counts, bd$design, lib.size = colSums(d$counts)))
            } else {
                sum(edgeR::filterByExpr(edgeR::DGEList(d$counts), bd$design))
            }

            cbind(key, min_cells = min(d$nCells), median_cells = stats::median(d$nCells),
                  gini_cells = gini(d$nCells), frac_below_min = frac_below_min,
                  cells_ratio = cells_ratio,
                  cells_cond_cor = point_biserial(log(d$nCells), cond),
                  libsize_cond_cor = point_biserial(log(colSums(d$counts)), cond),
                  n_testable = n_testable,
                  frac_testable = if (is.na(n_testable)) NA_real_ else n_testable / nrow(d$counts))
        }))
    }))
}

#' Correlation of a numeric vector with a 0/1 condition, NA when undefined
#' @noRd
point_biserial <- function(v, cond) {
    if (length(unique(cond)) < 2 || stats::sd(v) == 0) return(NA_real_)
    suppressWarnings(stats::cor(v, cond))
}
