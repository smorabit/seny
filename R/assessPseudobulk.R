#' Layer 2: pseudobulk quality (minimal)
#'
#' For every comparison, reports cells per sample (minimum and median) and the
#' number of testable genes (passing [edgeR::filterByExpr] under the comparison's
#' design). The testable-gene count is the denominator for the sensitivity
#' fraction in [estimateSensitivity()]. The full layer 2 (Gini, group imbalance,
#' library-size and outlier metrics) is added later.
#'
#' @param x A [S4Vectors::SimpleList] from [pbSetup()].
#' @return A `data.frame`, one row per comparison, with `min_cells`,
#'   `median_cells` and `n_testable`.
#' @export
#' @importFrom edgeR DGEList filterByExpr
assessPseudobulk <- function(x) {
    do.call(rbind, lapply(x, function(se) {
        cmp <- enumerate_comparisons(se)
        do.call(rbind, lapply(seq_len(nrow(cmp)), function(i) {
            d <- comparison_data(se, cmp$cluster[i], cmp$A[i], cmp$B[i])
            key <- data.frame(comparison = cmp$comparison[i], stringsAsFactors = FALSE)
            if (is.null(d) || length(d$group) < 2) {
                return(cbind(key, min_cells = NA_integer_, median_cells = NA_real_,
                             n_testable = 0L))
            }
            bd <- build_design(d$group, d$covariates)
            n_testable <- if (bd$estimable) {
                y <- edgeR::DGEList(d$counts)
                sum(edgeR::filterByExpr(y, bd$design))
            } else NA_integer_
            cbind(key, min_cells = min(d$nCells), median_cells = stats::median(d$nCells),
                  n_testable = n_testable)
        }))
    }))
}
