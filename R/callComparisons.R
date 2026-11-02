#' Assemble the comparison table and make the adequacy call
#'
#' Merges the layer metrics into one tidy table, one row per comparison, and
#' assigns each comparison a call: **Invalid** (the contrast is not estimable, or
#' too few samples per group), **Caution** (estimable but the condition is
#' confounded with a covariate, cell-count imbalance tracks the condition, or a
#' subtype shift within the cluster would by itself make enough genes look DE;
#' see [assessComposition()]),
#' **Underpowered** (valid but the expected recall at the target effect is below
#' `thresholds$sensitivity`, or too few testable genes), or **Adequate**. Severity order is
#' Invalid > Caution > Underpowered > Adequate; the `reason` column lists every
#' flag so the single call never hides a second problem.
#'
#' @param design,pseudobulk,variance,sensitivity The outputs of [assessDesign()],
#'   [assessPseudobulk()], [estimateVariance()] (its `$summary`) and
#'   [estimateSensitivity()].
#' @param composition Optional `$summary` of [assessComposition()]. `NULL`
#'   (default) skips the composition columns and flag.
#' @param thresholds A list from [defaultThresholds()].
#' @return A [S4Vectors::DataFrame], one row per comparison, with all metrics,
#'   the `call` and a `reason` string.
#' @export
#' @importFrom S4Vectors DataFrame
callComparisons <- function(design, pseudobulk, variance, sensitivity,
                            thresholds = defaultThresholds(), composition = NULL) {
    parts <- list(design, pseudobulk, variance, sensitivity)
    if (!is.null(composition)) parts <- c(parts, list(as.data.frame(composition)))
    tab <- Reduce(function(a, b) merge(a, b, by = 'comparison', all = TRUE), parts)
    th <- utils::modifyList(defaultThresholds(), thresholds)
    call <- character(nrow(tab)); reason <- character(nrow(tab))

    for (i in seq_len(nrow(tab))) {
        r <- tab[i, ]
        flags <- character(0)
        confounded <- (!is.na(r$max_cramerV) && r$max_cramerV >= th$cramerV) ||
                      (!is.na(r$max_vif) && r$max_vif >= th$vif)
        cells_cut <- cells_cond_cutoff(r$nA + r$nB, th$cellsCondCor, th$cellsCondAlpha)
        imbalanced <- !is.na(r$cells_cond_cor) && abs(r$cells_cond_cor) >= cells_cut
        few_testable <- !is.na(r$n_testable) && r$n_testable < th$minTestable
        underpowered <- !is.na(r$expected_recall) && r$expected_recall < th$sensitivity
        too_few <- is.na(r$estimable) || !r$estimable ||
                   min(r$nA, r$nB) < th$minSamples
        comp_shift <- !is.null(r$comp_min_fdr) && !is.na(r$comp_min_fdr) &&
                      r$comp_min_fdr < th$compFDR && !is.na(r$frac_comp_sensitive) &&
                      r$frac_comp_sensitive >= th$compFrac

        if (confounded) flags <- c(flags, sprintf('confounded (Cramer V %.2f, VIF %.1f)',
                                                  r$max_cramerV, r$max_vif))
        if (imbalanced) flags <- c(flags, sprintf('cell-count imbalance correlated with condition (r %.2f, cutoff %.2f at %d samples; consider sample weights)',
                                                  r$cells_cond_cor, cells_cut, r$nA + r$nB))
        if (comp_shift) flags <- c(flags, sprintf('composition shift within the cluster (%s %.0f%% -> %.0f%%, FDR %.2g): %d genes (%.1f%%) would look DE from the shift alone',
                                                  r$comp_child, 100 * r$comp_share_A, 100 * r$comp_share_B,
                                                  r$comp_min_fdr, as.integer(r$n_comp_sensitive),
                                                  100 * r$frac_comp_sensitive))
        if (few_testable) flags <- c(flags, sprintf('too few testable genes (%d)', r$n_testable))
        if (underpowered) flags <- c(flags, sprintf('underpowered (expected recall %.0f%% at log2FC %.2g)',
                                                    100 * r$expected_recall, r$target_lfc))
        if (too_few) flags <- c(flags, if (isTRUE(r$estimable)) 'too few samples per group'
                                       else 'design not estimable (rank-deficient or nested)')

        call[i] <- if (too_few) 'Invalid'
                   else if (confounded || imbalanced || comp_shift) 'Caution'
                   else if (underpowered || few_testable) 'Underpowered'
                   else 'Adequate'
        reason[i] <- if (length(flags)) paste(flags, collapse = '; ') else 'adequate'
    }
    tab$call <- factor(call, levels = c('Adequate', 'Underpowered', 'Caution', 'Invalid'))
    tab$reason <- reason
    S4Vectors::DataFrame(tab)
}
