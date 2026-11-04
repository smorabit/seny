#' Assemble the comparison table and make the adequacy call
#'
#' Merges the layer metrics into one tidy table, one row per comparison, and
#' assigns each comparison a call: **Invalid** (the contrast is not estimable, or
#' too few samples per group), **Caution** (estimable but the condition is
#' confounded with a covariate, cell-count imbalance tracks the condition under
#' DESeq2 (whose FDR control it breaks; with edgeR or voom the imbalance is
#' reported in the reason but does not change the call), or a
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
#'   the `call`, a `reason` string and `flags`, the same flags as short codes
#'   separated by `;` (`confounded`, `imbalance`, `composition`,
#'   `few_testable`, `underpowered`, `too_few_samples`, `not_estimable`; empty
#'   when adequate).
#' @examples
#' sce <- readRDS(system.file("extdata", "example_sce.rds", package = "seny"))
#' x <- pbSetup(sce, sample = "donor", group = "condition",
#'              annotations = "cell_type", thetaReps = 0)
#' variance <- estimateVariance(x)
#' sens <- estimateSensitivity(variance, targetLFC = 1)
#' res <- callComparisons(assessDesign(x), assessPseudobulk(x), variance$summary, sens)
#' res[, c("comparison", "call", "reason")]
#' @export
#' @importFrom S4Vectors DataFrame
callComparisons <- function(design, pseudobulk, variance, sensitivity,
                            thresholds = defaultThresholds(), composition = NULL) {
    parts <- list(design, pseudobulk, variance, sensitivity)
    if (!is.null(composition)) parts <- c(parts, list(as.data.frame(composition)))
    tab <- Reduce(function(a, b) merge(a, b, by = 'comparison', all = TRUE), parts)
    th <- utils::modifyList(defaultThresholds(), thresholds)
    call <- character(nrow(tab)); reason <- character(nrow(tab)); codes <- character(nrow(tab))

    for (i in seq_len(nrow(tab))) {
        r <- tab[i, ]
        flags <- character(0)
        confounded <- (!is.na(r$max_cramerV) && r$max_cramerV >= th$cramerV) ||
                      (!is.na(r$max_vif) && r$max_vif >= th$vif)
        cells_cut <- cells_cond_cutoff(r$nA + r$nB, th$cellsCondCor, th$cellsCondAlpha)
        imbalanced <- !is.na(r$cells_cond_cor) && abs(r$cells_cond_cor) >= cells_cut
        # imbalance breaks DESeq2's FDR control; edgeR and voom stay calibrated
        # (threshold audit on E7), so for them it is reported, not a Caution
        method <- if (is.null(r$method) || is.na(r$method)) NA_character_ else as.character(r$method)
        imbalance_caution <- imbalanced && (is.na(method) || method == 'deseq2')
        few_testable <- !is.na(r$n_testable) && r$n_testable < th$minTestable
        underpowered <- !is.na(r$expected_recall) && r$expected_recall < th$sensitivity
        too_few <- is.na(r$estimable) || !r$estimable ||
                   min(r$nA, r$nB) < th$minSamples
        comp_shift <- !is.null(r$comp_min_fdr) && !is.na(r$comp_min_fdr) &&
                      r$comp_min_fdr < th$compFDR && !is.na(r$frac_comp_expected) &&
                      r$frac_comp_expected >= th$compFrac

        if (confounded) flags <- c(flags, sprintf('confounded (Cramer V %.2f, VIF %.1f)',
                                                  r$max_cramerV, r$max_vif))
        if (imbalanced) flags <- c(flags, sprintf('cell-count imbalance correlated with condition (r %.2f, cutoff %.2f at %d samples; %s)',
                                                  r$cells_cond_cor, cells_cut, r$nA + r$nB,
                                                  if (imbalance_caution) 'DESeq2 loses FDR control here; consider edgeR or voom'
                                                  else 'reported only: edgeR and voom stay calibrated'))
        if (comp_shift) flags <- c(flags, sprintf('composition shift within the cluster (%s %.0f%% -> %.0f%%, FDR %.2g): ~%.0f genes (%.1f%%) expected to look DE from the shift alone',
                                                  r$comp_child, 100 * r$comp_share_A, 100 * r$comp_share_B,
                                                  r$comp_min_fdr, r$exp_comp_deg,
                                                  100 * r$frac_comp_expected))
        if (few_testable) flags <- c(flags, sprintf('too few testable genes (%d)', r$n_testable))
        if (underpowered) flags <- c(flags, sprintf('underpowered (expected recall %.0f%% at log2FC %.2g)',
                                                    100 * r$expected_recall, r$target_lfc))
        if (too_few) flags <- c(flags, if (isTRUE(r$estimable)) 'too few samples per group'
                                       else 'design not estimable (rank-deficient or nested)')

        codes[i] <- paste(c(if (confounded) 'confounded', if (imbalanced) 'imbalance',
                            if (comp_shift) 'composition', if (few_testable) 'few_testable',
                            if (underpowered) 'underpowered',
                            if (too_few) (if (isTRUE(r$estimable)) 'too_few_samples' else 'not_estimable')),
                          collapse = ';')
        call[i] <- if (too_few) 'Invalid'
                   else if (confounded || imbalance_caution || comp_shift) 'Caution'
                   else if (underpowered || few_testable) 'Underpowered'
                   else 'Adequate'
        reason[i] <- if (length(flags)) paste(flags, collapse = '; ') else 'adequate'
    }
    tab$call <- factor(call, levels = c('Adequate', 'Underpowered', 'Caution', 'Invalid'))
    tab$reason <- reason
    tab$flags <- codes
    S4Vectors::DataFrame(tab)
}
