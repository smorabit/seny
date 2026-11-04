#' Default decision thresholds for the adequacy call
#'
#' Values used by [callComparisons()], each from a validation experiment or a
#' stated convention (threshold audit before 1.0.0); a user can override any.
#'
#' @return A named list: `cramerV` (condition-covariate association above which a
#'   comparison is flagged Caution; on null splits with the confounder ignored,
#'   false discoveries rise above the base rate from 0.3), `vif` (variance
#'   inflation factor of the condition given the covariates at or above which the
#'   comparison is flagged Caution; a convention, R^2 = 0.8: adjusted fits stay
#'   calibrated and the power cost is already in the expected recall, but the
#'   estimate rests on the few samples that break the collinearity),
#'   `sensitivity` (expected recall at the target effect at or above which a
#'   comparison is adequately powered; 0.5 matches replication failure defined as
#'   recall below 0.5), `minSamples` (fewer samples per group is Invalid),
#'   `cellsCondCor` (absolute correlation of log cells per sample with condition
#'   above which cell-count imbalance is reported, and flagged Caution under
#'   DESeq2, whose FDR control it breaks), `cellsCondAlpha`
#'   (the correlation must also exceed the critical value at this level for the
#'   comparison's sample size, so chance correlations in small designs are not
#'   flagged: the cutoff is `max(cellsCondCor, r_crit(n))`, about 0.81 at 6
#'   samples and `cellsCondCor` from ~16 samples), and `minTestable` (fewer
#'   testable genes is Underpowered; a convention guarding the genome-wide
#'   estimates behind the call, pi1 from BH discoveries and the dispersion prior;
#'   power stays calibrated down to 100-200 testable genes). The composition values (see
#'   [assessComposition()]) were calibrated on imposed subtype shifts in
#'   OneK1K (E6): `compFDR` (propeller FDR below which a child's share within
#'   its parent is called shifted), `compFrac` (expected DEGs from that shift
#'   alone, as a share of testable genes, at or above which the comparison is
#'   flagged Caution; 0.01 maximised Youden's J for predicting artefacts in at
#'   least 1% of genes, with no flags without a shift) and `compExplained`
#'   (share of a parent DEG's observed log2 fold change the predicted shift-only
#'   change must explain for the gene to be attributed to composition; 0.5:
#'   on natural age shifts in OneK1K (E5) it labels 3% of real within-subtype
#'   DEGs composition-driven, against 8% at 0.25, the best value on imposed
#'   shifts, while 0.99 of its composition-driven labels are artefacts).
#' @examples
#' th <- defaultThresholds()
#' th$sensitivity
#' # override one threshold, keep the rest
#' str(utils::modifyList(th, list(sensitivity = 0.7)))
#' @export
defaultThresholds <- function() {
    list(cramerV = 0.3, vif = 5, sensitivity = 0.5, minSamples = 2,
         cellsCondCor = 0.5, cellsCondAlpha = 0.05, minTestable = 200,
         compFDR = 0.05, compFrac = 0.01, compExplained = 0.5)
}

#' Cutoff for the cell-count-vs-condition correlation at n samples: the larger of
#' the fixed floor and the two-sided critical Pearson r at level alpha
#' @noRd
cells_cond_cutoff <- function(n, floor, alpha = NULL) {
    if (is.null(alpha) || is.na(n)) return(floor)
    if (n < 3) return(Inf)
    tc <- stats::qt(1 - alpha / 2, n - 2)
    max(floor, tc / sqrt(n - 2 + tc^2))
}

#' One line on where each default threshold comes from (the threshold audit
#' before 0.99.0), for the report's appendix
#' @noRd
threshold_sources <- function() {
    c(cramerV = 'null splits: false discoveries rise above the base rate from 0.3 (E1)',
      vif = 'convention (R^2 = 0.8); adjusted fits stay calibrated, power cost is in expected recall',
      sensitivity = 'matches replication failure, recall < 0.5 (E4c)',
      minSamples = 'a group variance needs at least 2 samples',
      cellsCondCor = 'imbalance breaks DESeq2 FDR control, not edgeR or voom (E7); Caution under DESeq2 only',
      cellsCondAlpha = 'chance correlations in small designs are not flagged (E4)',
      minTestable = 'convention guarding the genome-wide estimates (pi1, dispersion prior)',
      compFDR = 'FDR convention for the subtype shift test',
      compFrac = 'imposed subtype shifts: Youden J 0.94 (E6)',
      compExplained = 'natural age shifts: 3% of real DEGs mislabelled (E5)')
}
