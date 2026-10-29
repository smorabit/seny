#' Cells per sample by condition
#'
#' The cells contributed per sample, split by condition, as a beeswarm over a
#' boxplot and faceted by cluster. When cell counts track the condition the
#' pseudobulks are heteroscedastic, which is the risk behind a cell-count
#' imbalance Caution.
#'
#' @param res The result of [assessAll()] (its `pb` setup object is required).
#' @return A `ggplot`, faceted by cluster (and annotation set when there is more
#'   than one).
#' @export
#' @importFrom ggplot2 .data
plotCells <- function(res) {
    pb <- require_attr(res, 'pb')
    cd <- do.call(rbind, lapply(as.list(pb), function(se) {
        d <- as.data.frame(SummarizedExperiment::colData(se))
        d$annotation <- S4Vectors::metadata(se)$annotation
        d[, c('annotation', 'cluster', 'group', 'nCells')]
    }))
    cd$group <- factor(cd$group)
    cols <- if (nlevels(cd$group) <= 2) seny_cond_colors[seq_len(nlevels(cd$group))]
            else scales::hue_pal()(nlevels(cd$group))
    p <- ggplot2::ggplot(cd, ggplot2::aes(x = .data$group, y = .data$nCells, color = .data$group)) +
        ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5) +
        ggbeeswarm::geom_quasirandom(width = 0.25, size = 1, alpha = 0.6) +
        ggplot2::scale_color_manual(values = cols, guide = 'none') +
        ggplot2::labs(x = 'condition', y = 'cells per sample', title = 'Cells per sample by condition') +
        seny_theme()
    if (length(unique(cd$annotation)) > 1) {
        p + ggplot2::facet_grid(annotation ~ cluster, scales = 'free')
    } else {
        p + ggplot2::facet_wrap(~ cluster)
    }
}
