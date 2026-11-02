#' Composition of one comparison: subtype shares and composition-driven DE
#'
#' Two panels for a comparison whose cluster has a finer annotation nested in
#' it (see [assessComposition()]). Left: each subtype's share of the cluster's
#' cells per sample, by condition; a shift here makes the subtypes' markers look
#' DE in the cluster's pseudobulk. Right: per gene, the observed log2 fold change
#' of the cluster's own DE test against the change the subtype shift alone
#' predicts. Genes the test calls DE are coloured by their attribution; genes on
#' the diagonal are explained by composition.
#'
#' @param res The result of [assessAll()] run with two nested annotation levels.
#' @param comparison A comparison key from `res$comparison`.
#' @param nLabel Number of composition-driven genes to label, largest observed
#'   change first. Default 10.
#' @return A `patchwork` of two `ggplot`s.
#' @export
#' @importFrom ggplot2 .data
plotComposition <- function(res, comparison, nLabel = 10) {
    comp <- require_attr(res, 'composition')
    row <- as.data.frame(res)[as.data.frame(res)$comparison == comparison, ]
    if (nrow(row) == 0) stop('unknown comparison: "', comparison, '"; see res$comparison', call. = FALSE)
    if (row$composition != 'assessed') {
        stop('composition ', row$composition, ' for "', comparison, '"; it needs a finer ',
             'annotation nested in this one (pass both levels to assessAll)', call. = FALSE)
    }
    pr <- comp$proportions[comp$proportions$comparison == comparison, ]
    pr$group <- factor(pr$group, levels = c(row$A, row$B))
    left <- ggplot2::ggplot(pr, ggplot2::aes(x = .data$child, y = .data$share, colour = .data$group)) +
        ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, position = ggplot2::position_dodge(0.75)) +
        ggbeeswarm::geom_quasirandom(dodge.width = 0.75, width = 0.12, size = 0.9, alpha = 0.6) +
        ggplot2::scale_colour_manual(values = seny_cond_colors, name = NULL) +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::labs(x = NULL, y = paste('share of', row$cluster, 'cells'), title = 'Subtype shares') +
        seny_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                       legend.position = 'top')

    g <- comp$genes[comp$genes$comparison == comparison & !is.na(comp$genes$logfc), ]
    g$label <- factor(ifelse(is.na(g$label), 'not DE', g$label),
                      levels = c('composition-driven', 'mixed', 'within-child DE', 'unresolved', 'not DE'))
    g <- g[order(g$label, decreasing = TRUE), ]
    lab <- g[g$label == 'composition-driven', ]
    lab <- utils::head(lab[order(-abs(lab$logfc)), ], nLabel)
    attr_cols <- c(`composition-driven` = '#D55E00', mixed = '#CC79A7', `within-child DE` = '#0072B2',
                   unresolved = '#999999', `not DE` = 'grey85')
    right <- ggplot2::ggplot(g, ggplot2::aes(x = .data$logfc, y = .data$lfc_comp, colour = .data$label)) +
        ggplot2::geom_hline(yintercept = 0, colour = 'grey70') +
        ggplot2::geom_vline(xintercept = 0, colour = 'grey70') +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'grey40') +
        ggplot2::geom_point(size = 1, alpha = 0.7) +
        ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = .data$gene), size = 3,
                                 colour = 'black', max.overlaps = Inf, min.segment.length = 0) +
        ggplot2::scale_colour_manual(values = attr_cols, name = NULL, drop = FALSE) +
        ggplot2::labs(x = 'observed log2FC', y = 'predicted from the subtype shift alone',
                      title = sprintf('%d of %d DEGs composition-driven', as.integer(row$n_comp_driven),
                                      as.integer(row$n_parent_deg))) +
        seny_theme() +
        ggplot2::theme(legend.position = 'top')
    patchwork::wrap_plots(left, right, nrow = 1, widths = c(1, 1.2)) +
        patchwork::plot_annotation(title = paste0(row$cluster, ' : ', row$B, ' vs ', row$A,
                                                  '  \u2014  ', as.character(row$call)))
}
