#' Variance attribution bars
#'
#' Stacked depth / cell / donor variance shares per comparison: where each
#' comparison's per-sample noise comes from. A donor-dominated bar means more
#' cells or reads will not help, only more donors. Comparisons fit with the voom
#' method have no negative-binomial decomposition and are dropped.
#'
#' @param res The result of [assessAll()] or [callComparisons()].
#' @return A `ggplot`, faceted by annotation set when there is more than one.
#' @export
#' @importFrom ggplot2 .data
plotVariance <- function(res) {
    df <- as.data.frame(res)
    if (any(is.na(df$share_donor))) {
        message('dropping ', sum(is.na(df$share_donor)),
                ' comparison(s) with no variance decomposition (voom method)')
        df <- df[!is.na(df$share_donor), , drop = FALSE]
    }
    if (nrow(df) == 0) stop('no variance decomposition to plot', call. = FALSE)
    long <- do.call(rbind, lapply(c('depth', 'cell', 'donor'), function(t) {
        data.frame(annotation = df$annotation, cluster = df$cluster, term = t,
                   share = df[[paste0('share_', t)]], stringsAsFactors = FALSE)
    }))
    long$term <- factor(long$term, levels = c('depth', 'cell', 'donor'))
    long$cluster <- factor(long$cluster, levels = rev(sort(unique(long$cluster))))
    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$share, y = .data$cluster, fill = .data$term)) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(values = seny_share_colors, name = 'variance term') +
        ggplot2::scale_x_continuous(labels = function(x) paste0(x * 100, '%')) +
        ggplot2::labs(x = 'share of per-sample variance', y = NULL, title = 'Variance attribution') +
        seny_theme()
    facet_by_annotation(p, long)
}
