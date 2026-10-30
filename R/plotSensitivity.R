#' Detectable effect vs expression
#'
#' For every gene, the minimum detectable log2 fold-change against its expression,
#' with the target effect drawn in and one panel per cluster. The share of genes
#' below the line is the sensitivity fraction; genes above it are out of reach at
#' the target.
#'
#' @param res The result of [assessAll()] (its per-gene `pergene` table is
#'   required).
#' @return A `ggplot`, faceted by cluster (and annotation set when there is more
#'   than one).
#' @export
#' @importFrom ggplot2 .data
plotSensitivity <- function(res) {
    pg <- require_attr(res, 'pergene')
    meta <- as.data.frame(res)[, c('comparison', 'cluster', 'annotation', 'target_lfc', 'alpha')]
    pg <- merge(pg, meta, by = 'comparison')
    pg$mde <- mde_from_se(pg$se_log2, pg$df, pg$alpha)
    target <- stats::median(pg$target_lfc)
    ymax <- stats::quantile(pg$mde, 0.98, na.rm = TRUE)
    p <- ggplot2::ggplot(pg, ggplot2::aes(x = .data$acpm, y = .data$mde)) +
        ggplot2::geom_point(size = 0.5, alpha = 0.25, colour = '#4C72B0') +
        ggplot2::geom_hline(yintercept = target, linetype = 'dashed', colour = '#D73027') +
        ggplot2::coord_cartesian(ylim = c(0, ymax)) +
        ggplot2::labs(x = 'expression (aveLogCPM)', y = 'minimum detectable |log2FC|',
                      title = 'Detectable effect vs expression') +
        seny_theme()
    if (length(unique(pg$annotation)) > 1) {
        p + ggplot2::facet_grid(annotation ~ cluster)
    } else {
        p + ggplot2::facet_wrap(~ cluster)
    }
}
