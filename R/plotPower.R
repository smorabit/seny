#' Power vs effect size
#'
#' Mean power across a sweep of target log2 fold-changes, one line per cluster,
#' with the 0.8 reference. Shows where the call would change at a smaller or larger
#' target effect.
#'
#' @param res The result of [assessAll()] (its per-gene `pergene` table is
#'   required).
#' @param targets Numeric vector of target absolute log2 fold-changes to sweep.
#' @return A `ggplot`, faceted by annotation set when there is more than one.
#' @export
#' @importFrom ggplot2 .data
plotPower <- function(res, targets = seq(0.25, 2, by = 0.05)) {
    pg <- require_attr(res, 'pergene')
    meta <- as.data.frame(res)[, c('comparison', 'cluster', 'annotation')]
    pg <- merge(pg, meta, by = 'comparison')
    curve <- do.call(rbind, lapply(targets, function(t) {
        pg$pw <- power_from_se(pg$se_log2, pg$df, t)
        agg <- stats::aggregate(pw ~ cluster + annotation, data = pg, FUN = mean)
        agg$lfc <- t
        agg
    }))
    p <- ggplot2::ggplot(curve, ggplot2::aes(x = .data$lfc, y = .data$pw, color = .data$cluster)) +
        ggplot2::geom_hline(yintercept = 0.8, linetype = 'dotted', colour = 'grey40') +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::labs(x = 'target |log2FC|', y = 'mean power', title = 'Power vs effect size') +
        seny_theme()
    facet_by_annotation(p, curve)
}
