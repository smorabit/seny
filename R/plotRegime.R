#' Variance regime ternary
#'
#' The depth / cell / donor decomposition as a single simplex: each comparison is
#' a point among the three variance terms, coloured by its call and sized by its
#' median detectable effect. Position shows where the noise is, size shows how hard
#' detection is. The axes are gridded so the shares can be read off directly.
#' Comparisons fit with the voom method have no decomposition and are dropped.
#'
#' @param res The result of [assessAll()] or [callComparisons()]. Comparisons from
#'   all annotation sets are pooled into one triangle; subset `res` first for a
#'   single annotation.
#' @return A `ggplot`.
#' @export
#' @importFrom ggplot2 .data
plotRegime <- function(res) {
    df <- as.data.frame(res)
    drop <- is.na(df$share_donor)
    if (any(drop)) {
        warning('dropping ', sum(drop), ' comparison(s) with no variance decomposition (voom method)',
                call. = FALSE)
        df <- df[!drop, , drop = FALSE]
    }
    if (nrow(df) == 0) stop('no variance decomposition to plot', call. = FALSE)

    xy <- ternary_xy(df$share_depth, df$share_cell, df$share_donor)
    df$x <- xy$x; df$y <- xy$y
    df$call <- factor(as.character(df$call), levels = names(seny_call_colors))
    df <- df[order(-df$median_mde), ]        # big points first, small ones on top
    g <- ternary_geom()

    ggplot2::ggplot() +
        ggplot2::geom_segment(data = g$grid,
            ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
            colour = 'grey82', linewidth = 0.4) +
        ggplot2::geom_text(data = g$ticks,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab, hjust = .data$hj, vjust = .data$vj),
            size = 2.8, colour = 'grey45') +
        ggplot2::geom_path(data = g$edges, ggplot2::aes(x = .data$x, y = .data$y),
            colour = 'grey40', linewidth = 0.6) +
        ggplot2::geom_point(data = df,
            ggplot2::aes(x = .data$x, y = .data$y, fill = .data$call, size = .data$median_mde),
            shape = 21, colour = 'white', stroke = 1.1) +
        ggplot2::scale_size_continuous(range = c(4, 11), name = 'median MDE',
            guide = ggplot2::guide_legend(override.aes = list(shape = 21, fill = 'grey55', colour = 'white'))) +
        ggrepel::geom_text_repel(data = df,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$cluster),
            size = 4, box.padding = 1, point.padding = 0.5) +
        ggplot2::geom_text(data = g$corners,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab, vjust = .data$vj),
            fontface = 'bold', size = 4.5) +
        ggplot2::scale_fill_manual(values = seny_call_colors, drop = FALSE, name = 'call',
            guide = ggplot2::guide_legend(override.aes = list(size = 5))) +
        ggplot2::coord_equal(clip = 'off') +
        ggplot2::theme_void(base_size = 13) +
        ggplot2::theme(plot.margin = ggplot2::margin(26, 38, 26, 38))
}
