#' Traffic-light matrix of adequacy calls
#'
#' Renders the comparison table from [callComparisons()] (or [assessAll()]) as a
#' clusters-by-contrasts grid, one plot per annotation set, each tile coloured by
#' the comparison's call (Adequate / Underpowered / Caution / Invalid). The call
#' reasons stay in the table; the matrix is the at-a-glance view.
#'
#' @param res The [S4Vectors::DataFrame] returned by [callComparisons()] or
#'   [assessAll()], one row per comparison, with `annotation`, `cluster`, `A`,
#'   `B` and `call`.
#' @return A named list of `ggplot` objects, one per annotation set (a
#'   single-level input returns a one-element list).
#' @export
#' @importFrom ggplot2 .data
plotAdequacy <- function(res) {
    df <- as.data.frame(res)
    df$contrast <- paste0(df$B, '_vs_', df$A)
    df$call <- factor(as.character(df$call), levels = names(seny_call_colors))

    sets <- split(df, df$annotation)
    out <- lapply(names(sets), function(ann) plot_one_matrix(sets[[ann]], ann))
    names(out) <- names(sets)
    out
}

#' One annotation set's tile grid
#' @noRd
plot_one_matrix <- function(d, ann) {
    d$cluster <- factor(d$cluster, levels = rev(sort(unique(d$cluster))))
    d$contrast <- factor(d$contrast, levels = sort(unique(d$contrast)))
    ggplot2::ggplot(d, ggplot2::aes(x = .data$contrast, y = .data$cluster, fill = .data$call)) +
        ggplot2::geom_tile(colour = 'white', linewidth = 0.6) +
        ggplot2::scale_fill_manual(values = seny_call_colors,
            limits = names(seny_call_colors), drop = FALSE,
            na.value = 'grey90', name = 'Call') +
        ggplot2::coord_equal() +
        ggplot2::labs(x = 'contrast', y = 'cluster', title = ann) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(linewidth = 1, colour = 'black', fill = NA),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            plot.title = ggplot2::element_text(hjust = 0.5))
}
