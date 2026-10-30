#---------------------------------------------------------#
# Shared plotting internals: palettes, theme, ternary geometry
#---------------------------------------------------------#

# traffic-light fill scale (established CVD-safe RdYlGn status order)
seny_call_colors <- c(Adequate = '#1A9850', Underpowered = '#FEE08B',
                      Caution = '#FDAE61', Invalid = '#D73027')

# variance-term fill (Seaborn-deep blue/orange/green)
seny_share_colors <- c(depth = '#4C72B0', cell = '#DD8452', donor = '#55A868')

# two-condition palette (Okabe-Ito)
seny_cond_colors <- c('#0072B2', '#E69F00')

# the house theme every seny plot shares
#' @noRd
seny_theme <- function(base_size = 12) {
    ggplot2::theme_minimal(base_size = base_size) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(linewidth = 0.8, colour = 'grey45', fill = NA),
            plot.title = ggplot2::element_text(face = 'bold'),
            strip.text = ggplot2::element_text(face = 'bold'))
}

# fetch an attribute assessAll() attaches, with a clear message when it is absent
#' @noRd
require_attr <- function(res, which) {
    a <- attr(res, which)
    if (is.null(a)) {
        stop('this plot needs `attr(res, "', which, '")` from assessAll(); pass an ',
             'assessAll() result, not a bare callComparisons() table.', call. = FALSE)
    }
    a
}

# per-gene minimum detectable effect from the stored SE and df (shared by the
# sensitivity plot and the card). The MDE does not depend on the target effect.
#' @noRd
mde_from_se <- function(se_log2, df, power = 0.8, fdr = 0.05, pi1 = 0.05, alphaFixed = 0.001) {
    alpha <- sensitivity_alpha(power, fdr, pi1, alphaFixed)[['call']]
    (stats::qt(1 - alpha / 2, df) + stats::qt(power, df)) * se_log2
}

# add a per-annotation facet only when there is more than one annotation set
#' @noRd
facet_by_annotation <- function(p, df) {
    if (length(unique(df$annotation)) > 1) {
        p + ggplot2::facet_wrap(~ annotation, scales = 'free_y')
    } else {
        p
    }
}

#' @noRd
power_from_se <- function(se_log2, df, target, power = 0.8, fdr = 0.05, pi1 = 0.05,
                          alphaFixed = 0.001) {
    alpha <- sensitivity_alpha(power, fdr, pi1, alphaFixed)[['call']]
    ncp <- target / se_log2
    tc <- stats::qt(1 - alpha / 2, df)
    stats::pt(tc, df, ncp, lower.tail = FALSE) + stats::pt(-tc, df, ncp, lower.tail = TRUE)
}

#---------------------------------------------------------#
# Ternary (3-simplex) geometry, hand-rolled so no ternary package is needed.
# corners: depth (0,0), cell (1,0), donor (0.5, sqrt(3)/2).
#---------------------------------------------------------#

# map depth/cell/donor shares to x,y in the triangle
#' @noRd
ternary_xy <- function(depth, cell, donor) {
    s <- depth + cell + donor
    data.frame(x = cell / s + 0.5 * donor / s, y = (sqrt(3) / 2) * (donor / s))
}

# the triangle outline, interior gridlines and edge ticks at the given levels
#' @noRd
ternary_geom <- function(levels = c(0.2, 0.4, 0.6, 0.8)) {
    sq <- sqrt(3) / 2
    L <- levels
    grid <- rbind(
        data.frame(x = L / 2,       y = sq * L, xend = 1 - L / 2,     yend = sq * L),         # const donor
        data.frame(x = L,           y = 0,      xend = 0.5 + L / 2,   yend = sq * (1 - L)),   # const cell
        data.frame(x = 1 - L,       y = 0,      xend = 0.5 * (1 - L), yend = sq * (1 - L)))   # const depth
    ticks <- rbind(
        data.frame(x = L,             y = 0,            lab = L * 100, hj = 0.5,  vj = 1.9),   # cell (bottom)
        data.frame(x = 1 - L / 2,     y = sq * L,       lab = L * 100, hj = -0.5, vj = 0.3),   # donor (right)
        data.frame(x = 0.5 * (1 - L), y = sq * (1 - L), lab = L * 100, hj = 1.5,  vj = 0.3))   # depth (left)
    corners <- data.frame(x = c(0, 1, 0.5), y = c(0, 0, sq),
                          lab = c('depth', 'cell', 'donor'), vj = c(1.9, 1.9, -0.9))
    edges <- corners[c(1, 2, 3, 1), c('x', 'y')]
    list(grid = grid, ticks = ticks, corners = corners, edges = edges)
}
