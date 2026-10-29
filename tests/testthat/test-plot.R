test_that('plotAdequacy returns one ggplot per annotation set', {
    sce <- sim_sce(n_per_group = 8, n_genes = 150, cells_per = 50, seed = 9)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type')
    plots <- plotAdequacy(res)
    expect_type(plots, 'list')
    expect_named(plots, 'cell_type')
    expect_true(inherits(plots[['cell_type']], 'ggplot'))
})

test_that('plotAdequacy covers every comparison as a tile', {
    sce <- sim_sce(n_per_group = 8, n_genes = 150, cells_per = 50, seed = 9)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type')
    p <- plotAdequacy(res)[['cell_type']]
    built <- ggplot2::ggplot_build(p)
    expect_equal(nrow(built$data[[1]]), nrow(as.data.frame(res)))
})
