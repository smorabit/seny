test_that('the plotX family returns ggplots on an assessAll result', {
    sce <- sim_sce(n_per_group = 10, n_genes = 150, cells_per = 50, seed = 21)
    res <- assessAll(sce, sample = 'donor', group = 'condition', annotations = 'cell_type')
    expect_true(inherits(plotConfounding(res), 'ggplot'))
    expect_true(inherits(plotCells(res), 'ggplot'))
    expect_true(inherits(plotVariance(res), 'ggplot'))
    expect_true(inherits(plotRegime(res), 'ggplot'))
    expect_true(inherits(plotSensitivity(res), 'ggplot'))
    expect_true(inherits(plotPower(res), 'ggplot'))
    cmp <- as.data.frame(res)$comparison[1]
    expect_true(inherits(plotComparison(res, cmp), 'patchwork'))
})

test_that('per-gene data is attached for the sensitivity plots', {
    sce <- sim_sce(n_per_group = 8, n_genes = 120, cells_per = 40, seed = 20)
    res <- assessAll(sce, sample = 'donor', group = 'condition', annotations = 'cell_type')
    pg <- attr(res, 'pergene')
    expect_true(all(c('comparison', 'gene', 'acpm', 'se_log2', 'df') %in% colnames(pg)))
    expect_true(all(is.finite(pg$se_log2)))
})

test_that('attr-dependent plots error on a bare callComparisons result', {
    sce <- sim_sce(n_per_group = 8, n_genes = 120, cells_per = 40, seed = 22)
    x <- pbSetup(sce, sample = 'donor', group = 'condition', annotations = 'cell_type')
    v <- estimateVariance(x)
    bare <- callComparisons(assessDesign(x), assessPseudobulk(x), v$summary, estimateSensitivity(v))
    cmp <- as.data.frame(bare)$comparison[1]
    expect_error(plotCells(bare), 'assessAll')
    expect_error(plotSensitivity(bare), 'assessAll')
    expect_error(plotComparison(bare, cmp), 'assessAll')
})

test_that('plotComparison errors on an unknown comparison', {
    sce <- sim_sce(n_per_group = 8, n_genes = 120, cells_per = 40, seed = 23)
    res <- assessAll(sce, sample = 'donor', group = 'condition', annotations = 'cell_type')
    expect_error(plotComparison(res, 'not-a-comparison'), 'unknown comparison')
})

test_that('NB-only plots have nothing to show under the voom method', {
    sce <- sim_sce(n_per_group = 10, n_genes = 150, cells_per = 50, seed = 24)
    res <- assessAll(sce, sample = 'donor', group = 'condition', annotations = 'cell_type',
                     method = 'voom')
    expect_true(all(is.na(as.data.frame(res)$share_donor)))
    expect_error(suppressMessages(plotVariance(res)), 'no variance decomposition')
    expect_error(suppressWarnings(plotRegime(res)), 'no variance decomposition')
    # the method-agnostic plots still work under voom
    expect_true(inherits(plotSensitivity(res), 'ggplot'))
    expect_true(inherits(plotPower(res), 'ggplot'))
})
