test_that('assessAll runs end to end and returns a tidy comparison table', {
    sce <- sim_sce(n_per_group = 10, n_genes = 200, cells_per = 50, seed = 1)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type')

    expect_s4_class(res, 'DataFrame')
    # 2 cell types x 1 contrast (B vs A)
    expect_equal(nrow(res), 2L)
    expect_true(all(c('comparison', 'cluster', 'call', 'reason', 'median_phi_donor',
                      'sensitivity_fraction', 'expected_recall', 'regime') %in% colnames(res)))
    expect_true(is.factor(res$call))
    expect_true(all(as.character(res$call) %in%
                    c('Adequate', 'Underpowered', 'Caution', 'Invalid')))
})

test_that('variance shares sum to one and phi_donor recovers the planted value', {
    sce <- sim_sce(n_per_group = 12, phi_donor = 0.05, n_genes = 200,
                   cells_per = 60, seed = 2)
    x <- pbSetup(sce, sample = 'donor', group = 'condition', annotations = 'cell_type')
    v <- estimateVariance(x)
    s <- v$summary
    expect_true(all(abs(s$share_depth + s$share_cell + s$share_donor - 1) < 1e-8))
    # planted phi_donor 0.05; moment estimate should land in a sane range
    expect_true(all(s$median_phi_donor > 0.005 & s$median_phi_donor < 0.2))
})
