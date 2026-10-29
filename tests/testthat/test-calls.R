test_that('a well-powered clean design is Adequate', {
    sce <- sim_sce(n_per_group = 20, phi_donor = 0.01, theta = 0.2,
                   n_genes = 200, cells_per = 80, seed = 3)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', targetLFC = 1)
    expect_true(all(as.character(res$call) == 'Adequate'))
})

test_that('few donors with high donor variance is Underpowered', {
    sce <- sim_sce(n_per_group = 3, phi_donor = 0.25, theta = 0.5,
                   n_genes = 200, cells_per = 30, seed = 4)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', targetLFC = 0.5)
    expect_true(any(as.character(res$call) == 'Underpowered'))
    expect_false(any(as.character(res$call) == 'Adequate'))
})

test_that('a batch nested in condition is Invalid', {
    sce <- sim_sce(n_per_group = 10, n_genes = 200, cells_per = 50,
                   batch = 'nested', seed = 5)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', covariates = 'batch')
    expect_true(all(as.character(res$call) == 'Invalid'))
})

test_that('a partially confounded batch is Caution', {
    sce <- sim_sce(n_per_group = 15, phi_donor = 0.01, n_genes = 200,
                   cells_per = 60, batch = 'partial', seed = 6)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', covariates = 'batch', targetLFC = 1)
    expect_true(all(as.character(res$call) == 'Caution'))
})

test_that('cell-count imbalance correlated with condition is Caution', {
    sce <- sim_sce(n_per_group = 15, phi_donor = 0.01, theta = 0.2, n_genes = 200,
                   cells_per = 80, cells_imbalance = 0.2, seed = 7)
    res <- assessAll(sce, sample = 'donor', group = 'condition',
                     annotations = 'cell_type', targetLFC = 1)
    expect_true(all(as.character(res$call) == 'Caution'))
    expect_true(all(grepl('imbalance', res$reason)))
})

test_that('too few testable genes is Underpowered', {
    cmp <- 'a:c:B_vs_A'
    design <- data.frame(comparison = cmp, estimable = TRUE, nA = 10, nB = 10,
                         max_cramerV = 0, max_vif = 1, stringsAsFactors = FALSE)
    pseudobulk <- data.frame(comparison = cmp, cells_cond_cor = 0, n_testable = 50,
                             stringsAsFactors = FALSE)
    variance <- data.frame(comparison = cmp, stringsAsFactors = FALSE)
    sensitivity <- data.frame(comparison = cmp, sensitivity_fraction = 0.9,
                              target_lfc = 1, stringsAsFactors = FALSE)
    out <- callComparisons(design, pseudobulk, variance, sensitivity)
    expect_equal(as.character(out$call), 'Underpowered')
    expect_true(grepl('testable', out$reason))
})
