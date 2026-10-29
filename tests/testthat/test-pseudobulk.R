test_that('gini is 0 for equal values and near 1 for one-dominant', {
    expect_equal(gini(rep(5, 10)), 0)
    expect_true(gini(c(rep(0, 99), 100)) > 0.9)
    expect_equal(gini(numeric(0)), 0)
    expect_equal(gini(7), 0)
})

test_that('point_biserial sign tracks the group with more cells, NA without variation', {
    cond <- rep(0:1, each = 5)
    more_in_b <- c(rep(10, 5), rep(40, 5))
    more_in_a <- c(rep(40, 5), rep(10, 5))
    expect_true(point_biserial(log(more_in_b), cond) > 0)
    expect_true(point_biserial(log(more_in_a), cond) < 0)
    expect_true(is.na(point_biserial(rep(3, 10), cond)))
})

test_that('assessPseudobulk reports the full layer-2 columns in range', {
    sce <- sim_sce(n_per_group = 8, n_genes = 150, cells_per = 50, seed = 8)
    x <- pbSetup(sce, sample = 'donor', group = 'condition', annotations = 'cell_type')
    pb <- assessPseudobulk(x)
    expect_true(all(c('gini_cells', 'frac_below_min', 'cells_ratio', 'cells_cond_cor',
                      'libsize_cond_cor', 'n_testable', 'frac_testable') %in% colnames(pb)))
    expect_true(all(pb$frac_testable >= 0 & pb$frac_testable <= 1))
    expect_true(all(pb$cells_ratio >= 0 & pb$cells_ratio <= 1))
    expect_true(all(pb$frac_below_min >= 0 & pb$frac_below_min <= 1))
    expect_true(all(pb$gini_cells >= 0 & pb$gini_cells <= 1))
})

test_that('frac_below_min catches samples dropped below minCells', {
    # thin one condition-B donor down to 3 cells so it is dropped at minCells = 10
    sce <- sim_sce(n_per_group = 6, n_genes = 120, cells_per = 40, seed = 10)
    d7 <- which(sce$donor == 'd7')                    # a condition-B donor
    sce2 <- sce[, -d7[-(1:3)]]
    x <- pbSetup(sce2, sample = 'donor', group = 'condition',
                 annotations = 'cell_type', minCells = 10)
    pb <- assessPseudobulk(x)
    expect_true(any(pb$frac_below_min > 0))
})
