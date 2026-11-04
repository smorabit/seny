skip_report <- function() {
    skip_if_not_installed('rmarkdown')
    skip_if_not_installed('DT')
    skip_if_not(rmarkdown::pandoc_available(), 'pandoc not available')
}

report_res <- function() {
    sce <- readRDS(system.file('extdata', 'example_sce.rds', package = 'seny'))
    assessAll(sce, sample = 'donor', group = 'condition',
              annotations = c('cell_type', 'cell_type_fine'), targetLFC = 1, thetaReps = 0)
}

test_that('the report renders every comparison, the warning and capped cards', {
    skip_report()
    res <- report_res()
    out <- reportSeny(res, file = tempfile(fileext = '.html'), maxCards = 2, open = FALSE)
    expect_true(file.exists(out))
    html <- paste(readLines(out, warn = FALSE), collapse = '\n')
    expect_match(html, sprintf('%d comparisons were assessed', nrow(res)))
    for (cur in unique(as.character(res$cluster))) expect_match(html, cur, fixed = TRUE)
    expect_equal(lengths(regmatches(html, gregexpr('id="card-', html))), 2)
    expect_match(html, sprintf('%d comparison\\(s\\) left out', nrow(res) - 2))
})

test_that('the report needs an assessAll() result', {
    skip_report()
    res <- report_res()
    bare <- res
    attr(bare, 'pb') <- NULL
    expect_error(reportSeny(bare, file = tempfile(fileext = '.html'), open = FALSE), 'assessAll')
})

test_that('settings and flags are attached for the report', {
    res <- report_res()
    s <- attr(res, 'settings')
    expect_equal(s$method, 'edgeR')
    expect_equal(s$targetLFC, 1)
    expect_true(all(c('cramerV', 'compFrac') %in% names(s$thresholds)))
    # under edgeR an imbalance is reported without changing an Adequate call
    expect_true(all(res$flags[res$call == 'Adequate'] %in% c('', 'imbalance')))
    expect_true(all(grepl('composition', res$flags[grepl('composition shift', res$reason)])))
})

test_that('the template can be drafted for editing', {
    skip_if_not_installed('rmarkdown')
    dir <- tempfile()
    dir.create(dir)
    f <- file.path(dir, 'my_report.Rmd')
    rmarkdown::draft(f, template = 'seny_report', package = 'seny', edit = FALSE)
    expect_true(file.exists(f))
})
