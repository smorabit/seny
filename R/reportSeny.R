#' Shareable HTML report of an assessment
#'
#' Renders an [assessAll()] result into one self-contained HTML file: a summary
#' of the calls and the reasons behind them, a reminder that every comparison
#' run should be reported, the traffic-light matrices, a sortable and
#' searchable table of every comparison, the design and variance plots, the
#' composition flags, a card per comparison ([plotComparison()]), and the
#' thresholds used with their sources.
#'
#' Cards go to flagged comparisons first (Invalid, Caution, Underpowered, then
#' Adequate; by expected recall within each), up to `maxCards`; the table
#' always lists every comparison. The report is an R Markdown template: for an
#' editable copy use
#' `rmarkdown::draft("my_report.Rmd", template = "seny_report", package = "seny")`
#' and set its `res_path` parameter to a saved [assessAll()] result.
#'
#' Needs the rmarkdown and DT packages and pandoc (bundled with RStudio;
#' elsewhere install it, e.g. from conda-forge).
#'
#' @param res The result of [assessAll()].
#' @param file Path of the HTML file to write. Default `"seny_report.html"`.
#' @param title Report title.
#' @param maxCards Maximum number of comparison cards. Default 100.
#' @param open Open the report in a browser when done. Default
#'   `interactive()`.
#' @return The path of the report, invisibly.
#' @examples
#' if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available() &&
#'     requireNamespace("DT", quietly = TRUE)) {
#'     sce <- readRDS(system.file("extdata", "example_sce.rds", package = "seny"))
#'     res <- assessAll(sce, sample = "donor", group = "condition",
#'                      annotations = "cell_type", targetLFC = 1, thetaReps = 0)
#'     out <- reportSeny(res, file = tempfile(fileext = ".html"), open = FALSE)
#'     file.exists(out)
#' }
#' @export
reportSeny <- function(res, file = 'seny_report.html', title = 'seny report', maxCards = 100,
                       open = interactive()) {
    require_attr(res, 'pb')
    for (pkg in c('rmarkdown', 'DT')) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop('reportSeny() needs the "', pkg, '" package; install it first', call. = FALSE)
        }
    }
    if (!rmarkdown::pandoc_available()) {
        stop('reportSeny() needs pandoc (bundled with RStudio; otherwise install it, e.g. ',
             'from conda-forge)', call. = FALSE)
    }
    template <- system.file('rmarkdown', 'templates', 'seny_report', 'skeleton', 'skeleton.Rmd',
                            package = 'seny')
    work <- tempfile('seny_report_')
    dir.create(work)
    on.exit(unlink(work, recursive = TRUE), add = TRUE)
    rmd <- file.path(work, 'seny_report.Rmd')
    file.copy(template, rmd)
    res_path <- file.path(work, 'result.rds')
    saveRDS(res, res_path)
    out <- normalizePath(file, mustWork = FALSE)
    rmarkdown::render(rmd, output_file = basename(out), output_dir = dirname(out),
                      intermediates_dir = work, knit_root_dir = work,
                      params = list(res_path = res_path, max_cards = maxCards, title = title),
                      envir = new.env(), quiet = TRUE)
    if (open) utils::browseURL(out)
    invisible(out)
}
