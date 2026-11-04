# Shareable HTML report of an assessment

Renders an
[`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
result into one self-contained HTML file: a summary of the calls and the
reasons behind them, a reminder that every comparison run should be
reported, the traffic-light matrices, a sortable and searchable table of
every comparison, the design and variance plots, the composition flags,
a card per comparison
([`plotComparison()`](https://smorabit.github.io/seny/reference/plotComparison.md)),
and the thresholds used with their sources.

## Usage

``` r
reportSeny(
  res,
  file = "seny_report.html",
  title = "seny report",
  maxCards = 100,
  open = interactive()
)
```

## Arguments

- res:

  The result of
  [`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md).

- file:

  Path of the HTML file to write. Default `"seny_report.html"`.

- title:

  Report title.

- maxCards:

  Maximum number of comparison cards. Default 100.

- open:

  Open the report in a browser when done. Default
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

## Value

The path of the report, invisibly.

## Details

Cards go to flagged comparisons first (Invalid, Caution, Underpowered,
then Adequate; by expected recall within each), up to `maxCards`; the
table always lists every comparison. The report is an R Markdown
template: for an editable copy use
`rmarkdown::draft("my_report.Rmd", template = "seny_report", package = "seny")`
and set its `res_path` parameter to a saved
[`assessAll()`](https://smorabit.github.io/seny/reference/assessAll.md)
result.

Needs the rmarkdown and DT packages and pandoc (bundled with RStudio;
elsewhere install it, e.g. from conda-forge).

## Examples

``` r
if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available() &&
    requireNamespace("DT", quietly = TRUE)) {
    sce <- readRDS(system.file("extdata", "example_sce.rds", package = "seny"))
    res <- assessAll(sce, sample = "donor", group = "condition",
                     annotations = "cell_type", targetLFC = 1, thetaReps = 0)
    out <- reportSeny(res, file = tempfile(fileext = ".html"), open = FALSE)
    file.exists(out)
}
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> [1] TRUE
```
