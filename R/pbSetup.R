#' Set up pseudobulk comparisons from a single-cell dataset
#'
#' Aggregates a single-cell dataset to pseudobulk once per annotation level and
#' returns the object the rest of the package operates on: a named
#' [S4Vectors::SimpleList] holding one pseudobulk [SummarizedExperiment] per
#' annotation set. Each pseudobulk is a gene-by-(sample x cluster) matrix with
#' per-pseudobulk metadata in `colData`, and the comparison definitions, column
#' names and per-cluster cell-level overdispersion (`theta`) in `metadata()`.
#'
#' Cell-level `theta` (within-donor cell-to-cell overdispersion, needed by
#' [estimateVariance()]) is estimated here, while the cells are in hand, so the
#' later steps never need the single-cell object again.
#'
#' @param sce A [SingleCellExperiment::SingleCellExperiment] or
#'   [SummarizedExperiment] with a raw-count assay (non-negative integers).
#' @param sample Character scalar; `colData` column giving the biological sample
#'   (donor) each cell belongs to.
#' @param group Character scalar; `colData` column giving the condition compared.
#' @param annotations Character vector of `colData` columns, each a clustering
#'   (annotation level) to pseudobulk, e.g. `c("cell_type_broad")`.
#' @param covariates Character vector of `colData` columns to carry as candidate
#'   confounders. Default none.
#' @param contrasts Optional list of length-2 character vectors `c(A, B)` naming
#'   `group` levels to compare (B vs A). Default `NULL`: all non-reference levels
#'   versus the first level.
#' @param pairing Optional `colData` column for paired designs (reserved; not
#'   used in the MVP).
#' @param minCells Integer; pseudobulks with fewer contributing cells are
#'   dropped. Default 10.
#' @param aggregation Either `"sum"` (counts, for the negative-binomial model) or
#'   `"mean"` (reserved for the limma path). Default `"sum"`.
#' @param assayName Name of the raw-count assay in `sce`. Default `"counts"`.
#' @param thetaCells,thetaDonors,thetaReps Subsample size (cells per donor,
#'   donors) and number of draws for the cell-level `theta` ensemble.
#' @param BPPARAM A `BiocParallelParam` object for the per-cluster theta
#'   estimation. Default [BiocParallel::SerialParam()].
#'
#' @return A named [S4Vectors::SimpleList] of [SummarizedExperiment], one per
#'   annotation level.
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames colData
#' @importFrom S4Vectors SimpleList metadata metadata<-
pbSetup <- function(sce, sample, group, annotations, covariates = character(0),
                    contrasts = NULL, pairing = NULL, minCells = 10,
                    aggregation = c('sum', 'mean'), assayName = 'counts',
                    thetaCells = 50, thetaDonors = 50, thetaReps = 3,
                    BPPARAM = BiocParallel::SerialParam()) {
    aggregation <- match.arg(aggregation)
    cd <- as.data.frame(SummarizedExperiment::colData(sce))
    need <- c(sample, group, annotations, covariates)
    miss <- setdiff(need, colnames(cd))
    if (length(miss)) stop('columns not found in colData: ', paste(miss, collapse = ', '))

    if (!assayName %in% SummarizedExperiment::assayNames(sce)) {
        assayName <- SummarizedExperiment::assayNames(sce)[1]
    }
    counts <- SummarizedExperiment::assay(sce, assayName)
    nzv <- if (methods::is(counts, 'sparseMatrix')) counts@x else counts[counts != 0]
    if (any(nzv < 0)) stop('counts contain negative values')
    if (any(nzv != round(nzv))) stop('counts are not integers; aggregate raw counts, never round')

    group_vec <- as.character(cd[[group]])
    contrasts <- resolve_contrasts(group_vec, contrasts)

    sets <- BiocParallel::bplapply(annotations, function(ann) {
        build_annotation_set(counts, cd, sample, group, ann, covariates, contrasts,
                             minCells, aggregation, thetaCells, thetaDonors, thetaReps)
    }, BPPARAM = BPPARAM)
    names(sets) <- annotations
    S4Vectors::SimpleList(sets)
}

#' Default contrasts: every non-reference level versus the first level
#' @noRd
resolve_contrasts <- function(group_vec, contrasts) {
    lv <- levels(factor(group_vec))
    if (!is.null(contrasts)) return(contrasts)
    if (length(lv) < 2) stop('group has fewer than two levels; nothing to compare')
    lapply(lv[-1], function(b) c(lv[1], b))
}

#' Aggregate one annotation level and attach per-cluster theta
#' @noRd
build_annotation_set <- function(counts, cd, sample, group, ann, covariates,
                                 contrasts, minCells, aggregation, thetaCells,
                                 thetaDonors, thetaReps) {
    cluster <- as.character(cd[[ann]])
    keep_cell <- !is.na(cluster)
    counts <- counts[, keep_cell]
    cd <- cd[keep_cell, , drop = FALSE]
    cluster <- cluster[keep_cell]

    pbid <- interaction(cd[[sample]], cluster, drop = TRUE, sep = '__')
    pb <- pb_sum(counts, pbid)
    ncells <- as.integer(table(pbid)[colnames(pb)])
    # mean aggregation (non-integer) pairs with the limma path; sum stays integer
    if (aggregation == 'mean') pb <- sweep(pb, 2, ncells, '/')

    # per-pseudobulk metadata: sample, cluster, group, covariates (constant within pb)
    first_of <- function(col) unname(vapply(split(cd[[col]], pbid),
        function(x) as.character(x[1]), character(1))[colnames(pb)])
    meta <- data.frame(sample = first_of(sample), cluster = first_of(ann),
                       group = first_of(group), nCells = ncells,
                       row.names = colnames(pb), stringsAsFactors = FALSE)
    for (cov in covariates) meta[[cov]] <- first_of(cov)

    # pre-filter per-sample cell counts, so layer 2 can report the fraction of
    # samples dropped for falling below minCells (the drop happens next)
    cell_counts <- meta[, c('cluster', 'group', 'nCells')]

    keep_pb <- meta$nCells >= minCells
    pb <- pb[, keep_pb, drop = FALSE]
    meta <- meta[keep_pb, , drop = FALSE]

    # per-cluster cell-level theta, from the cells of this annotation level
    thetas <- lapply(unique(meta$cluster), function(cl) {
        sel <- cluster == cl
        estimate_cell_theta(counts[, sel], cd[[sample]][sel], n_rep = thetaReps,
                            max_cells = thetaCells, max_donors = thetaDonors)
    })
    names(thetas) <- unique(meta$cluster)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = pb),
        colData = S4Vectors::DataFrame(meta))
    S4Vectors::metadata(se) <- list(sample = sample, group = group,
        covariates = covariates, contrasts = contrasts, minCells = minCells,
        aggregation = aggregation, theta = thetas, annotation = ann,
        cell_counts = cell_counts)
    se
}
