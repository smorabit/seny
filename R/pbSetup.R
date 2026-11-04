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
#'   [SummarizedExperiment] with a raw-count assay (non-negative integers). The
#'   assay may be a base or sparse matrix, or an out-of-memory matrix such as a
#'   `DelayedMatrix` or an HDF5-backed matrix (for atlases beyond one sparse
#'   matrix's 2^31 nonzeros): it is then checked and summed a block of cells at
#'   a time (`options(seny.block_cells = ...)`, default 1e5).
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
#'   donors) and number of draws for the cell-level `theta` ensemble. The
#'   default of 20 donors per draw gives the same variance shares and regime
#'   as 50, about 5 times faster (each draw's Cox-Reid-adjusted fit grows
#'   steeply with the number of donors).
#'   The default `thetaReps = 0` skips the theta fit, the slowest step (about
#'   45 s per cluster with 3 reps): the adequacy call is the same either way,
#'   since power uses the total pseudobulk dispersion, but the depth/cell/donor
#'   variance decomposition from [estimateVariance()] (shares, regime,
#'   [plotVariance()], [plotRegime()]) needs theta. Use `thetaReps = 3` to get
#'   it, i.e. to see whether more donors or more cells would help.
#' @param BPPARAM A `BiocParallelParam` object for the per-cluster theta
#'   estimation. Default [BiocParallel::SerialParam()].
#'
#' @return A named [S4Vectors::SimpleList] of [SummarizedExperiment], one per
#'   annotation level.
#' @examples
#' sce <- readRDS(system.file("extdata", "example_sce.rds", package = "seny"))
#' x <- pbSetup(sce, sample = "donor", group = "condition",
#'              annotations = "cell_type", thetaReps = 0)
#' x
#' names(x)
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames colData
#' @importFrom S4Vectors SimpleList metadata metadata<-
pbSetup <- function(sce, sample, group, annotations, covariates = character(0),
                    contrasts = NULL, pairing = NULL, minCells = 10,
                    aggregation = c('sum', 'mean'), assayName = 'counts',
                    thetaCells = 50, thetaDonors = 20, thetaReps = 0,
                    BPPARAM = BiocParallel::SerialParam()) {
    aggregation <- match.arg(aggregation)
    if (length(thetaReps) != 1 || is.na(thetaReps) || thetaReps < 0 || thetaReps != round(thetaReps)) {
        stop('thetaReps must be a non-negative integer (0 skips the theta fit)')
    }
    cd <- as.data.frame(SummarizedExperiment::colData(sce))
    need <- c(sample, group, annotations, covariates)
    miss <- setdiff(need, colnames(cd))
    if (length(miss)) stop('columns not found in colData: ', paste(miss, collapse = ', '))

    if (!assayName %in% SummarizedExperiment::assayNames(sce)) {
        assayName <- SummarizedExperiment::assayNames(sce)[1]
    }
    counts <- SummarizedExperiment::assay(sce, assayName)
    check_counts(counts)

    group_vec <- as.character(cd[[group]])
    contrasts <- resolve_contrasts(group_vec, contrasts)

    sets <- BiocParallel::bplapply(annotations, function(ann) {
        build_annotation_set(counts, cd, sample, group, ann, covariates, contrasts,
                             minCells, aggregation, thetaCells, thetaDonors, thetaReps)
    }, BPPARAM = BPPARAM)
    names(sets) <- annotations
    out <- S4Vectors::SimpleList(sets)
    if (length(annotations) >= 2) {
        S4Vectors::metadata(out)$hierarchy <- infer_hierarchy(cd, sample, annotations)
    }
    out
}

#' Parent-child map between annotation levels, from the cells
#'
#' For every pair of annotation levels, the coarser one (fewer clusters) is the
#' parent when its clusters nest the finer ones: each child goes to the parent
#' holding most of its cells, and the pair is kept when at least 80% of the
#' cells sit in their child's majority parent. Returns the child-to-parent map
#' (with each child's purity) and the per-sample cell counts of every
#' parent-child cell set, which give the children's shares within a parent.
#' @noRd
infer_hierarchy <- function(cd, sample, annotations, min_purity = 0.8) {
    n_cl <- vapply(annotations, function(a) length(unique(stats::na.omit(cd[[a]]))), numeric(1))
    maps <- list(); cells <- list()
    for (pa in annotations) for (ch in annotations) {
        if (n_cl[[pa]] >= n_cl[[ch]]) next
        ok <- !is.na(cd[[pa]]) & !is.na(cd[[ch]])
        tab <- table(as.character(cd[[ch]][ok]), as.character(cd[[pa]][ok]))
        if (sum(apply(tab, 1, max)) / sum(tab) < min_purity) next
        maps[[length(maps) + 1]] <- data.frame(
            parent_annotation = pa, child_annotation = ch, child = rownames(tab),
            parent = colnames(tab)[max.col(tab, ties.method = 'first')],
            purity = apply(tab, 1, max) / rowSums(tab), n_cells = as.numeric(rowSums(tab)),
            stringsAsFactors = FALSE)
        jt <- as.data.frame(table(sample = as.character(cd[[sample]][ok]),
                                  parent = as.character(cd[[pa]][ok]),
                                  child = as.character(cd[[ch]][ok])),
                            stringsAsFactors = FALSE)
        jt <- jt[jt$Freq > 0, ]
        cells[[length(cells) + 1]] <- data.frame(parent_annotation = pa, child_annotation = ch,
            jt[, c('sample', 'parent', 'child')], nCells = jt$Freq, stringsAsFactors = FALSE)
    }
    if (length(maps) == 0) return(NULL)
    list(map = do.call(rbind, maps), cells = do.call(rbind, cells))
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
    # per-cluster gene totals over all its cells, before minCells: the per-cell
    # expression of each cluster that the composition layer mixes
    pb_cluster <- unname(vapply(split(cluster, pbid), function(x) x[1], character(1))[colnames(pb)])
    totals <- cluster_totals(pb, ncells, pb_cluster, cd, pbid, group)
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
    cell_counts <- meta[, c('sample', 'cluster', 'group', 'nCells')]

    keep_pb <- meta$nCells >= minCells
    pb <- pb[, keep_pb, drop = FALSE]
    meta <- meta[keep_pb, , drop = FALSE]

    # per-cluster cell-level theta, from the cells of this annotation level.
    # thetaReps = 0 skips it: the call does not use theta, only the decomposition
    thetas <- NULL
    if (thetaReps > 0) {
        thetas <- lapply(unique(meta$cluster), function(cl) {
            sel <- cluster == cl
            estimate_cell_theta(counts[, sel], cd[[sample]][sel], n_rep = thetaReps,
                                max_cells = thetaCells, max_donors = thetaDonors)
        })
        names(thetas) <- unique(meta$cluster)
    }

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = pb),
        colData = S4Vectors::DataFrame(meta))
    S4Vectors::metadata(se) <- list(sample = sample, group = group,
        covariates = covariates, contrasts = contrasts, minCells = minCells,
        aggregation = aggregation, theta = thetas, annotation = ann,
        cell_counts = cell_counts, cluster_totals = totals)
    se
}

#' Gene totals and cell numbers per cluster, overall and per group, from the
#' (sum) pseudobulks before any minCells filter
#' @noRd
cluster_totals <- function(pb, ncells, pb_cluster, cd, pbid, group) {
    pb_group <- unname(vapply(split(as.character(cd[[group]]), pbid),
                              function(x) x[1], character(1))[colnames(pb)])
    sum_by <- function(sel) {
        f <- factor(pb_cluster[sel])
        ind <- Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
                                    dims = c(length(f), nlevels(f)),
                                    dimnames = list(NULL, levels(f)))
        list(counts = as.matrix(pb[, sel, drop = FALSE] %*% ind),
             cells = stats::setNames(as.numeric(tapply(ncells[sel], f, sum)), levels(f)))
    }
    by_group <- lapply(stats::setNames(nm = unique(pb_group)), function(g) sum_by(pb_group == g))
    c(sum_by(rep(TRUE, length(pb_cluster))), list(by_group = by_group))
}
