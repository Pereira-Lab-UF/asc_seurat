#' Detect Doublets with scDblFinder
#'
#' @param obj A Seurat object.
#' @param sample_col Optional metadata column identifying separate samples.
#' @return The Seurat object with scDblFinder metadata columns.
#' @keywords internal
detect_doublets <- function(obj, sample_col = NULL) {
    if (!requireNamespace("scDblFinder", quietly = TRUE)) {
        stop(
            "Package 'scDblFinder' is not installed. Install with: BiocManager::install('scDblFinder')",
            call. = FALSE
        )
    }
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
        !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop(
            "Doublet detection requires SingleCellExperiment and SummarizedExperiment.",
            call. = FALSE
        )
    }

    counts <- get_counts_matrix(obj)
    counts <- counts[, colnames(obj), drop = FALSE]
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = obj@meta.data[colnames(obj), , drop = FALSE]
    )
    sample_values <- NULL
    if (!is.null(sample_col) && sample_col %in% colnames(obj@meta.data)) {
        sample_values <- obj[[sample_col]][, 1]
    }

    sce <- scDblFinder::scDblFinder(
        sce,
        samples = sample_values
    )

    obj$scDblFinder_class <- sce$scDblFinder.class
    obj$scDblFinder_score <- sce$scDblFinder.score
    obj
}
