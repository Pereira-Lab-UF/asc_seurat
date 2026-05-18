#' Run QC on a Seurat Object
#'
#' Adds percentage mitochondrial content and filters cells.
#'
#' @param obj A Seurat object.
#' @param mito_pattern Character. Regex to identify mitochondrial genes.
#' @param min_features Integer. Min genes per cell.
#' @param max_features Integer. Max genes per cell.
#' @param max_mito_pct Numeric. Max mitochondrial percentage.
#' @return A filtered Seurat object.
#' @keywords internal
run_qc_filter <- function(obj, mito_pattern = "^MT-",
                          min_features = 200, max_features = Inf,
                          max_mito_pct = 5) {
    # Calculate mito percentage
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mito_pattern)

    # Filter
    obj <- subset(obj,
                  subset = nFeature_RNA > min_features &
                      nFeature_RNA < max_features &
                      percent.mt < max_mito_pct)

    message("After filtering: ", ncol(obj), " cells remain.")
    return(obj)
}


#' Run Normalization Pipeline
#'
#' Normalizes data using either LogNormalize or SCTransform v2.
#'
#' @param obj A Seurat object.
#' @param method Character. "lognorm" or "sctransform".
#' @param nfeatures Integer. Number of variable features.
#' @param scale_factor Numeric. Scale factor for LogNormalize.
#' @param vars_to_regress Optional character vector of metadata columns to regress.
#' @param score_cell_cycle Logical. Whether to score cell-cycle phase.
#' @param cc_organism Character. "human" or "mouse" for cell-cycle genes.
#' @param regress_cc Logical. Whether to regress out cell-cycle scores.
#' @return A normalized Seurat object with PCA computed.
#' @keywords internal
run_normalization <- function(obj, method = "lognorm",
                              nfeatures = 3000, scale_factor = 10000,
                              vars_to_regress = NULL,
                              score_cell_cycle = FALSE,
                              cc_organism = "human",
                              regress_cc = FALSE) {
    score_cc_metadata <- function(target_obj) {
        s_genes <- Seurat::cc.genes.updated.2019$s.genes
        g2m_genes <- Seurat::cc.genes.updated.2019$g2m.genes

        if (identical(cc_organism, "mouse")) {
            to_title_case <- function(values) {
                lower <- tolower(values)
                paste0(toupper(substring(lower, 1, 1)), substring(lower, 2))
            }
            s_genes <- to_title_case(s_genes)
            g2m_genes <- to_title_case(g2m_genes)
        }

        s_genes <- intersect(s_genes, rownames(target_obj))
        g2m_genes <- intersect(g2m_genes, rownames(target_obj))
        if (!length(s_genes) || !length(g2m_genes)) {
            target_obj$S.Score <- 0
            target_obj$G2M.Score <- 0
            target_obj$Phase <- "Undetected"
            return(target_obj)
        }

        Seurat::CellCycleScoring(
            target_obj,
            s.features = s_genes,
            g2m.features = g2m_genes,
            ctrl = min(100, max(1, floor(nrow(target_obj) / 100))),
            set.ident = FALSE
        )
    }

    if (method == "sctransform") {
        message("Running SCTransform v2...")
        if (isTRUE(score_cell_cycle)) {
            temp_obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                                      scale.factor = scale_factor, verbose = FALSE)
            temp_obj <- score_cc_metadata(temp_obj)
            obj$S.Score <- temp_obj$S.Score
            obj$G2M.Score <- temp_obj$G2M.Score
            obj$Phase <- temp_obj$Phase
        }
        if (isTRUE(regress_cc)) {
            vars_to_regress <- unique(c(vars_to_regress, "S.Score", "G2M.Score"))
        }
        obj <- SCTransform(obj, vst.flavor = "v2",
                           variable.features.n = nfeatures,
                           vars.to.regress = vars_to_regress,
                           verbose = FALSE)
        obj <- prepare_rna_assay_for_de(
            obj,
            nfeatures = nfeatures,
            scale_factor = scale_factor
        )
        SeuratObject::DefaultAssay(obj) <- "SCT"
    } else {
        message("Running LogNormalize pipeline...")
        if ("RNA" %in% names(obj@assays)) {
            SeuratObject::DefaultAssay(obj) <- "RNA"
        }
        obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                            scale.factor = scale_factor, verbose = FALSE)
        obj <- FindVariableFeatures(obj, selection.method = "vst",
                                    nfeatures = nfeatures, verbose = FALSE)
        if (isTRUE(score_cell_cycle)) {
            obj <- score_cc_metadata(obj)
        }
        if (isTRUE(regress_cc)) {
            vars_to_regress <- unique(c(vars_to_regress, "S.Score", "G2M.Score"))
        }
        obj <- ScaleData(obj, vars.to.regress = vars_to_regress, verbose = FALSE)
    }

    # PCA
    message("Running PCA...")
    obj <- RunPCA(obj, verbose = FALSE)

    return(obj)
}


#' Run Clustering Pipeline
#'
#' Performs neighbor finding, clustering, and dimensionality reduction.
#'
#' @param obj A Seurat object with PCA computed.
#' @param dims Integer vector. PCA dimensions to use.
#' @param resolution Numeric. Clustering resolution.
#' @param run_tsne Logical. Whether to also run tSNE.
#' @param reduction Character. Dimensionality reduction to use for neighbors
#'   and UMAP. Defaults to "pca".
#' @return A clustered Seurat object with UMAP (and optionally tSNE).
#' @keywords internal
run_clustering <- function(obj, dims = 1:30, resolution = 0.6,
                           run_tsne = FALSE, reduction = "pca") {
    if (!reduction %in% names(obj@reductions)) {
        stop("Reduction not found in Seurat object: ", reduction, call. = FALSE)
    }
    dims <- valid_reduction_dims(obj, reduction, dims)

    message("Finding neighbors...")
    obj <- FindNeighbors(obj, reduction = reduction, dims = dims, verbose = FALSE)

    message("Clustering at resolution ", resolution, "...")
    obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)

    message("Running UMAP...")
    obj <- RunUMAP(obj, reduction = reduction, dims = dims, verbose = FALSE)

    if (run_tsne) {
        message("Running tSNE...")
        obj <- RunTSNE(obj, reduction = reduction, dims = dims)
    }

    n_clusters <- length(levels(Idents(obj)))
    message("Found ", n_clusters, " clusters.")
    attr(obj, "ascseurat_clustering_reduction") <- reduction

    return(obj)
}


#' Run Integration with RPCA or Harmony
#'
#' Integrates multiple samples using either RPCA or Harmony.
#'
#' @param obj A merged Seurat object with multiple samples.
#' @param method Character. "rpca" or "harmony".
#' @param dims Integer vector. Dimensions for integration.
#' @param sample_col Character. Column in metadata identifying samples.
#' @param normalization_method Character. "lognorm" or "sctransform".
#' @param nfeatures Integer. Number of variable features.
#' @param scale_factor Numeric. Scale factor for LogNormalize.
#' @return An integrated Seurat object.
#' @keywords internal
run_integration <- function(obj, method = "rpca", dims = 1:30,
                            sample_col = "orig.ident",
                            normalization_method = "lognorm",
                            nfeatures = 3000,
                            scale_factor = 10000) {
    if (!sample_col %in% colnames(obj@meta.data)) {
        stop("Metadata column not found for integration: ", sample_col, call. = FALSE)
    }

    sample_values <- as.character(obj[[sample_col]][, 1])
    sample_values <- sample_values[!is.na(sample_values) & nzchar(sample_values)]
    if (length(unique(sample_values)) < 2) {
        stop("Integration requires at least two distinct samples in '", sample_col, "'.",
             call. = FALSE)
    }
    min_sample_cells <- min(table(sample_values))

    memory_guard <- local_future_globals_max_size()
    on.exit(options(memory_guard$old), add = TRUE)

    message(
        "Setting future.globals.maxSize to ",
        format_bytes(memory_guard$target),
        " (available memory: ",
        format_bytes(memory_guard$available),
        ")."
    )

    if (method == "harmony") {
        message("Integrating with Harmony...")
        if (!requireNamespace("harmony", quietly = TRUE)) {
            stop("Package 'harmony' required. Install with: install.packages('harmony')",
                 call. = FALSE)
        }
        obj <- run_normalization(
            obj,
            method = normalization_method,
            nfeatures = nfeatures,
            scale_factor = scale_factor
        )
        obj <- RunPCA(obj, verbose = FALSE)
        obj <- harmony::RunHarmony(obj, group.by.vars = sample_col)
        integration_dims <- valid_reduction_dims(obj, "harmony", dims)
        obj <- FindNeighbors(obj, reduction = "harmony",
                             dims = integration_dims, verbose = FALSE)
        obj <- RunUMAP(obj, reduction = "harmony",
                       dims = integration_dims, verbose = FALSE)
        attr(obj, "ascseurat_default_reduction") <- "harmony"
    } else {
        message("Integrating with RPCA via IntegrateLayers...")
        obj <- run_normalization(
            obj,
            method = normalization_method,
            nfeatures = nfeatures,
            scale_factor = scale_factor
        )
        obj <- RunPCA(obj, verbose = FALSE)

        integration_assay <- if (identical(normalization_method, "sctransform") &&
            "SCT" %in% names(obj@assays)) {
            "SCT"
        } else {
            SeuratObject::DefaultAssay(obj)
        }
        normalization_label <- if (identical(normalization_method, "sctransform")) {
            "SCT"
        } else {
            "LogNormalize"
        }
        pca_dims <- valid_reduction_dims(obj, "pca", dims)
        obj <- Seurat::IntegrateLayers(
            object = obj,
            method = Seurat::RPCAIntegration,
            orig.reduction = "pca",
            new.reduction = "integrated.rpca",
            assay = integration_assay,
            normalization.method = normalization_label,
            dims = pca_dims,
            k.weight = min(100, max(1, floor(as.integer(min_sample_cells) / 4))),
            verbose = FALSE
        )
        integration_dims <- valid_reduction_dims(obj, "integrated.rpca", dims)
        obj <- FindNeighbors(obj, reduction = "integrated.rpca",
                             dims = integration_dims, verbose = FALSE)
        obj <- RunUMAP(obj, reduction = "integrated.rpca",
                       dims = integration_dims, verbose = FALSE)
        attr(obj, "ascseurat_default_reduction") <- "integrated.rpca"
    }

    return(obj)
}


#' Prepare RNA Assay for Differential Expression and Visualization
#'
#' SCTransform is used through PCA/clustering, while RNA-normalized values are
#' kept available for DE and expression plots.
#'
#' @param obj A Seurat object.
#' @param nfeatures Integer. Number of variable features.
#' @param scale_factor Numeric. Scale factor for LogNormalize.
#' @return A Seurat object with RNA data and scale.data prepared when possible.
#' @keywords internal
prepare_rna_assay_for_de <- function(obj, nfeatures = 3000,
                                     scale_factor = 10000) {
    if (!"RNA" %in% names(obj@assays)) {
        return(obj)
    }

    current_assay <- SeuratObject::DefaultAssay(obj)
    on.exit(SeuratObject::DefaultAssay(obj) <- current_assay, add = TRUE)

    SeuratObject::DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(
        obj,
        normalization.method = "LogNormalize",
        scale.factor = scale_factor,
        verbose = FALSE
    )
    obj <- FindVariableFeatures(
        obj,
        selection.method = "vst",
        nfeatures = nfeatures,
        verbose = FALSE
    )
    obj <- ScaleData(obj, verbose = FALSE)
    obj
}


valid_reduction_dims <- function(obj, reduction, dims) {
    max_dims <- ncol(Seurat::Embeddings(obj, reduction))
    dims <- dims[dims <= max_dims]
    if (!length(dims)) {
        stop("No valid dimensions available for reduction: ", reduction, call. = FALSE)
    }
    dims
}
