test_that("normalization and clustering work on the demo object", {
    obj <- demo_fixture()

    norm_obj <- ascseurat:::run_normalization(
        obj,
        method = "lognorm",
        nfeatures = 500,
        vars_to_regress = "percent.mt",
        score_cell_cycle = TRUE,
        cc_organism = "human",
        regress_cc = TRUE
    )

    expect_true("pca" %in% names(norm_obj@reductions))
    expect_true(all(c("S.Score", "G2M.Score", "Phase") %in% colnames(norm_obj@meta.data)))
    expect_identical(
        norm_obj@misc$ascseurat_scale_regress,
        c("percent.mt", "S.Score", "G2M.Score")
    )
    # Replayed when a subset is reanalyzed.
    expect_identical(
        norm_obj@misc$ascseurat_normalization,
        list(
            method = "lognorm",
            nfeatures = 500,
            scale_factor = 10000,
            vars_to_regress = "percent.mt",
            score_cell_cycle = TRUE,
            cc_organism = "human",
            regress_cc = TRUE
        )
    )
    expect_identical(norm_obj@misc$ascseurat_scale_cells, as.integer(ncol(obj)))

    clust_obj <- ascseurat:::run_clustering(norm_obj, dims = 1:10, resolution = 0.4, run_tsne = TRUE)
    expect_true("umap" %in% names(clust_obj@reductions))
    expect_true("tsne" %in% names(clust_obj@reductions))
    expect_gt(length(levels(SeuratObject::Idents(clust_obj))), 0)
    expect_equal(attr(clust_obj, "ascseurat_clustering_reduction"), "pca")
})

test_that("SCTransform keeps RNA prepared for downstream DE and plots", {
    obj <- demo_fixture()

    norm_obj <- ascseurat:::run_normalization(
        obj,
        method = "sctransform",
        nfeatures = 500
    )

    expect_true("SCT" %in% names(norm_obj@assays))
    expect_true("RNA" %in% names(norm_obj@assays))
    expect_equal(SeuratObject::DefaultAssay(norm_obj), "SCT")
    expect_identical(norm_obj@misc$ascseurat_scale_regress, character(0))

    SeuratObject::DefaultAssay(norm_obj) <- "RNA"
    expect_gt(ncol(Seurat::GetAssayData(norm_obj, layer = "data")), 0)

    # RNA scaled values are computed on demand for plots.
    genes <- head(SeuratObject::VariableFeatures(norm_obj), 3)
    scaled <- ascseurat:::fetch_expression_layer(norm_obj, genes, "scale.data")
    expect_identical(colnames(scaled), genes)
})

test_that("clustering can use an integrated reduction", {
    obj <- ascseurat:::run_normalization(demo_fixture(), method = "lognorm", nfeatures = 500)
    obj[["integrated.rpca"]] <- obj[["pca"]]

    clust_obj <- ascseurat:::run_clustering(
        obj,
        dims = 1:10,
        resolution = 0.4,
        reduction = "integrated.rpca"
    )

    expect_equal(attr(clust_obj, "ascseurat_clustering_reduction"), "integrated.rpca")
})

test_that("normalization and clustering work on small subsets", {
    obj <- demo_fixture()

    small <- suppressWarnings(ascseurat:::run_normalization(
        obj[, 1:30],
        method = "lognorm",
        nfeatures = 500
    ))
    expect_equal(ncol(Seurat::Embeddings(small, "pca")), 29)
    clust_small <- ascseurat:::run_clustering(small, dims = 1:10, resolution = 0.4, run_tsne = TRUE)
    expect_true(all(c("umap", "tsne") %in% names(clust_small@reductions)))

    # 60 cells support 50 PCs but not the default tSNE perplexity of 30.
    medium <- suppressWarnings(ascseurat:::run_normalization(
        obj[, 1:60],
        method = "lognorm",
        nfeatures = 500
    ))
    clust_medium <- ascseurat:::run_clustering(medium, dims = 1:10, resolution = 0.4, run_tsne = TRUE)
    expect_true("tsne" %in% names(clust_medium@reductions))
})

test_that("Harmony integration uses the requested PCs and records its parameters", {
    obj <- demo_fixture()
    obj$samples <- rep(c("A", "B"), length.out = ncol(obj))
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$samples)

    int_obj <- suppressWarnings(ascseurat:::run_integration(
        obj,
        method = "harmony",
        dims = 1:10,
        sample_col = "samples",
        nfeatures = 500
    ))

    expect_equal(ncol(Seurat::Embeddings(int_obj, "harmony")), 10)
    expected_record <- list(
        method = "harmony",
        dims = 1:10,
        normalization_method = "lognorm",
        sample_col = "samples"
    )
    expect_identical(int_obj@misc$ascseurat_integration, expected_record)
    expect_identical(attr(int_obj, "ascseurat_default_reduction", exact = TRUE), "harmony")
    expect_identical(int_obj@misc$ascseurat_scale_cells, as.integer(ncol(obj)))
    # The record has to survive subsetting for subset re-integration.
    sub_obj <- subset(int_obj, cells = colnames(int_obj)[1:100])
    expect_identical(sub_obj@misc$ascseurat_integration, expected_record)

    # Normalizing again (e.g. an integrated .rds in the Single Sample tab)
    # replaces the batch correction with a fresh PCA.
    renorm <- suppressWarnings(ascseurat:::run_normalization(int_obj, nfeatures = 500))
    expect_false("harmony" %in% names(renorm@reductions))
    expect_null(renorm@misc$ascseurat_integration)
    expect_null(attr(renorm, "ascseurat_default_reduction", exact = TRUE))
})

test_that("RPCA integration keeps its reduction and records its parameters", {
    obj <- demo_fixture()
    obj$samples <- rep(c("A", "B"), length.out = ncol(obj))
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$samples)

    int_obj <- suppressWarnings(ascseurat:::run_integration(
        obj,
        method = "rpca",
        dims = 1:10,
        sample_col = "samples",
        nfeatures = 500
    ))

    expect_true("integrated.rpca" %in% names(int_obj@reductions))
    expect_identical(attr(int_obj, "ascseurat_default_reduction", exact = TRUE), "integrated.rpca")
    expect_identical(int_obj@misc$ascseurat_integration$method, "rpca")
    expect_identical(int_obj@misc$ascseurat_scale_cells, as.integer(ncol(obj)))

    renorm <- suppressWarnings(ascseurat:::run_normalization(int_obj, nfeatures = 500))
    expect_false("integrated.rpca" %in% names(renorm@reductions))
    expect_null(renorm@misc$ascseurat_integration)
    expect_null(attr(renorm, "ascseurat_default_reduction", exact = TRUE))
})
