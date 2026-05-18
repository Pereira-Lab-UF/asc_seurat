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

    SeuratObject::DefaultAssay(norm_obj) <- "RNA"
    expect_gt(ncol(Seurat::GetAssayData(norm_obj, layer = "data")), 0)
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
