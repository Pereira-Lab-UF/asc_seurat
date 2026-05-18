test_that("clustering module returns a clustered Seurat object", {
    obj <- ascseurat:::run_normalization(demo_fixture(), method = "lognorm", nfeatures = 500)

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_clustering = 1)
            expect_s4_class(active_obj(), "Seurat")
            expect_true("tsne" %in% names(active_obj()@reductions))
            expect_equal(attr(active_obj(), "ascseurat_clustering_n_pcs", exact = TRUE), 10L)

            session$setInputs(rename_1 = "Renamed cluster", apply_rename = 1)
            expect_true(any(grepl("Renamed cluster", levels(SeuratObject::Idents(active_obj())))))
        }
    )
})
