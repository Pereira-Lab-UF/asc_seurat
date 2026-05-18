test_that("QC module can pass the loaded object through with no filtering", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(session$setInputs(filter_mode = "none"))
            expect_s4_class(qc_result(), "Seurat")
            expect_equal(ncol(qc_result()), ncol(obj))
            expect_identical(
                attr(qc_result(), "ascseurat_filter_mode", exact = TRUE),
                "none"
            )
        }
    )
})

test_that("QC module applies valid filtering thresholds", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(session$setInputs(
                filter_mode = "filter",
                min_genes = 0,
                max_genes = max(obj$nFeature_RNA, na.rm = TRUE) + 1,
                max_mito = 100,
                detect_doublets = FALSE,
                apply_filter = 1
            ))

            filtered <- suppressWarnings(qc_result())
            expect_s4_class(filtered, "Seurat")
            expect_equal(ncol(filtered), ncol(obj))
            expect_identical(
                attr(filtered, "ascseurat_filter_mode", exact = TRUE),
                "filter"
            )
        }
    )
})

test_that("QC module can apply no metric filtering when requested", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(session$setInputs(
                filter_mode = "none",
                detect_doublets = FALSE,
                apply_filter = 1
            ))

            filtered <- suppressWarnings(qc_result())
            expect_s4_class(filtered, "Seurat")
            expect_equal(ncol(filtered), ncol(obj))
            expect_identical(
                attr(filtered, "ascseurat_filter_mode", exact = TRUE),
                "none"
            )
        }
    )
})

test_that("QC module blocks thresholds that remove every cell", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(session$setInputs(
                filter_mode = "filter",
                min_genes = max(obj$nFeature_RNA, na.rm = TRUE) + 1,
                max_genes = max(obj$nFeature_RNA, na.rm = TRUE) + 2,
                max_mito = 100,
                detect_doublets = FALSE,
                apply_filter = 1
            ))

            expect_null(filtered_obj())
            expect_error(suppressWarnings(qc_result()), class = "shiny.silent.error")
        }
    )
})
