test_that("SingleR module runs when cached references are available", {
    skip_if_not_installed("SingleR")
    skip_if_not_installed("celldex")

    cache_file <- file.path(tools::R_user_dir("ascseurat", "cache"), "celldex", "hpca.rds")
    skip_if_not(file.exists(cache_file), "SingleR reference cache is not available")

    shiny::testServer(
        ascseurat::mod_celltype_annotation_server,
        args = list(
            seurat_obj_single = shiny::reactive(clustered_fixture()),
            seurat_obj_multi = shiny::reactive(clustered_fixture())
        ),
        {
            session$setInputs(
                data_source = "single",
                reference = "hpca",
                granularity = "label.main",
                run_singler = 1
            )
            expect_true("singler_label" %in% colnames(annotated_obj()@meta.data))
        }
    )
})

test_that("SingleR expression prep joins split assay layers", {
    skip_if_not("JoinLayers" %in% getNamespaceExports("SeuratObject"))

    counts <- matrix(rpois(10 * 8, 5), nrow = 10)
    rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
    colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))

    obj <- Seurat::CreateSeuratObject(counts)
    obj$sample <- rep(c("A", "B"), each = 4)
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)

    expect_error(
        Seurat::GetAssayData(obj, assay = "RNA", layer = "data"),
        "multiple layers"
    )

    prepared <- ascseurat:::prepare_singler_input(obj)
    expect_equal(dim(prepared$expr_matrix), c(nrow(obj), ncol(obj)))
    expect_true("data" %in% SeuratObject::Layers(prepared$obj[["RNA"]]))
    expect_false(any(grepl("^data[.]", SeuratObject::Layers(prepared$obj[["RNA"]]))))
})
