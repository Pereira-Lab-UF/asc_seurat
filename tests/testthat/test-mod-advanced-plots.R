test_that("Advanced Plots ordering is populated after plot generation", {
    obj <- clustered_fixture()
    genes <- head(rownames(obj), 4)
    marker_file <- tempfile("advanced_markers_", fileext = ".csv")
    utils::write.csv(data.frame(gene = genes), marker_file, row.names = FALSE)

    shiny::testServer(
        ascseurat::mod_advanced_plots_server,
        args = list(
            seurat_obj_single = shiny::reactive(obj),
            seurat_obj_multi = shiny::reactive(obj)
        ),
        {
            session$setInputs(
                data_source = "single",
                plot_type = "dotplot",
                gene_count_mode = "all",
                gene_list = list(datapath = marker_file),
                make_plot = 1
            )

            expect_equal(gene_order(), genes)
            expect_equal(applied_gene_order(), genes)
            expect_equal(selected_genes(), genes)

            session$setInputs(
                order_move = list(
                    type = "gene",
                    direction = "down",
                    index = 1,
                    nonce = 1
                )
            )
            expect_equal(gene_order()[1:2], genes[c(2, 1)])
            expect_equal(selected_genes(), genes)

            session$setInputs(apply_order = 1)
            expect_equal(applied_gene_order()[1:2], genes[c(2, 1)])
            expect_equal(selected_genes()[1:2], genes[c(2, 1)])

            initial_clusters <- cluster_order()
            initial_applied_clusters <- applied_cluster_order()
            if (length(initial_clusters) > 1) {
                session$setInputs(
                    order_move = list(
                        type = "cluster",
                        direction = "down",
                        index = 1,
                        nonce = 2
                    )
                )
                expect_equal(cluster_order()[1:2], initial_clusters[c(2, 1)])
                expect_equal(applied_cluster_order(), initial_applied_clusters)

                session$setInputs(apply_order = 2)
                expect_equal(applied_cluster_order()[1:2], initial_clusters[c(2, 1)])
            }
        }
    )
})

advanced_plot_state <- function(obj, genes) {
    marker_file <- tempfile("advanced_markers_", fileext = ".csv")
    utils::write.csv(data.frame(gene = genes), marker_file, row.names = FALSE)
    state <- NULL
    shiny::testServer(
        ascseurat::mod_advanced_plots_server,
        args = list(seurat_obj_single = shiny::reactive(obj)),
        {
            session$setInputs(
                data_source = "single",
                plot_type = "dotplot",
                gene_count_mode = "all",
                gene_list = list(datapath = marker_file),
                make_plot = 1
            )
            expect_false(session$isClosed())
            expect_s3_class(adv_plot(), "ggplot")
            state <<- plot_state()
        }
    )
    state
}

unnormalized_fixture <- function() {
    obj <- clustered_fixture()
    raw <- SeuratObject::CreateSeuratObject(
        counts = SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
    )
    SeuratObject::Idents(raw) <- SeuratObject::Idents(obj)[colnames(raw)]
    raw
}

test_that("Advanced Plots normalize an RNA assay that has only counts", {
    notifications <- character()
    local_mocked_bindings(
        showNotification = function(ui, ...) {
            notifications <<- c(notifications, paste(ui, collapse = " "))
            invisible("id")
        },
        .package = "ascseurat"
    )
    raw <- unnormalized_fixture()
    expect_identical(SeuratObject::Layers(raw[["RNA"]]), "counts")
    genes <- head(SeuratObject::VariableFeatures(clustered_fixture()), 3)

    state <- advanced_plot_state(raw, genes)
    expect_identical(SeuratObject::DefaultAssay(state$obj), "RNA")
    expected <- Seurat::NormalizeData(raw, verbose = FALSE)
    expect_equal(
        SeuratObject::LayerData(state$obj, assay = "RNA", layer = "data"),
        SeuratObject::LayerData(expected, assay = "RNA", layer = "data")
    )
    expect_true(any(grepl("log-normalized", notifications, fixed = TRUE)))

    # v4-style objects keep raw counts in the data layer until normalized.
    notifications <- character()
    v4_style <- raw
    SeuratObject::LayerData(v4_style, assay = "RNA", layer = "data") <-
        SeuratObject::LayerData(raw, assay = "RNA", layer = "counts")
    state <- advanced_plot_state(v4_style, genes)
    expect_equal(
        SeuratObject::LayerData(state$obj, assay = "RNA", layer = "data"),
        SeuratObject::LayerData(expected, assay = "RNA", layer = "data")
    )
    expect_true(any(grepl("log-normalized", notifications, fixed = TRUE)))
})

test_that("Advanced Plots use SCT when the RNA assay is not normalized", {
    notifications <- character()
    local_mocked_bindings(
        showNotification = function(ui, ...) {
            notifications <<- c(notifications, paste(ui, collapse = " "))
            invisible("id")
        },
        .package = "ascseurat"
    )
    obj <- clustered_fixture()
    sct <- unnormalized_fixture()
    # A data-only stand-in for SCT; Seurat warns that it has no counts layer.
    suppressWarnings(sct[["SCT"]] <- SeuratObject::CreateAssayObject(
        data = SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
    ))
    genes <- head(SeuratObject::VariableFeatures(obj), 3)

    state <- advanced_plot_state(sct, genes)
    expect_identical(SeuratObject::DefaultAssay(state$obj), "SCT")
    expect_identical(SeuratObject::Layers(state$obj[["RNA"]]), "counts")
    expect_true(any(grepl("SCT assay", notifications, fixed = TRUE)))
})

test_that("Advanced Plots keep the normalized RNA assay of app objects", {
    notifications <- character()
    local_mocked_bindings(
        showNotification = function(ui, ...) {
            notifications <<- c(notifications, paste(ui, collapse = " "))
            invisible("id")
        },
        .package = "ascseurat"
    )
    obj <- clustered_fixture()
    genes <- head(SeuratObject::VariableFeatures(obj), 3)

    state <- advanced_plot_state(obj, genes)
    expect_identical(SeuratObject::DefaultAssay(state$obj), "RNA")
    expect_identical(
        SeuratObject::LayerData(state$obj, assay = "RNA", layer = "data"),
        SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
    )

    # SCTransform objects from the app also keep a normalized RNA assay.
    with_sct <- obj
    suppressWarnings(with_sct[["SCT"]] <- SeuratObject::CreateAssayObject(
        data = SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
    ))
    SeuratObject::DefaultAssay(with_sct) <- "SCT"
    state <- advanced_plot_state(with_sct, genes)
    expect_identical(SeuratObject::DefaultAssay(state$obj), "RNA")
    expect_length(notifications, 0)
})
