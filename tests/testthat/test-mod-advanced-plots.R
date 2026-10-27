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

renamed_fixture <- function() {
    obj <- clustered_fixture()
    SeuratObject::Idents(obj) <- paste0("c", as.character(SeuratObject::Idents(obj)))
    obj
}

test_that("Advanced Plots are outdated once the source tab replaces or clears its object", {
    local_mocked_bindings(
        showNotification = function(ui, ...) invisible("id"),
        .package = "ascseurat"
    )
    obj <- clustered_fixture()
    renamed <- renamed_fixture()
    genes <- head(SeuratObject::VariableFeatures(obj), 3)
    marker_file <- tempfile("advanced_markers_", fileext = ".csv")
    utils::write.csv(data.frame(gene = genes), marker_file, row.names = FALSE)
    current_obj <- shiny::reactiveVal(obj)
    # Re-runs the upstream reactive without changing its object.
    refresh <- shiny::reactiveVal(0)

    shiny::testServer(
        ascseurat::mod_advanced_plots_server,
        args = list(seurat_obj_single = shiny::reactive({
            refresh()
            shiny::req(current_obj())
            current_obj()
        })),
        {
            session$setInputs(
                data_source = "single",
                plot_type = "dotplot",
                gene_count_mode = "all",
                gene_list = list(datapath = marker_file),
                make_plot = 1
            )
            expect_s3_class(adv_plot(), "ggplot")

            # The same object again keeps the plot.
            refresh(refresh() + 1)
            session$flushReact()
            expect_identical(output$plot_ready, "true")
            expect_s3_class(adv_plot(), "ggplot")

            # Renamed or re-clustered object in the Single Sample tab.
            current_obj(renamed)
            session$flushReact()
            expect_false(session$isClosed())
            expect_error(plot_state(), "has changed or is no longer available")
            expect_error(output$plot_and_ordering, "has changed or is no longer available")
            expect_identical(output$plot_ready, "false")

            # Generating the plot again uses the current object.
            session$setInputs(make_plot = 2)
            expect_identical(output$plot_ready, "true")
            expect_identical(
                levels(SeuratObject::Idents(plot_state()$obj)),
                levels(SeuratObject::Idents(renamed))
            )

            # Choosing another source does not clear a plot that is current.
            session$setInputs(data_source = "upload")
            expect_identical(output$plot_ready, "true")

            # The Single Sample object is no longer available.
            current_obj(NULL)
            session$flushReact()
            expect_false(session$isClosed())
            expect_identical(output$plot_ready, "false")
        }
    )
})

test_that("Advanced Plots of an uploaded object ignore the other tabs", {
    local_mocked_bindings(
        showNotification = function(ui, ...) invisible("id"),
        .package = "ascseurat"
    )
    obj <- clustered_fixture()
    genes <- head(SeuratObject::VariableFeatures(obj), 3)
    marker_file <- tempfile("advanced_markers_", fileext = ".csv")
    utils::write.csv(data.frame(gene = genes), marker_file, row.names = FALSE)
    current_obj <- shiny::reactiveVal(obj)

    shiny::testServer(
        ascseurat::mod_advanced_plots_server,
        args = list(seurat_obj_single = shiny::reactive({
            shiny::req(current_obj())
            current_obj()
        })),
        {
            session$setInputs(
                data_source = "upload",
                rds_file = list(datapath = clustered_fixture_path()),
                plot_type = "dotplot",
                gene_count_mode = "all",
                gene_list = list(datapath = marker_file),
                make_plot = 1
            )
            expect_identical(output$plot_ready, "true")
            current_obj(NULL)
            session$flushReact()
            expect_identical(output$plot_ready, "true")
            expect_s3_class(adv_plot(), "ggplot")
        }
    )
})

test_that("Advanced Plots downloads are hidden until there is a current plot", {
    skip_on_cran()
    skip_if_not_installed("shinytest2")
    skip_if_not_installed("chromote")
    skip_if(is.null(suppressMessages(chromote::find_chrome())), "Chrome is not available")
    root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
    skip_if_not(file.exists(file.path(root, "DESCRIPTION")), "Package sources are not available")

    genes <- head(SeuratObject::VariableFeatures(clustered_fixture()), 3)
    marker_file <- tempfile("advanced_markers_", fileext = ".csv")
    utils::write.csv(data.frame(gene = genes), marker_file, row.names = FALSE)
    app_dir <- tempfile("advanced_app_")
    dir.create(app_dir)
    writeLines(c(
        sprintf("pkgload::load_all(%s, quiet = TRUE)", deparse(root)),
        sprintf("obj <- readRDS(%s)", deparse(clustered_fixture_path())),
        "renamed <- obj",
        "SeuratObject::Idents(renamed) <- paste0('c', as.character(SeuratObject::Idents(obj)))",
        "ui <- bslib::page_fluid(shiny::actionButton('swap', 'Swap'), mod_advanced_plots_ui('adv'))",
        "server <- function(input, output, session) {",
        "    current_obj <- shiny::reactiveVal(obj)",
        "    shiny::observeEvent(input$swap, current_obj(renamed))",
        "    mod_advanced_plots_server('adv', seurat_obj_single = current_obj)",
        "}",
        "shiny::shinyApp(ui, server)"
    ), file.path(app_dir, "app.R"))

    app <- shinytest2::AppDriver$new(
        app_dir,
        load_timeout = 120000,
        timeout = 60000,
        shiny_args = list(port = httpuv::randomPort(min = 5000, max = 5900))
    )
    on.exit(app$stop(), add = TRUE)
    visible <- "$('#adv-download_per_gene').is(':visible')"

    expect_false(app$get_js(visible))
    app$set_inputs(`adv-data_source` = "single", `adv-plot_type` = "dotplot")
    app$upload_file(`adv-gene_list` = marker_file)
    app$click("adv-make_plot")
    app$wait_for_js(visible, timeout = 60000)
    expect_true(app$get_js(visible))

    app$click("swap")
    app$wait_for_js(sprintf("!%s", visible), timeout = 60000)
    expect_false(app$get_js(visible))
    expect_match(app$get_text("#adv-plot_and_ordering"), "has changed or is no longer available")
})
