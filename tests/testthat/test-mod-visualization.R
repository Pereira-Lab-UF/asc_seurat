test_that("visualization helpers fetch the requested expression layer", {
    obj <- suppressWarnings(clustered_fixture())
    genes <- head(rownames(obj), 5)

    normalized <- ascseurat:::fetch_expression_layer(obj, genes, "data")
    scaled <- ascseurat:::fetch_expression_layer(obj, genes, "scale.data")

    expect_equal(colnames(normalized), colnames(scaled))
    expect_false(isTRUE(all.equal(normalized, scaled)))
})

test_that("clustered expression dot plot is generated from selected layer", {
    obj <- suppressWarnings(clustered_fixture())
    genes <- head(rownames(obj), 5)

    plot <- ascseurat:::clustered_expression_dot_plot(obj, genes, "scale.data")

    expect_s3_class(plot, "ggplot")
})

test_that("clustered expression dot plot can split clusters by sample", {
    obj <- suppressWarnings(clustered_fixture())
    obj$samples <- rep(c("sample_A", "sample_B"), length.out = ncol(obj))
    genes <- head(rownames(obj), 5)

    expect_identical(ascseurat:::visualization_sample_col(obj), "samples")
    plot <- ascseurat:::clustered_expression_dot_plot(
        obj,
        genes,
        "data",
        sample_col = "samples"
    )

    split_labels <- levels(plot$data$cluster)
    expect_true(any(grepl("sample_A", split_labels, fixed = TRUE)))
    expect_true(any(grepl("sample_B", split_labels, fixed = TRUE)))
    expect_identical(plot$labels$x, "Cluster | Sample")
})

test_that("feature plots split by sample work with more genes than samples", {
    obj <- suppressWarnings(clustered_fixture())
    obj$samples <- rep(c("sample_A", "sample_B"), length.out = ncol(obj))
    genes <- head(SeuratObject::VariableFeatures(obj), 3)
    markers <- data.frame(gene = genes, cluster = "0", p_val_adj = c(0.001, 0.002, 0.003))

    shiny::testServer(
        ascseurat::mod_visualization_server,
        args = list(
            seurat_obj = shiny::reactive(obj),
            markers_de = shiny::reactive(markers)
        ),
        {
            session$setInputs(
                marker_source = "de_results",
                top_n_genes = 3,
                slot_selection = "data",
                split_by_sample = TRUE,
                show_expression = 1
            )
            expect_identical(requested_sample_split_col(), "samples")
            expect_s3_class(feature_plot(), "patchwork")

            # Scaled split panels keep negative z-scores instead of being
            # clipped to the normalized-data range.
            session$setInputs(slot_selection = "scale.data", show_expression = 2)
            plotted <- unlist(lapply(seq_along(feature_plot()), function(i) {
                panel_data <- feature_plot()[[i]]$data
                unlist(panel_data[intersect(genes, colnames(panel_data))])
            }))
            expect_lt(min(plotted, na.rm = TRUE), 0)
        }
    )
})

test_that("invalid Top N genes is validated instead of closing the session", {
    obj <- suppressWarnings(clustered_fixture())
    genes <- head(SeuratObject::VariableFeatures(obj), 3)
    markers <- data.frame(gene = genes, cluster = "0", p_val_adj = c(0.001, 0.002, 0.003))
    top_n_message <- "Top N genes must be a whole number of at least 1."
    notifications <- character()
    local_mocked_bindings(
        showNotification = function(ui, ...) {
            notifications <<- c(notifications, paste(ui, collapse = " "))
            invisible("id")
        },
        .package = "ascseurat"
    )

    for (top_n in list(NA, -5, 0, 2.5)) {
        notifications <- character()
        shiny::testServer(
            ascseurat::mod_visualization_server,
            args = list(
                seurat_obj = shiny::reactive(obj),
                markers_de = shiny::reactive(markers)
            ),
            {
                session$setInputs(
                    marker_source = "de_results",
                    de_cluster = "0",
                    top_n_genes = top_n,
                    slot_selection = "data"
                )
                expect_error(output$de_genes_preview, top_n_message, fixed = TRUE)

                session$setInputs(show_expression = 1)
                expect_false(session$isClosed())
                expect_null(plot_request())
                expect_true(top_n_message %in% notifications)

                session$setInputs(top_n_genes = 2, show_expression = 2)
                expect_false(session$isClosed())
                expect_identical(plot_request()$genes, genes[1:2])
            }
        )
    }
})

test_that("clustered expression dot plot works with a single gene", {
    obj <- suppressWarnings(clustered_fixture())
    gene <- head(SeuratObject::VariableFeatures(obj), 1)
    n_clusters <- length(unique(SeuratObject::Idents(obj)))

    plot <- ascseurat:::clustered_expression_dot_plot(obj, gene, "data")
    expect_s3_class(plot, "ggplot")
    expect_equal(nrow(plot$data), n_clusters)
    expect_identical(unique(as.character(plot$data$gene)), gene)

    obj$samples <- rep(c("sample_A", "sample_B"), length.out = ncol(obj))
    plot <- ascseurat:::clustered_expression_dot_plot(
        obj,
        gene,
        "scale.data",
        sample_col = "samples"
    )
    expect_s3_class(plot, "ggplot")
    expect_equal(nrow(plot$data), n_clusters * 2)
})

test_that("scaled values do not depend on the other selected genes", {
    obj <- suppressWarnings(clustered_fixture())
    obj <- suppressWarnings(Seurat::ScaleData(
        obj,
        vars.to.regress = "nCount_RNA",
        verbose = FALSE
    ))
    obj@misc$ascseurat_scale_regress <- "nCount_RNA"
    obj@misc$ascseurat_scale_cells <- ncol(obj)
    scaled_genes <- rownames(SeuratObject::LayerData(obj, layer = "scale.data"))
    scaled_gene <- scaled_genes[[1]]
    unscaled <- setdiff(rownames(obj), scaled_genes)
    counts <- SeuratObject::LayerData(obj, layer = "counts")[unscaled, , drop = FALSE]
    unscaled_gene <- unscaled[[which.max(Matrix::rowSums(counts > 0))]]

    alone <- ascseurat:::fetch_expression_layer(obj, scaled_gene, "scale.data")
    together <- ascseurat:::fetch_expression_layer(
        obj,
        c(scaled_gene, unscaled_gene),
        "scale.data"
    )
    expect_equal(together[, scaled_gene], alone[, scaled_gene])

    # The missing gene is scaled with the regression used at normalization.
    expected <- suppressWarnings(Seurat::ScaleData(
        obj,
        features = unscaled_gene,
        vars.to.regress = "nCount_RNA",
        verbose = FALSE
    ))
    expected <- SeuratObject::LayerData(expected, layer = "scale.data")
    expect_equal(
        unname(together[, unscaled_gene]),
        unname(expected[unscaled_gene, rownames(together)])
    )

    # Only the requested rows are kept, and the existing one is unchanged.
    ensured <- ascseurat:::ensure_expression_layer(
        obj,
        c(scaled_gene, unscaled_gene),
        "scale.data"
    )
    ensured_scaled <- SeuratObject::LayerData(ensured, layer = "scale.data")
    expect_setequal(rownames(ensured_scaled), c(scaled_gene, unscaled_gene))
    expect_equal(
        ensured_scaled[scaled_gene, ],
        SeuratObject::LayerData(obj, layer = "scale.data")[scaled_gene, ]
    )
})

test_that("on-demand scaling ignores absent or unknown regression variables", {
    obj <- suppressWarnings(clustered_fixture())
    scaled_genes <- rownames(SeuratObject::LayerData(obj, layer = "scale.data"))
    gene <- setdiff(rownames(obj), scaled_genes)[[1]]
    expected <- suppressWarnings(Seurat::ScaleData(obj, features = gene, verbose = FALSE))
    expected <- SeuratObject::LayerData(expected, layer = "scale.data")[gene, ]

    obj@misc$ascseurat_scale_regress <- NULL
    without_misc <- ascseurat:::fetch_expression_layer(obj, gene, "scale.data")
    expect_equal(unname(without_misc[, gene]), unname(expected[rownames(without_misc)]))

    obj@misc$ascseurat_scale_regress <- "not_a_metadata_column"
    unknown_var <- ascseurat:::fetch_expression_layer(obj, gene, "scale.data")
    expect_equal(unknown_var, without_misc)
})

expected_scaled_values <- function(obj, genes, regress = NULL) {
    scaled <- suppressWarnings(Seurat::ScaleData(
        obj,
        features = genes,
        vars.to.regress = regress,
        verbose = FALSE
    ))
    t(as.matrix(SeuratObject::LayerData(scaled, layer = "scale.data")[genes, , drop = FALSE]))
}

test_that("scaled values of a subset that kept its clusters come from the subset cells", {
    obj <- suppressWarnings(clustered_fixture())
    obj@misc$ascseurat_scale_regress <- "nCount_RNA"
    obj@misc$ascseurat_scale_cells <- ncol(obj)
    scaled_genes <- rownames(SeuratObject::LayerData(obj, layer = "scale.data"))
    unscaled <- setdiff(rownames(obj), scaled_genes)
    counts <- SeuratObject::LayerData(obj, layer = "counts")[unscaled, , drop = FALSE]
    genes <- c(scaled_genes[[1]], unscaled[[which.max(Matrix::rowSums(counts > 0))]])
    # "Keep current reduction and clusters" subsets the object as is.
    keep <- as.character(SeuratObject::Idents(obj)) %in% c("0", "1")
    sub <- subset(obj, cells = colnames(obj)[keep])

    # scale.data was computed on other cells, so every requested gene is
    # rescaled on the subset with the recorded regression.
    for (requested in list(genes[[1]], genes)) {
        values <- ascseurat:::fetch_expression_layer(sub, requested, "scale.data")
        expected <- expected_scaled_values(sub, requested, "nCount_RNA")
        expect_equal(values, expected[rownames(values), colnames(values), drop = FALSE])
    }
    plot <- ascseurat:::clustered_expression_dot_plot(sub, genes, "scale.data")
    expect_s3_class(plot, "ggplot")

    # Without a record, scale.data is used as is when it has every gene...
    sub@misc$ascseurat_scale_cells <- NULL
    values <- ascseurat:::fetch_expression_layer(sub, genes[[1]], "scale.data")
    existing <- SeuratObject::LayerData(sub, layer = "scale.data")[genes[[1]], ]
    expect_equal(values[, genes[[1]]], existing[rownames(values)])

    # ...and all requested genes are rescaled together when one is missing.
    values <- ascseurat:::fetch_expression_layer(sub, genes, "scale.data")
    expected <- expected_scaled_values(sub, genes, "nCount_RNA")
    expect_equal(values, expected[rownames(values), colnames(values), drop = FALSE])
})

test_that("expression plots are cleared when the object or the DE results go away", {
    obj <- suppressWarnings(clustered_fixture())
    genes <- head(SeuratObject::VariableFeatures(obj), 3)
    markers <- data.frame(gene = genes, cluster = "0", p_val_adj = c(0.001, 0.002, 0.003))
    renamed <- obj
    SeuratObject::Idents(renamed) <- paste0("c", as.character(SeuratObject::Idents(obj)))
    current_obj <- shiny::reactiveVal(obj)
    current_markers <- shiny::reactiveVal(markers)
    # Re-runs the upstream reactive without changing its object.
    refresh <- shiny::reactiveVal(0)
    notifications <- character()
    local_mocked_bindings(
        showNotification = function(ui, ...) {
            notifications <<- c(notifications, paste(ui, collapse = " "))
            invisible("id")
        },
        .package = "ascseurat"
    )

    shiny::testServer(
        ascseurat::mod_visualization_server,
        args = list(
            seurat_obj = shiny::reactive({
                refresh()
                shiny::req(current_obj())
                current_obj()
            }),
            markers_de = shiny::reactive({
                shiny::req(current_markers())
                current_markers()
            })
        ),
        {
            session$setInputs(
                marker_source = "de_results",
                de_cluster = "0",
                top_n_genes = 3,
                slot_selection = "data",
                show_expression = 1
            )
            expect_identical(plot_request()$genes, genes)

            # The same object again keeps the plots.
            refresh(refresh() + 1)
            session$flushReact()
            expect_identical(plot_request()$genes, genes)

            # Renamed or re-clustered object: the plots are cleared and come
            # back only with a new click.
            current_obj(renamed)
            session$flushReact()
            expect_null(plot_request())
            expect_identical(output$plots_requested, "false")
            session$setInputs(show_expression = 2)
            expect_identical(plot_request()$genes, genes)
            expect_identical(output$plots_requested, "true")

            # Object no longer available, then available again.
            current_obj(NULL)
            session$flushReact()
            expect_null(plot_request())
            current_obj(obj)
            session$flushReact()
            expect_null(plot_request())
            session$setInputs(show_expression = 3)
            expect_identical(plot_request()$genes, genes)

            # The DE results the plots were built from are cleared.
            current_markers(NULL)
            session$flushReact()
            expect_null(plot_request())

            # Clicking with nothing to select warns instead of doing nothing.
            n_notes <- length(notifications)
            session$setInputs(show_expression = 4)
            expect_false(session$isClosed())
            expect_null(plot_request())
            expect_true(any(grepl(
                "No genes can be selected",
                notifications[-seq_len(n_notes)],
                fixed = TRUE
            )))
        }
    )
})

test_that("re-running DE with no markers clears the plots built from DE genes", {
    obj <- suppressWarnings(clustered_fixture())
    genes <- head(SeuratObject::VariableFeatures(obj), 3)
    markers <- data.frame(gene = genes, cluster = "0", p_val_adj = c(0.001, 0.002, 0.003))
    current_markers <- shiny::reactiveVal(markers)
    local_mocked_bindings(
        showNotification = function(ui, ...) invisible("id"),
        .package = "ascseurat"
    )

    shiny::testServer(
        ascseurat::mod_visualization_server,
        args = list(
            seurat_obj = shiny::reactive(obj),
            markers_de = shiny::reactive({
                shiny::req(current_markers())
                current_markers()
            })
        ),
        {
            session$setInputs(
                marker_source = "de_results",
                de_cluster = "0",
                top_n_genes = 3,
                slot_selection = "data",
                show_expression = 1
            )
            expect_identical(plot_request()$genes, genes)

            current_markers(markers[0, ])
            session$flushReact()
            expect_null(plot_request())
            expect_identical(output$plots_requested, "false")
        }
    )
})

test_that("an unreadable marker file warns instead of closing the session", {
    obj <- suppressWarnings(clustered_fixture())
    empty_file <- tempfile("viz_markers_", fileext = ".csv")
    file.create(empty_file)
    notifications <- character()
    local_mocked_bindings(
        showNotification = function(ui, ...) {
            notifications <<- c(notifications, paste(ui, collapse = " "))
            invisible("id")
        },
        .package = "ascseurat"
    )

    shiny::testServer(
        ascseurat::mod_visualization_server,
        args = list(
            seurat_obj = shiny::reactive(obj),
            markers_de = shiny::reactive(NULL)
        ),
        {
            session$setInputs(
                marker_source = "upload",
                markers_file = list(datapath = empty_file),
                header_opt = "Yes",
                load_markers = 1,
                slot_selection = "data"
            )
            expect_false(session$isClosed())
            session$setInputs(show_expression = 1)
            expect_false(session$isClosed())
            expect_null(plot_request())
            expect_true(any(grepl("No genes can be selected", notifications, fixed = TRUE)))
        }
    )
})

test_that("clearing DE results keeps plots of uploaded genes", {
    obj <- suppressWarnings(clustered_fixture())
    genes <- head(SeuratObject::VariableFeatures(obj), 3)
    marker_file <- tempfile("viz_markers_", fileext = ".csv")
    utils::write.csv(data.frame(gene = genes), marker_file, row.names = FALSE)
    current_markers <- shiny::reactiveVal(
        data.frame(gene = genes, cluster = "0", p_val_adj = c(0.001, 0.002, 0.003))
    )
    local_mocked_bindings(
        showNotification = function(ui, ...) invisible("id"),
        .package = "ascseurat"
    )

    shiny::testServer(
        ascseurat::mod_visualization_server,
        args = list(
            seurat_obj = shiny::reactive(obj),
            markers_de = shiny::reactive({
                shiny::req(current_markers())
                current_markers()
            })
        ),
        {
            session$setInputs(
                marker_source = "upload",
                markers_file = list(datapath = marker_file),
                header_opt = "Yes",
                load_markers = 1,
                slot_selection = "data",
                show_expression = 1
            )
            expect_identical(plot_request()$genes, genes)
            current_markers(NULL)
            session$flushReact()
            expect_identical(plot_request()$genes, genes)
        }
    )
})

test_that("expression plot downloads are hidden once the plots are cleared", {
    skip_on_cran()
    skip_if_not_installed("shinytest2")
    skip_if_not_installed("chromote")
    skip_if(is.null(suppressMessages(chromote::find_chrome())), "Chrome is not available")
    root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
    skip_if_not(file.exists(file.path(root, "DESCRIPTION")), "Package sources are not available")

    obj <- suppressWarnings(clustered_fixture())
    genes <- head(SeuratObject::VariableFeatures(obj), 3)
    app_dir <- tempfile("viz_app_")
    dir.create(app_dir)
    writeLines(c(
        sprintf("pkgload::load_all(%s, quiet = TRUE)", deparse(root)),
        sprintf("obj <- readRDS(%s)", deparse(clustered_fixture_path())),
        "renamed <- obj",
        "SeuratObject::Idents(renamed) <- paste0('c', as.character(SeuratObject::Idents(obj)))",
        sprintf(
            "markers <- data.frame(gene = %s, cluster = '0', p_val_adj = c(0.001, 0.002, 0.003))",
            paste(deparse(genes), collapse = "")
        ),
        "ui <- bslib::page_fluid(shiny::actionButton('swap', 'Swap'), mod_visualization_ui('viz'))",
        "server <- function(input, output, session) {",
        "    current_obj <- shiny::reactiveVal(obj)",
        "    shiny::observeEvent(input$swap, current_obj(renamed))",
        "    mod_visualization_server('viz', seurat_obj = current_obj,",
        "                             markers_de = shiny::reactive(markers))",
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
    visible <- "$('#viz-download_per_gene').is(':visible')"

    expect_false(app$get_js(visible))
    app$set_inputs(`viz-marker_source` = "de_results")
    app$click("viz-show_expression")
    app$wait_for_js(visible, timeout = 60000)
    expect_true(app$get_js(visible))

    app$click("swap")
    app$wait_for_js(sprintf("!%s", visible), timeout = 60000)
    expect_false(app$get_js(visible))
})
