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
