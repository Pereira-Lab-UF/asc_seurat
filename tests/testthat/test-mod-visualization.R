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
