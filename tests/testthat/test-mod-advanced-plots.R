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
