test_that("trajectory UI exposes supported methods and recommended start option", {
    html <- htmltools::renderTags(ascseurat::mod_trajectory_ui("traj"))$html

    expect_match(html, "Slingshot", fixed = TRUE)
    expect_match(html, "PAGA", fixed = TRUE)
    expect_false(grepl("PAGA \\(Scanpy\\)", html))
    expect_match(html, "Monocle 3", fixed = TRUE)
    expect_match(html, "Yes \\(strongly recommended\\)")
    expect_match(html, "PseudotimeDE-fast \\(recommended\\)")
    expect_match(html, "tradeSeq \\(long runtime\\)")
    expect_match(html, "It can take a long time on large datasets", fixed = TRUE)
    expect_match(html, "scMaSigPro identifies genes that change", fixed = TRUE)
    expect_match(html, "Top genes to display/plot", fixed = TRUE)
    expect_match(html, "All variable genes are tested", fixed = TRUE)
    expect_match(html, "Connected Cluster Gene Discovery", fixed = TRUE)
    expect_match(html, "Connected cluster markers", fixed = TRUE)
    expect_match(html, "Pseudotime-associated genes", fixed = TRUE)
    expect_match(html, "Trajectory-Variable Genes", fixed = TRUE)
    expect_match(html, "Find Gene Modules", fixed = TRUE)
    expect_match(html, "Gene Module Heatmap", fixed = TRUE)
    expect_false(grepl("Tests expression changes along Slingshot pseudotime", html,
                       fixed = TRUE))
    expect_false(grepl("For this method, use the saved trajectory object", html,
                       fixed = TRUE))
    expect_match(html, "Download Trajectory Object", fixed = TRUE)
    expect_false(grepl("Download Trajectory Object (RDS)", html, fixed = TRUE))
    expect_match(html, "Set root/start cluster", fixed = TRUE)
    expect_match(html, "value=\"upload\" checked")
    expect_match(html, "value=\"pseudotimede\" selected")
    expect_match(html, "value=\"1\" checked")
    expect_match(html, "traj-end_cluster_ui", fixed = TRUE)
    expect_match(html, "input[&#39;traj-trajectory_method&#39;] == &#39;slingshot&#39;",
                 fixed = TRUE)
    expect_lt(regexpr("Upload from browser", html, fixed = TRUE)[[1]],
              regexpr("Path on this computer/server", html, fixed = TRUE)[[1]])
    expect_lt(regexpr("Yes \\(strongly recommended\\)", html)[[1]],
              regexpr(">No<", html, fixed = TRUE)[[1]])
    expect_lt(regexpr("scMaSigPro identifies genes that change", html,
                      fixed = TRUE)[[1]],
              regexpr("Top genes to display/plot", html, fixed = TRUE)[[1]])
})

test_that("trajectory controls wait for uploaded data without missing-path errors", {
    expect_false(ascseurat:::trajectory_rds_available(
        "path",
        "RDS_files/seurat_clustered.rds",
        NULL
    ))
    expect_false(ascseurat:::trajectory_rds_available("upload", "", NULL))

    shiny::testServer(
        ascseurat::mod_trajectory_server,
        {
            session$setInputs(
                rds_source = "upload",
                trajectory_method = "paga",
                set_start = 1
            )

            expect_match(output$start_cluster_ui$html,
                         "Load a processed Seurat RDS", fixed = TRUE)
            session$setInputs(run_ti = 1)
            expect_null(ti_result())
        }
    )
})

test_that("trajectory module produces pseudotime on a clustered object", {
    skip_if_not_installed("DelayedMatrixStats")

    shiny::testServer(
        ascseurat::mod_trajectory_server,
        {
            session$setInputs(
                rds_source = "upload",
                rds_file = list(datapath = clustered_fixture_path()),
                trajectory_method = "slingshot",
                set_start = 1,
                run_ti = 1
            )
            res <- ti_result()
            expect_true("pseudotime" %in% colnames(res$obj@meta.data))
            expect_true("trajectory_method" %in% colnames(res$obj@meta.data))
            expect_equal(unique(res$obj$trajectory_method), "Slingshot")
            mst_layout <- ascseurat:::slingshot_mst_layout(res)
            expect_true(nrow(mst_layout$nodes) > 0)
            expect_true(nrow(mst_layout$segments) > 0)
            expect_true(nrow(mst_layout$lineage_paths) > 0)
            expect_s3_class(ascseurat:::plot_slingshot_mst(res), "ggplot")
        }
    )
})

test_that("PAGA trajectory smoke test produces standardized metadata", {
    skip_if_not(ascseurat:::paga_available(), "PAGA Python stack is not available")

    res <- ascseurat:::run_paga_trajectory(clustered_fixture())
    expect_equal(res$method, "paga")
    expect_true(all(c("trajectory_method", "pseudotime", "trajectory_path") %in%
        colnames(res$obj@meta.data)))
    expect_true(any(is.finite(res$obj$pseudotime)))
})

test_that("Monocle 3 trajectory smoke test produces standardized metadata", {
    skip_if_not_installed("monocle3")

    res <- ascseurat:::run_monocle3_trajectory(clustered_fixture())
    expect_equal(res$method, "monocle3")
    expect_true(all(c("trajectory_method", "pseudotime", "trajectory_path") %in%
        colnames(res$obj@meta.data)))
})

test_that("PAGA strongest connection summary is compact and sorted", {
    connectivities <- matrix(
        c(
            0, 0.1, 0.4,
            0.1, 0, 0.02,
            0.4, 0.02, 0
        ),
        nrow = 3,
        dimnames = list(c("0", "1", "2"), c("0", "1", "2"))
    )

    edges <- ascseurat:::strongest_paga_connections(connectivities, n = 2)
    expect_equal(edges$from, c("0", "0"))
    expect_equal(edges$to, c("2", "1"))
    expect_equal(edges$connectivity, c(0.4, 0.1))

    all_edges <- ascseurat:::paga_connections_table(connectivities, min_weight = 0.05)
    expect_equal(all_edges$edge, c("0||2", "0||1"))
    expect_equal(all_edges$source_cluster, c("0", "0"))
    expect_equal(all_edges$target_cluster, c("2", "1"))
    expect_equal(all_edges$connectivity, c(0.4, 0.1))
})

test_that("PAGA connected-cluster marker mode returns unified gene table", {
    obj <- clustered_fixture()
    clusters <- ascseurat:::cluster_levels(obj)
    connectivities <- matrix(
        0,
        nrow = length(clusters),
        ncol = length(clusters),
        dimnames = list(clusters, clusters)
    )
    connectivities[1, 2] <- connectivities[2, 1] <- 0.9
    res <- list(
        obj = obj,
        method = "paga",
        paga = list(connectivities = connectivities)
    )

    edges <- ascseurat:::selected_paga_edges(res, scope = "top", top_n = 1)
    result <- ascseurat:::run_paga_edge_markers(
        res,
        edges,
        gene_subset = head(rownames(obj), 10)
    )

    expect_equal(result$mode, "edge_markers")
    expect_true(all(c(
        "edge", "source_cluster", "target_cluster", "gene", "direction",
        "avg_log2FC", "p_val_adj", "pct_source", "pct_target", "connectivity"
    ) %in% colnames(result$table)))
})

test_that("PAGA pseudotime gene mode handles uninformative pseudotime", {
    obj <- clustered_fixture()
    obj$pseudotime <- 1
    obj$paga_cluster <- as.character(SeuratObject::Idents(obj))
    clusters <- ascseurat:::cluster_levels(obj)
    connectivities <- matrix(
        0,
        nrow = length(clusters),
        ncol = length(clusters),
        dimnames = list(clusters, clusters)
    )
    connectivities[1, 2] <- connectivities[2, 1] <- 0.9
    res <- list(
        obj = obj,
        method = "paga",
        paga = list(connectivities = connectivities)
    )

    edges <- ascseurat:::selected_paga_edges(res, scope = "top", top_n = 1)
    result <- ascseurat:::run_paga_pseudotime_genes(
        res,
        edges,
        gene_subset = head(rownames(obj), 10)
    )

    expect_equal(result$mode, "pseudotime")
    expect_equal(nrow(result$table), 0)
    expect_true(all(c("gene", "rho", "p_val_adj") %in% colnames(result$table)))
})

test_that("Monocle 3 graph-test wrapper reports missing dependency clearly", {
    if (requireNamespace("monocle3", quietly = TRUE)) {
        skip("monocle3 is installed; trajectory smoke tests cover the dependency path")
    }

    expect_error(
        ascseurat:::run_monocle3_graph_test(list(cds = NULL)),
        "Monocle 3 gene discovery requires the monocle3 package",
        fixed = TRUE
    )
})

test_that("Monocle 3 module heatmap plots module genes across pseudotime", {
    obj <- clustered_fixture()
    obj$pseudotime <- seq_len(ncol(obj))
    genes <- head(rownames(obj), 6)
    module_result <- list(
        modules = data.frame(
            gene_id = genes,
            module = rep(c("1", "2"), each = 3),
            stringsAsFactors = FALSE
        )
    )
    res <- list(obj = obj)

    plot <- ascseurat:::plot_monocle_module_heatmap(res, module_result,
                                                    genes_per_module = 2,
                                                    bins = 4)
    expect_s3_class(plot, "ggplot")
})

test_that("PseudotimeDE-fast smoke test ranks Slingshot trajectory genes when installed", {
    skip_if_not_installed("PseudotimeDE")

    res <- ascseurat:::run_slingshot_trajectory(clustered_fixture())
    de <- ascseurat:::run_pseudotimede_fast(
        res$sce,
        gene_subset = head(rownames(res$obj), 10),
        cores = 1
    )

    expect_equal(de$engine, "PseudotimeDE-fast")
    expect_true(all(c("gene", "p_value", "adjusted_p_value") %in%
        colnames(de$table)))
})

test_that("PseudotimeDE-fast missing dependency message is friendly", {
    msg <- ascseurat:::pseudotimede_install_message()

    expect_match(msg, "PseudotimeDE-fast is not installed", fixed = TRUE)
    expect_match(msg, "remotes::install_git", fixed = TRUE)
    expect_false(grepl("Package 'PseudotimeDE'", msg, fixed = TRUE))

    if (!ascseurat:::pseudotimede_available()) {
        error_message <- tryCatch(
            ascseurat:::assert_pseudotimede_available(),
            error = conditionMessage
        )
        expect_match(error_message, "PseudotimeDE-fast is not installed",
                     fixed = TRUE)
        expect_false(grepl("Package 'PseudotimeDE'", error_message, fixed = TRUE))
    }
})
