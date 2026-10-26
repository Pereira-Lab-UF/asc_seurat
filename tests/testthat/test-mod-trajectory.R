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

    # One root node: the node most start-cluster cells project to.
    roots <- res$cds@principal_graph_aux[["UMAP"]]$root_pr_nodes
    expect_length(roots, 1)
    # Recompute each start-cluster cell's nearest principal point from the
    # UMAP coordinates, without Monocle's closest-vertex index.
    principal_points <- res$cds@principal_graph_aux[["UMAP"]]$dp_mst
    umap <- SingleCellExperiment::reducedDims(res$cds)[["UMAP"]][res$root_cells, , drop = FALSE]
    nearest <- apply(umap, 1, function(cell) {
        colnames(principal_points)[which.min(colSums((principal_points - cell)^2))]
    })
    expect_equal(unname(roots), names(which.max(table(nearest))))
    expect_equal(res$root_pr_nodes, unname(roots))
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

test_that("trajectory downloads stop with a message when no result exists", {
    shiny::testServer(function(input, output, session) {}, {
        expect_error(
            ascseurat:::download_content_or_stop(function() NULL, "the trajectory gene analysis"),
            "Nothing to download yet: run the trajectory gene analysis first"
        )
        expect_error(
            ascseurat:::download_content_or_stop(function() shiny::req(FALSE), "trajectory inference"),
            "Nothing to download yet: run trajectory inference first"
        )
        table <- data.frame(gene = "A")
        expect_identical(
            ascseurat:::download_content_or_stop(function() table, "x"),
            table
        )
    })
})

slingshot_fixture <- local({
    cached <- NULL
    function() {
        if (is.null(cached)) {
            cached <<- ascseurat:::run_slingshot_trajectory(clustered_fixture())
        }
        cached
    }
})

test_that("Slingshot per-lineage pseudotime is NA for cells off that lineage", {
    skip_if_not_installed("DelayedMatrixStats")

    res <- slingshot_fixture()
    for (i in seq_along(res$assignment$path_names)) {
        values <- res$obj[[paste0("pseudotime_", res$assignment$path_names[[i]])]][, 1]
        on_lineage <- res$assignment$weights[, i] > 0
        expect_true(any(!on_lineage))
        expect_true(all(is.na(values[!on_lineage])))
        expect_true(all(is.finite(values[on_lineage])))
    }
    # tradeSeq still receives the complete pseudotime matrix.
    expect_false(anyNA(res$assignment$pseudotime))
})

test_that("PseudotimeDE-fast tests each lineage on its own cells", {
    skip_if_not_installed("PseudotimeDE")

    res <- slingshot_fixture()
    genes <- head(Seurat::VariableFeatures(res$obj), 20)
    de <- ascseurat:::run_pseudotimede_fast(res$sce, gene_subset = genes, cores = 1)
    table <- de$table

    expect_true("lineage" %in% colnames(table))
    expect_setequal(unique(table$lineage), res$assignment$path_names)
    expect_equal(nrow(table), length(genes) * length(res$assignment$path_names))
    expect_equal(table$adjusted_p_value, stats::p.adjust(table$p_value, method = "BH"))
    # Ties (many p-values are exactly 0) are broken by effect size, not by name.
    expect_true(sum(table$adjusted_p_value == 0) > 1)
    expect_equal(
        order(table$adjusted_p_value, -table$taustar, table$gene, table$lineage),
        seq_len(nrow(table))
    )

    path <- res$assignment$path_names[[2]]
    on_lineage <- res$assignment$weights[, path] > 0
    # Two genes: runTauStarDE drops a one-gene matrix to a vector.
    direct <- PseudotimeDE::runTauStarDE(
        gene.vec = genes[1:2],
        ori.tbl = data.frame(
            cell = colnames(res$sce)[on_lineage],
            pseudotime = res$assignment$pseudotime[on_lineage, path]
        ),
        mat = SingleCellExperiment::counts(res$sce),
        mode = "asymptotic",
        mc.cores = 1
    )
    row <- table[table$gene == genes[[1]] & table$lineage == path, ]
    expect_true(is.finite(row$taustar))
    expect_equal(row$taustar, direct$taustar[direct$gene == genes[[1]]])
})

test_that("tradeSeq is reproducible, adjusts p-values and uses full library sizes", {
    skip_if_not_installed("tradeSeq")

    res <- slingshot_fixture()
    genes <- head(Seurat::VariableFeatures(res$obj), 12)
    run_with_workers <- function(workers) {
        testthat::with_mocked_bindings(
            ascseurat:::run_tradeseq(res$sce, res$assignment,
                                     gene_subset = genes, nknots = 4),
            multicoreWorkers = function(...) workers,
            .package = "BiocParallel"
        )
    }

    set.seed(42)
    seed_before <- .Random.seed
    serial <- run_with_workers(1L)
    expect_identical(.Random.seed, seed_before)
    parallel <- run_with_workers(2L)
    expect_identical(parallel$table, serial$table)
    expect_identical(run_with_workers(1L)$table, serial$table)

    table <- serial$table
    expect_equal(table$adjusted_p_value, stats::p.adjust(table$p_value, method = "BH"))
    expect_true(sum(table$adjusted_p_value == 0) > 1)
    expect_equal(
        order(table$adjusted_p_value, -table$waldStat, table$gene),
        seq_len(nrow(table))
    )

    design <- SummarizedExperiment::colData(serial$fit)$tradeSeq$dm
    counts <- SingleCellExperiment::counts(res$sce)[, res$assignment$keep_cells]
    expect_equal(unname(design[, "offset(offset)"]),
                 unname(log(Matrix::colSums(counts))))
})

test_that("scMaSigPro shows internal failures instead of no significant genes", {
    skip_if_not_installed("scMaSigPro")

    res <- slingshot_fixture()
    genes <- head(Seurat::VariableFeatures(res$obj), 20)

    expect_error(
        testthat::with_mocked_bindings(
            ascseurat:::run_scmasigpro(res$sce, gene_subset = genes),
            sc.t.fit = function(...) stop("simulated scMaSigPro failure"),
            .package = "scMaSigPro"
        ),
        "simulated scMaSigPro failure"
    )

    # sc.p.vector returns the object unchanged when no gene is significant.
    none <- testthat::with_mocked_bindings(
        ascseurat:::run_scmasigpro(res$sce, gene_subset = genes),
        sc.p.vector = function(scmpObj, ...) scmpObj,
        .package = "scMaSigPro"
    )
    expect_equal(nrow(none$table), 0)
    expect_match(none$message, "no significant branching genes", fixed = TRUE)
})

test_that("downstream results only belong to the current trajectory run", {
    expect_null(ascseurat:::current_run_result(NULL, list(run_id = 1L)))
    expect_null(ascseurat:::current_run_result(list(run_id = 1L), NULL))
    expect_null(ascseurat:::current_run_result(list(run_id = 1L), list(run_id = 2L)))
    result <- list(run_id = 2L, table = data.frame(gene = "A"))
    expect_identical(ascseurat:::current_run_result(result, list(run_id = 2L)), result)
})

test_that("trajectory gene results go stale after re-running or loading new data", {
    skip_if_not_installed("DelayedMatrixStats")
    skip_if_not_installed("PseudotimeDE")

    obj <- clustered_fixture()
    Seurat::VariableFeatures(obj) <- head(Seurat::VariableFeatures(obj), 20)
    first_path <- tempfile("ascseurat_traj_a_", fileext = ".rds")
    second_path <- tempfile("ascseurat_traj_b_", fileext = ".rds")
    saveRDS(obj, first_path)
    saveRDS(obj, second_path)

    shiny::testServer(
        ascseurat::mod_trajectory_server,
        {
            session$setInputs(
                rds_source = "upload",
                rds_file = list(datapath = first_path),
                trajectory_method = "slingshot",
                set_start = 1,
                run_ti = 1
            )
            expect_false(is.null(ti_result()))

            # Root cells at pseudotime 0 are coloured, not drawn as missing.
            built <- ggplot2::ggplot_build(pseudotime_plot())
            expect_false(any(built$data[[1]]$colour == "lightgray"))

            session$setInputs(traj_de_method = "pseudotimede", n_top_genes = 10,
                              run_traj_de = 1)
            expect_false(is.null(traj_de_result()))
            expect_equal(traj_de_result()$run_id, ti_result()$run_id)

            session$setInputs(run_ti = 2)
            expect_false(is.null(ti_result()))
            expect_null(traj_de_result())
            expect_error(output$download_traj_genes, "Nothing to download yet")

            session$setInputs(rds_file = list(datapath = second_path))
            expect_null(ti_result())
            expect_error(output$download_ti_rds, "Nothing to download yet")
            expect_false(session$isClosed())
        }
    )
})

test_that("Monocle 3 gene results go stale after re-running trajectory inference", {
    skip_if_not_installed("monocle3")

    obj <- clustered_fixture()
    Seurat::VariableFeatures(obj) <- head(Seurat::VariableFeatures(obj), 20)
    rds_path <- tempfile("ascseurat_traj_monocle_", fileext = ".rds")
    saveRDS(obj, rds_path)

    shiny::testServer(
        ascseurat::mod_trajectory_server,
        {
            session$setInputs(
                rds_source = "upload",
                rds_file = list(datapath = rds_path),
                trajectory_method = "monocle3",
                set_start = 1,
                run_ti = 1
            )
            expect_false(is.null(ti_result()))
            expect_no_error(ggplot2::ggplot_build(pseudotime_plot()))

            session$setInputs(monocle_genes_to_test = "variable", monocle_top_genes = 10,
                              run_monocle_genes = 1)
            expect_false(is.null(monocle_gene_result()))

            session$setInputs(run_ti = 2)
            expect_false(is.null(ti_result()))
            expect_null(monocle_gene_result())
            expect_error(output$download_monocle_genes, "Nothing to download yet")
            expect_false(session$isClosed())
        }
    )
})
