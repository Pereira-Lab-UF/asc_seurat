test_that("paga availability helper returns a logical", {
    expect_type(ascseurat:::paga_available(), "logical")
    expect_length(ascseurat:::paga_available(), 1)
})

test_that("PAGA Python requirements include the supported Scanpy stack", {
    expect_true(all(c("scanpy", "anndata", "numpy", "scipy", "pandas", "leidenalg", "igraph") %in%
        ascseurat:::paga_python_packages()))
})

paga_test_object <- function(n_harmony = 20) {
    obj <- clustered_fixture()
    embeddings <- SeuratObject::Embeddings(obj, "pca")[, seq_len(n_harmony)] * 2
    colnames(embeddings) <- paste0("harmony_", seq_len(n_harmony))
    obj[["harmony"]] <- SeuratObject::CreateDimReducObject(
        embeddings = embeddings,
        key = "harmony_",
        assay = SeuratObject::DefaultAssay(obj)
    )
    obj
}

test_that("PAGA uses the clustering reduction and dimensions", {
    obj <- paga_test_object()
    attr(obj, "ascseurat_clustering_reduction") <- "harmony"
    attr(obj, "ascseurat_clustering_n_pcs") <- 10L
    expect_equal(
        ascseurat:::paga_reduction_dims(obj),
        list(reduction = "harmony", dims = 1:10)
    )

    attr(obj, "ascseurat_clustering_n_pcs") <- 99L
    expect_equal(ascseurat:::paga_reduction_dims(obj)$dims, 1:20)

    attr(obj, "ascseurat_clustering_reduction") <- NULL
    attr(obj, "ascseurat_clustering_n_pcs") <- NULL
    expect_equal(
        ascseurat:::paga_reduction_dims(obj),
        list(reduction = "harmony", dims = 1:20)
    )

    attr(obj, "ascseurat_clustering_reduction") <- "missing_reduction"
    expect_equal(ascseurat:::paga_reduction_dims(obj)$reduction, "harmony")

    plain <- clustered_fixture()
    attr(plain, "ascseurat_clustering_reduction") <- NULL
    expect_equal(ascseurat:::paga_reduction_dims(plain)$reduction, "pca")
})

test_that("PAGA passes the clustering embedding to AnnData", {
    obj <- paga_test_object()
    attr(obj, "ascseurat_clustering_reduction") <- "harmony"
    attr(obj, "ascseurat_clustering_n_pcs") <- 10L

    # Python is mocked: stop at AnnData creation and capture its X matrix.
    local_mocked_bindings(
        check_trajectory_dependencies = function(...) invisible(TRUE),
        .package = "ascseurat"
    )
    local_mocked_bindings(
        import = function(module, ...) {
            list(
                DataFrame = function(...) NULL,
                AnnData = function(X, ...) {
                    stop(structure(
                        class = c("anndata_x", "error", "condition"),
                        list(message = "captured", call = NULL, X = X)
                    ))
                }
            )
        },
        r_to_py = function(x, ...) x,
        dict = function(...) list(...),
        .package = "reticulate"
    )

    x <- tryCatch(
        ascseurat:::run_paga_trajectory(obj),
        anndata_x = function(e) e$X
    )
    expected <- SeuratObject::Embeddings(obj, "harmony")[, 1:10]
    expect_equal(x, unname(expected))
})

test_that("PAGA top-edge gene tables are ranked across all edges", {
    obj <- clustered_fixture()
    obj$pseudotime <- SeuratObject::Embeddings(obj, "pca")[, 1]
    obj$paga_cluster <- as.character(SeuratObject::Idents(obj))
    clusters <- ascseurat:::cluster_levels(obj)
    connectivities <- matrix(
        0,
        nrow = length(clusters),
        ncol = length(clusters),
        dimnames = list(clusters, clusters)
    )
    connectivities[1, 2] <- connectivities[2, 1] <- 0.9
    connectivities[1, 4] <- connectivities[4, 1] <- 0.5
    res <- list(
        obj = obj,
        method = "paga",
        paga = list(connectivities = connectivities)
    )
    edges <- ascseurat:::selected_paga_edges(res, scope = "top", top_n = 2)
    genes <- head(Seurat::VariableFeatures(obj), 30)

    markers <- ascseurat:::run_paga_edge_markers(res, edges, gene_subset = genes)$table
    expect_setequal(unique(markers$edge), edges$edge)
    expect_equal(
        order(markers$p_val_adj, -abs(markers$avg_log2FC), markers$gene),
        seq_len(nrow(markers))
    )

    pseudotime_genes <- ascseurat:::run_paga_pseudotime_genes(
        res, edges, gene_subset = genes
    )$table
    expect_setequal(unique(pseudotime_genes$edge), edges$edge)
    expect_equal(
        order(pseudotime_genes$p_val_adj, -abs(pseudotime_genes$rho),
              pseudotime_genes$gene),
        seq_len(nrow(pseudotime_genes))
    )
})
