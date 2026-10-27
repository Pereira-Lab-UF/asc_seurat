test_that("clustering module returns a clustered Seurat object", {
    obj <- ascseurat:::run_normalization(demo_fixture(), method = "lognorm", nfeatures = 500)

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_tsne = TRUE, run_clustering = 1)
            expect_s4_class(active_obj(), "Seurat")
            expect_true("tsne" %in% names(active_obj()@reductions))
            expect_equal(attr(active_obj(), "ascseurat_clustering_n_pcs", exact = TRUE), 10L)

            session$setInputs(rename_1 = "Renamed cluster", apply_rename = 1)
            expect_true(any(grepl("Renamed cluster", levels(SeuratObject::Idents(active_obj())))))
        }
    )
})

capture_notifications <- function(env = parent.frame()) {
    notes <- character()
    testthat::local_mocked_bindings(
        showNotification = function(ui, ...) {
            notes <<- c(notes, paste(as.character(ui), collapse = " "))
            invisible(NULL)
        },
        .env = env
    )
    function() notes
}

two_sample_fixture <- local({
    cached <- NULL
    function() {
        if (is.null(cached)) {
            obj <- demo_fixture()
            obj$samples <- rep(c("A", "B"), length.out = ncol(obj))
            obj[["RNA"]] <- split(obj[["RNA"]], f = obj$samples)
            cached <<- suppressWarnings(ascseurat:::run_integration(
                obj,
                method = "harmony",
                dims = 1:10,
                sample_col = "samples",
                nfeatures = 500
            ))
        }
        cached
    }
})

test_that("reanalysis refuses empty and tiny subsets with a notification", {
    obj <- clustered_fixture()
    notes <- capture_notifications()

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            state_obj(obj)
            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "exclude",
                selected_clusters = levels(SeuratObject::Idents(obj)),
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_true(any(grepl("No cells remain", notes())))
            expect_equal(ncol(active_obj()), ncol(obj))

            tiny <- obj
            SeuratObject::Idents(tiny) <- factor(rep(c("tiny", "rest"), c(5, ncol(obj) - 5)))
            state_obj(tiny)
            session$setInputs(filter_mode = "select", selected_clusters = "tiny", reanalyze = 2)
            expect_false(session$isClosed())
            expect_true(any(grepl("only 5 cells", notes())))
            expect_null(pending_obj())
        }
    )
})

test_that("small subsets can be reanalyzed and clustered with tSNE", {
    obj <- clustered_fixture()
    SeuratObject::Idents(obj) <- factor(rep(c("small", "rest"), c(40, ncol(obj) - 40)))

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            state_obj(obj)
            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = "small",
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_equal(ncol(pending_obj()), 40)

            session$setInputs(resolution = 0.4, run_tsne = TRUE, run_clustering = 1)
            expect_false(session$isClosed())
            expect_equal(ncol(active_obj()), 40)
            expect_true(all(c("umap", "tsne") %in% names(active_obj()@reductions)))
        }
    )
})

test_that("clustering errors show a notification instead of closing the session", {
    obj <- clustered_fixture()
    notes <- capture_notifications()
    local_mocked_bindings(run_clustering = function(...) stop("boom"))

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_clustering = 1)
            expect_false(session$isClosed())
            expect_true(any(grepl("Clustering failed: boom", notes())))
            expect_null(state_obj())
        }
    )
})

test_that("clustering state resets when the upstream object becomes unavailable", {
    obj <- clustered_fixture()
    upstream <- shiny::reactiveVal(obj)

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive({
                shiny::req(upstream())
                upstream()
            }),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_tsne = FALSE, run_clustering = 1)
            expect_false(is.null(state_obj()))

            upstream(NULL)
            session$flushReact()
            expect_null(state_obj())
            expect_null(pending_obj())
            expect_false(pending_reanalysis())
        }
    )
})

test_that("later clustering runs stay on the reanalyzed subset until the upstream changes", {
    obj <- clustered_fixture()
    upstream <- shiny::reactiveVal(obj)

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive({
                shiny::req(upstream())
                upstream()
            }),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_tsne = FALSE, run_clustering = 1)
            n_full <- ncol(active_obj())

            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "keep",
                filter_mode = "exclude",
                selected_clusters = "0",
                reanalyze = 1
            )
            n_kept <- ncol(active_obj())
            expect_lt(n_kept, n_full)
            session$setInputs(run_clustering = 2)
            expect_equal(ncol(active_obj()), n_kept)

            session$setInputs(
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = c("0", "1"),
                reanalyze = 2
            )
            n_subset <- ncol(pending_obj())
            session$setInputs(run_clustering = 3)
            expect_equal(ncol(active_obj()), n_subset)
            session$setInputs(resolution = 0.8, run_clustering = 4)
            expect_equal(ncol(active_obj()), n_subset)

            upstream(subset(obj, cells = colnames(obj)[1:500]))
            session$flushReact()
            session$setInputs(run_clustering = 5)
            expect_equal(ncol(active_obj()), 500)
        }
    )
})

test_that("recomputing a subset of an integrated object repeats the batch correction", {
    int_obj <- two_sample_fixture()
    notes <- capture_notifications()

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(int_obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_clustering = 1)
            expect_equal(attr(active_obj(), "ascseurat_clustering_reduction"), "harmony")

            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = c("0", "1"),
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_true("harmony" %in% names(pending_obj()@reductions))
            n_subset <- ncol(pending_obj())

            session$setInputs(run_clustering = 2)
            expect_equal(attr(active_obj(), "ascseurat_clustering_reduction"), "harmony")
            expect_equal(ncol(active_obj()), n_subset)

            # A subset with a single sample falls back to plain normalization.
            one_sample <- active_obj()
            SeuratObject::Idents(one_sample) <- factor(one_sample$samples)
            state_obj(one_sample)
            session$setInputs(selected_clusters = "A", reanalyze = 2)
            expect_false(session$isClosed())
            expect_identical(names(pending_obj()@reductions), "pca")
            expect_true(any(grepl("Only one sample remains", notes())))
            # JoinLayers() alone keeps counts.A/data.A next to the joined layers.
            rna_layers <- SeuratObject::Layers(pending_obj()[["RNA"]])
            expect_true(all(c("counts", "data") %in% rna_layers))
            expect_false(any(grepl("^(counts|data)\\.", rna_layers)))
        }
    )
})

test_that("a failed subset re-integration falls back to reprocessing without batch correction", {
    int_obj <- two_sample_fixture()
    notes <- capture_notifications()
    local_mocked_bindings(run_integration = function(...) {
        stop("too few cells per sample for anchors")
    })

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(int_obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_clustering = 1)

            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = c("0", "1"),
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_false(any(grepl("Subset reanalysis failed", notes())))
            expect_true(any(grepl(
                paste0(
                    "Batch correction with harmony could not be repeated for this subset ",
                    "(too few cells per sample for anchors), so it was reprocessed without batch correction."
                ),
                notes(),
                fixed = TRUE
            )))
            expect_identical(names(pending_obj()@reductions), "pca")
            expect_null(pending_obj()@misc$ascseurat_integration)
        }
    )
})

test_that("a subset re-integration with missing values falls back to reprocessing without batch correction", {
    int_obj <- two_sample_fixture()
    notes <- capture_notifications()
    local_mocked_bindings(run_integration = function(obj, ...) {
        obj <- suppressWarnings(ascseurat:::run_normalization(obj, nfeatures = 500))
        embeddings <- SeuratObject::Embeddings(obj, "pca")
        embeddings[1:5, ] <- NA
        obj[["integrated.rpca"]] <- SeuratObject::CreateDimReducObject(
            embeddings = embeddings,
            key = "integratedrpca_",
            assay = "RNA"
        )
        attr(obj, "ascseurat_default_reduction") <- "integrated.rpca"
        obj
    })

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(int_obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_clustering = 1)

            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = c("0", "1"),
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_true(any(grepl(
                "could not be repeated for this subset (the integrated reduction has missing values)",
                notes(),
                fixed = TRUE
            )))
            expect_identical(names(pending_obj()@reductions), "pca")

            session$setInputs(run_clustering = 2)
            expect_false(any(grepl("Clustering failed", notes())))
            expect_equal(attr(active_obj(), "ascseurat_clustering_reduction"), "pca")
        }
    )
})

rpca_fixture <- local({
    cached <- NULL
    function() {
        if (is.null(cached)) {
            obj <- demo_fixture()
            obj$samples <- rep(c("A", "B"), length.out = ncol(obj))
            obj[["RNA"]] <- split(obj[["RNA"]], f = obj$samples)
            cached <<- suppressWarnings(ascseurat:::run_integration(
                obj,
                method = "rpca",
                dims = 1:10,
                sample_col = "samples",
                nfeatures = 500
            ))
        }
        cached
    }
})

test_that("an RPCA subset with too few cells in one sample is reprocessed without batch correction", {
    int_obj <- rpca_fixture()
    picked <- c(
        colnames(int_obj)[int_obj$samples == "A"][1:152],
        colnames(int_obj)[int_obj$samples == "B"][1:8]
    )
    SeuratObject::Idents(int_obj) <- factor(ifelse(colnames(int_obj) %in% picked, "picked", "rest"))
    notes <- capture_notifications()

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(int_obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            state_obj(int_obj)
            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = "picked",
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_false(any(grepl("Subset reanalysis failed", notes())))
            expect_equal(ncol(pending_obj()), 160)
            expect_identical(names(pending_obj()@reductions), "pca")
            expect_true(any(grepl(
                paste0(
                    "Batch correction with rpca could not be repeated for this subset ",
                    "(fewer than two samples have at least 31 cells)"
                ),
                notes(),
                fixed = TRUE
            )))

            session$setInputs(resolution = 0.4, run_clustering = 1)
            expect_equal(attr(active_obj(), "ascseurat_clustering_reduction"), "pca")
        }
    )
})

test_that("an RPCA subset with about 30 cells per sample can be clustered", {
    int_obj <- rpca_fixture()
    # RPCA on two samples of 31 cells returns rows of NA in integrated.rpca.
    picked <- c(
        colnames(int_obj)[int_obj$samples == "A"][1:31],
        colnames(int_obj)[int_obj$samples == "B"][1:31]
    )
    SeuratObject::Idents(int_obj) <- factor(ifelse(colnames(int_obj) %in% picked, "picked", "rest"))
    notes <- capture_notifications()

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(int_obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            state_obj(int_obj)
            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = "picked",
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_equal(ncol(pending_obj()), 62)

            session$setInputs(resolution = 0.4, run_clustering = 1)
            expect_false(any(grepl("Clustering failed", notes())))
            expect_identical(output$clustering_ready, "true")
            reduction <- attr(active_obj(), "ascseurat_clustering_reduction")
            expect_false(anyNA(SeuratObject::Embeddings(active_obj(), reduction)))
        }
    )
})

test_that("an SCTransform subset with a tiny sample falls back on joined RNA layers", {
    obj <- demo_fixture()[, 1:300]
    obj$samples <- rep(c("A", "B"), c(290, 10))
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$samples)
    int_obj <- suppressWarnings(ascseurat:::run_integration(
        obj,
        method = "harmony",
        dims = 1:10,
        sample_col = "samples",
        normalization_method = "sctransform",
        nfeatures = 500
    ))
    # SCTransform fails on a 2-cell layer, so RNA must not stay split by sample.
    picked <- c(colnames(int_obj)[1:100], colnames(int_obj)[291:292])
    SeuratObject::Idents(int_obj) <- factor(ifelse(colnames(int_obj) %in% picked, "picked", "rest"))
    notes <- capture_notifications()

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(int_obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            state_obj(int_obj)
            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = "picked",
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_false(any(grepl("Subset reanalysis failed", notes())))
            expect_true(any(grepl("fewer than two samples have at least 10 cells", notes(), fixed = TRUE)))
            expect_equal(ncol(pending_obj()), 102)
            expect_equal(SeuratObject::DefaultAssay(pending_obj()), "SCT")
        }
    )
})

test_that("RPCA subset re-integration groups small samples and caps the dims", {
    int_obj <- rpca_fixture()
    # Integrated with 40 PCs: more than the smallest usable sample has cells.
    int_obj@misc$ascseurat_integration$dims <- 1:40
    b_cells <- colnames(int_obj)[int_obj$samples == "B"]
    small_cells <- b_cells[36:43]
    int_obj@meta.data[small_cells, "samples"] <- "C"
    picked <- c(colnames(int_obj)[int_obj$samples == "A"][1:100], b_cells[1:35], small_cells)
    SeuratObject::Idents(int_obj) <- factor(ifelse(colnames(int_obj) %in% picked, "picked", "rest"))
    notes <- capture_notifications()

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(int_obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            state_obj(int_obj)
            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = "picked",
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_false(any(grepl("Subset reanalysis failed|could not be repeated", notes())))
            sub <- pending_obj()
            expect_true("integrated.rpca" %in% names(sub@reductions))
            expect_true(any(grepl(
                "Samples with fewer than 31 cells in the subset (C) were not batch-corrected separately",
                notes(),
                fixed = TRUE
            )))
            # Sample labels and the original settings are kept for later subsets.
            expect_equal(as.integer(table(sub$samples)[c("A", "B", "C")]), c(100L, 35L, 8L))
            expect_false("ascseurat_batch" %in% colnames(sub@meta.data))
            expect_identical(sub@misc$ascseurat_integration, int_obj@misc$ascseurat_integration)
        }
    )
})

test_that("recomputing a subset repeats the recorded normalization settings", {
    obj <- demo_fixture()
    obj$dropped_later <- obj$nCount_RNA
    obj <- suppressWarnings(ascseurat:::run_normalization(
        obj,
        method = "lognorm",
        nfeatures = 400,
        vars_to_regress = c("percent.mt", "dropped_later"),
        score_cell_cycle = TRUE,
        regress_cc = TRUE
    ))
    # Variables missing from the subset metadata are skipped.
    obj$dropped_later <- NULL
    obj <- ascseurat:::run_clustering(obj, dims = 1:10, resolution = 0.4)
    notes <- capture_notifications()

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            state_obj(obj)
            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = c("0", "1"),
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_false(any(grepl("Subset reanalysis failed", notes())))
            sub <- pending_obj()
            expect_length(SeuratObject::VariableFeatures(sub), 400)
            expect_identical(
                sub@misc$ascseurat_scale_regress,
                c("percent.mt", "S.Score", "G2M.Score")
            )
            expect_identical(sub@misc$ascseurat_scale_cells, as.integer(ncol(sub)))
            # The PC setting is in Step 1 of the Integration tab, not Step 3.
            expect_true(any(grepl("using the selected number of PCs", notes(), fixed = TRUE)))
            expect_false(any(grepl("Step 3", notes(), fixed = TRUE)))
        }
    )
})

test_that("subset Recompute skips regressed variables with a single value in the subset", {
    obj <- demo_fixture()
    obj$batch <- rep(c("b1", "b2"), length.out = ncol(obj))
    obj <- suppressWarnings(ascseurat:::run_normalization(
        obj,
        method = "lognorm",
        nfeatures = 500,
        vars_to_regress = c("batch", "percent.mt")
    ))
    SeuratObject::Idents(obj) <- factor(obj$batch)
    notes <- capture_notifications()

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            state_obj(obj)
            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = "b1",
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_false(any(grepl("Subset reanalysis failed", notes())))
            sub <- pending_obj()
            expect_equal(ncol(sub), sum(obj$batch == "b1"))
            expect_identical(sub@misc$ascseurat_scale_regress, "percent.mt")
            expect_true(any(grepl(
                "Not regressed in the subset because it has a single value: batch.",
                notes(),
                fixed = TRUE
            )))
        }
    )
})

test_that("subsets of objects without a normalization record use the defaults", {
    obj <- clustered_fixture()
    obj@misc$ascseurat_normalization <- NULL
    SeuratObject::Idents(obj) <- factor(rep(c("small", "rest"), c(60, ncol(obj) - 60)))

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            state_obj(obj)
            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = "small",
                reanalyze = 1
            )
            expect_false(session$isClosed())
            settings <- pending_obj()@misc$ascseurat_normalization
            expect_identical(settings$method, "lognorm")
            expect_equal(settings$nfeatures, 3000)
        }
    )
})

test_that("clustering results survive upstream re-runs that return the same object", {
    obj <- clustered_fixture()
    filter_mode <- shiny::reactiveVal("none")
    upstream <- shiny::reactiveVal(obj)

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            # Like active_integrated_obj(), which reads the filter radio button.
            obj = shiny::reactive({
                filter_mode()
                upstream()
            }),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_tsne = FALSE, run_clustering = 1)
            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "keep",
                filter_mode = "exclude",
                selected_clusters = "0",
                reanalyze = 1
            )
            kept <- state_obj()
            expect_false(is.null(subset_base()))

            filter_mode("filter")
            session$flushReact()
            expect_identical(state_obj(), kept)
            expect_identical(subset_base(), kept)

            upstream(subset(obj, cells = colnames(obj)[1:500]))
            session$flushReact()
            expect_null(state_obj())
            expect_null(subset_base())
        }
    )
})

test_that("tSNE runs only when the checkbox is ticked for a single sample", {
    obj <- clustered_fixture()
    multi <- obj
    multi$samples <- rep(c("A", "B"), length.out = ncol(multi))
    upstream <- shiny::reactiveVal(obj)
    notes <- capture_notifications()

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(upstream()),
            n_pcs = shiny::reactive(10)
        )),
        {
            # Integration tab: there is no tSNE checkbox.
            session$setInputs(resolution = 0.4, run_clustering = 1)
            expect_false("tsne" %in% names(active_obj()@reductions))

            # Single Sample tab with metadata that adds several samples.
            upstream(multi)
            session$flushReact()
            session$setInputs(run_tsne = TRUE, run_clustering = 2)
            expect_false("tsne" %in% names(active_obj()@reductions))
            expect_true(any(grepl("tSNE is skipped because the object has several samples", notes())))
        }
    )
})

test_that("a re-normalized integrated object is clustered on its new PCA", {
    obj <- suppressWarnings(ascseurat:::run_normalization(
        two_sample_fixture(),
        method = "lognorm",
        nfeatures = 500
    ))

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_clustering = 1)
            expect_equal(attr(active_obj(), "ascseurat_clustering_reduction"), "pca")
        }
    )
})

test_that("recomputing a subset re-integrates objects with a v3 RNA assay", {
    withr::local_options(Seurat.object.assay.version = "v3")
    counts <- SeuratObject::LayerData(demo_fixture(), layer = "counts")
    a <- SeuratObject::CreateSeuratObject(counts[, 1:1000], project = "A")
    b <- SeuratObject::CreateSeuratObject(counts[, 1001:2000], project = "B")
    merged <- merge(a, y = b, add.cell.ids = c("A", "B"))
    merged$samples <- merged$orig.ident
    expect_false(inherits(merged[["RNA"]], "StdAssay"))
    int_obj <- suppressWarnings(ascseurat:::run_integration(
        merged,
        method = "harmony",
        dims = 1:10,
        sample_col = "samples",
        nfeatures = 500
    ))
    notes <- capture_notifications()

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(int_obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_clustering = 1)
            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = c("0", "1"),
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_false(any(grepl("Subset reanalysis failed", notes())))
            expect_true("harmony" %in% names(pending_obj()@reductions))
        }
    )
})

test_that("capping the clustering PCs to the integrated reduction is reported", {
    int_obj <- two_sample_fixture()
    notes <- capture_notifications()

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(int_obj),
            n_pcs = shiny::reactive(25)
        )),
        {
            session$setInputs(resolution = 0.4, run_clustering = 1)
            expect_equal(attr(active_obj(), "ascseurat_clustering_n_pcs", exact = TRUE), 10L)
            expect_true(any(grepl(
                "used the 10 dimensions available in the 'harmony' reduction instead of the requested 25",
                notes(),
                fixed = TRUE
            )))
        }
    )
})
