test_that("clustering module returns a clustered Seurat object", {
    obj <- ascseurat:::run_normalization(demo_fixture(), method = "lognorm", nfeatures = 500)

    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(obj),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_clustering = 1)
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
        }
    )
})

test_that("a failed subset re-integration keeps the previous clustering", {
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
            clustered <- active_obj()

            session$setInputs(
                filter_clusters = 1,
                reanalyze_strategy = "recompute",
                filter_mode = "select",
                selected_clusters = c("0", "1"),
                reanalyze = 1
            )
            expect_false(session$isClosed())
            expect_true(any(grepl(
                "Subset reanalysis failed: re-integration with harmony failed (too few cells",
                notes(),
                fixed = TRUE
            )))
            expect_identical(active_obj(), clustered)
            expect_null(pending_obj())
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
