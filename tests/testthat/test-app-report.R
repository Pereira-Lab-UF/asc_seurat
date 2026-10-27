app_file <- function() {
    path <- system.file("app", "app.R", package = "ascseurat")
    testthat::skip_if(!nzchar(path), "The app directory is not available.")
    path
}

logged_actions <- function(events) {
    vapply(events, function(event) event$action, character(1))
}

write_report_10x <- function(dir, cells) {
    counts <- SeuratObject::LayerData(demo_fixture()[, cells], layer = "counts")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    Matrix::writeMM(counts, file.path(dir, "matrix.mtx"))
    writeLines(colnames(counts), file.path(dir, "barcodes.tsv"))
    writeLines(paste(rownames(counts), rownames(counts), sep = "\t"),
               file.path(dir, "genes.tsv"))
    dir
}

# Steps and their inputs, as listed under "## Steps run" in the Markdown report.
report_steps <- function(lines) {
    starts <- grep("^[0-9]+\\. \\*\\*", lines)
    ends <- c(starts[-1] - 1L, length(lines))
    steps <- lapply(seq_along(starts), function(i) {
        block <- lines[starts[[i]]:ends[[i]]]
        params <- grep("^   - ", block, value = TRUE)
        list(
            action = sub("^[0-9]+\\. \\*\\*(.*)\\*\\* \\(.*$", "\\1", block[[1]]),
            params = sub("^   - ", "", params)
        )
    })
    steps
}

test_that("a refused QC filter is not logged and an applied one is", {
    shiny::testServer(shiny::shinyAppFile(app_file()), {
        # The first flush only initialises the observers.
        session$setInputs(main_nav = "home")
        # The QC violin plots warn that the demo object has no "data" layer.
        suppressWarnings(session$setInputs(main_nav = "demo"))
        expect_s4_class(demo_data(), "Seurat")
        expect_equal(logged_actions(session_log$events), "Load demo data (demo tab)")

        session$setInputs(
            `demo_qc-filter_mode` = "filter",
            `demo_qc-min_genes` = NA,
            `demo_qc-max_genes` = 5000,
            `demo_qc-max_mito` = 5
        )
        session$setInputs(`demo_qc-apply_filter` = 1)
        expect_equal(logged_actions(session_log$events), "Load demo data (demo tab)")

        session$setInputs(`demo_qc-min_genes` = 200)
        suppressWarnings(session$setInputs(`demo_qc-apply_filter` = 2))
        expect_equal(
            logged_actions(session_log$events),
            c("Load demo data (demo tab)", "QC filter (demo)")
        )
        expect_equal(session_log$events[[2]]$params$min_genes, "200")
    })
})

test_that("integration post-filter apply and reset are logged when they change the object", {
    root <- withr::local_tempdir()
    sample_path <- write_report_10x(file.path(root, "s1"), 1:300)
    testthat::local_mocked_bindings(
        run_integration = function(obj, ...) {
            embeddings <- matrix(
                stats::rnorm(ncol(obj) * 10),
                nrow = ncol(obj),
                dimnames = list(colnames(obj), paste0("PC_", 1:10))
            )
            obj[["pca"]] <- SeuratObject::CreateDimReducObject(
                embeddings = embeddings,
                key = "PC_",
                assay = "RNA"
            )
            obj
        },
        .package = "ascseurat"
    )

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        session$setInputs(
            `multi_int-n_samples` = 1,
            `multi_int-sample_name_1` = "s1",
            `multi_int-sample_path_1` = sample_path,
            `multi_int-qc_mito_pattern_1` = "",
            `multi_int-qc_min_genes_1` = 0,
            `multi_int-qc_max_genes_1` = 1e6,
            `multi_int-qc_max_mito_1` = 100,
            `multi_int-int_method` = "harmony",
            `multi_int-norm_method` = "lognorm",
            `multi_int-n_dims` = 10,
            `multi_int-post_filter_mode` = "none"
        )
        suppressMessages(session$setInputs(`multi_int-run_integration` = 1))
        expect_s4_class(multi_int$obj(), "Seurat")
        expect_equal(logged_actions(session_log$events), "Run integration")

        # Refused: no cells pass the threshold.
        session$setInputs(
            `multi_int-post_filter_mode` = "filter",
            `multi_int-filter_min_features` = 1e9
        )
        session$setInputs(`multi_int-apply_int_filter` = 1)
        expect_equal(logged_actions(session_log$events), "Run integration")

        session$setInputs(`multi_int-filter_min_features` = 0)
        session$setInputs(`multi_int-apply_int_filter` = 2)
        session$setInputs(`multi_int-reset_int_filter` = 1)
        expect_equal(
            logged_actions(session_log$events),
            c(
                "Run integration",
                "Apply post-integration filter",
                "Reset post-integration filter"
            )
        )
        expect_equal(session_log$events[[2]]$params$filter_min_features, "0")
    })
})

test_that("refused trajectory and SingleR runs are not logged and a completed one is", {
    rds_path <- suppressWarnings(suppressMessages(clustered_fixture_path()))

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        session$setInputs(
            `traj-rds_source` = "upload",
            `traj-trajectory_method` = "slingshot",
            `traj-set_start` = 1,
            `celltype-data_source` = "single"
        )
        # No processed RDS and no clustered single sample yet.
        session$setInputs(`traj-run_ti` = 1)
        session$setInputs(`celltype-run_singler` = 1)
        expect_length(session_log$events, 0)

        session$setInputs(`traj-rds_file` = list(datapath = rds_path))
        # The pseudotime plot warns about plotting a meta.data feature.
        suppressWarnings(session$setInputs(`traj-run_ti` = 2))
        expect_false(is.null(traj$ti_result()))
        expect_equal(logged_actions(session_log$events), "Run trajectory inference")
    })
})

test_that("a refused advanced plot is not logged and a generated one is", {
    # An Advanced Plots module that returns its plot_state: NULL when the click
    # is refused, the plot state when a plot is made.
    testthat::local_mocked_bindings(
        mod_advanced_plots_server = function(id, ...) {
            shiny::moduleServer(id, function(input, output, session) {
                plot_state <- shiny::eventReactive(input$make_plot, {
                    if (is.null(input$gene_list)) {
                        return(NULL)
                    }
                    list(genes = input$gene_list, plot_type = input$plot_type)
                })
                list(plot_state = plot_state)
            })
        },
        .package = "ascseurat"
    )

    shiny::testServer(shiny::shinyAppFile(app_file()), {
        session$setInputs(`adv_plots-plot_type` = "violin")
        session$setInputs(`adv_plots-make_plot` = 1)
        expect_length(session_log$events, 0)

        session$setInputs(`adv_plots-gene_list` = "CD3E")
        session$setInputs(`adv_plots-make_plot` = 2)
        session$setInputs(`adv_plots-plot_type` = "dot")
        expect_equal(logged_actions(session_log$events), "Generate advanced plot")
        expect_equal(session_log$events[[1]]$params$plot_type, "violin")
    })
})

test_that("the session report lists completed demo steps only", {
    skip_on_cran()
    skip_if_not_installed("shinytest2")
    skip_if_not_installed("chromote")

    app_path <- app_file()
    wrapper_dir <- withr::local_tempdir()
    # Run the source version of the package when it is loaded with pkgload.
    setup <- if (pkgload::is_dev_package("ascseurat")) {
        sprintf(
            "pkgload::load_all(%s, quiet = TRUE, helpers = FALSE)",
            deparse(dirname(dirname(dirname(app_path))))
        )
    }
    writeLines(
        c(setup, sprintf("shiny::shinyAppFile(%s)", deparse(app_path))),
        file.path(wrapper_dir, "app.R")
    )

    # Chrome can take longer than the default 10 seconds to start on a busy machine.
    withr::local_options(chromote.timeout = 60)
    app <- shinytest2::AppDriver$new(
        wrapper_dir,
        name = "session-report",
        load_timeout = 120000,
        timeout = 120000,
        shiny_args = list(port = httpuv::randomPort(min = 5000L, max = 5900L))
    )
    withr::defer(app$stop())
    wait <- function() app$wait_for_idle(duration = 1000, timeout = 300000)

    app$set_inputs(main_nav = "demo", wait_ = FALSE)
    app$wait_for_value(output = "demo_data-data_ready", ignore = list(NULL, "false"))
    app$click("demo_norm-run_norm", wait_ = FALSE)
    app$wait_for_value(output = "demo_norm-norm_ready", ignore = list(NULL, "false"),
                       timeout = 300000)
    app$set_inputs(`demo_clust-run_tsne` = FALSE, wait_ = FALSE)
    app$click("demo_clust-run_clustering", wait_ = FALSE)
    app$wait_for_value(output = "demo_clust-clustering_ready", ignore = list(NULL, "false"),
                       timeout = 300000)
    wait()

    # Find Markers with a blank p-value cutoff is refused.
    app$set_inputs(`demo_de-de_type` = "one", wait_ = FALSE)
    wait()
    app$set_inputs(`demo_de-pval_cutoff` = NA, wait_ = FALSE)
    app$click("demo_de-run_de", wait_ = FALSE)
    wait()
    app$set_inputs(`demo_de-pval_cutoff` = 0.05, wait_ = FALSE)
    app$click("demo_de-run_de", wait_ = FALSE)
    wait()
    # Changed after the run, so not part of the logged step.
    app$set_inputs(`demo_de-pval_cutoff` = 0.2, wait_ = FALSE)
    wait()

    app$set_inputs(`demo_clust-rename_1` = "Renamed", wait_ = FALSE)
    app$click("demo_clust-apply_rename", wait_ = FALSE)
    wait()

    # Reanalyze Subset with no cluster selected is refused.
    app$set_inputs(`demo_clust-filter_clusters` = "1", wait_ = FALSE)
    app$set_inputs(`demo_clust-reanalyze_strategy` = "keep", wait_ = FALSE)
    wait()
    app$click("demo_clust-reanalyze", wait_ = FALSE)
    wait()
    app$set_inputs(`demo_clust-selected_clusters` = "Renamed", wait_ = FALSE)
    app$click("demo_clust-reanalyze", wait_ = FALSE)
    wait()
    # Recomputing the PCA completes by clearing the clustered object.
    app$set_inputs(`demo_clust-reanalyze_strategy` = "recompute", wait_ = FALSE)
    app$set_inputs(`demo_clust-selected_clusters` = "Renamed", wait_ = FALSE)
    app$click("demo_clust-reanalyze", wait_ = FALSE)
    app$wait_for_value(output = "demo_clust-reanalyze_pending", ignore = list(NULL, "false"),
                       timeout = 300000)
    wait()

    app$click("session_report", wait_ = FALSE)
    wait()
    report <- readLines(app$get_download("download_session_report"))
    steps <- report_steps(report)

    expect_equal(
        vapply(steps, function(step) step$action, character(1)),
        c(
            "Load demo data (demo tab)",
            "Normalize (demo)",
            "Cluster (demo)",
            "Run DE (demo)",
            "Rename clusters (demo)",
            "Reanalyze subset (demo)",
            "Reanalyze subset (demo)"
        )
    )
    expect_true("pval_cutoff: 0.05" %in% steps[[4]]$params)
    expect_true("rename_1: Renamed" %in% steps[[5]]$params)
    expect_true("selected_clusters: Renamed" %in% steps[[6]]$params)
    expect_true("reanalyze_strategy: keep" %in% steps[[6]]$params)
    expect_true("reanalyze_strategy: recompute" %in% steps[[7]]$params)
})
