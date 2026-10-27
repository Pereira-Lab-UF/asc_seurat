test_that("DE module finds markers for a clustered object", {
    obj <- clustered_fixture()

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(seurat_obj = shiny::reactive(obj), allow_conserved = FALSE),
        {
            session$setInputs(
                de_type = "all",
                test_use = "wilcox",
                pval_cutoff = 0.99,
                logfc_threshold = 0,
                min_pct = 0,
                only_pos = FALSE,
                run_de = 1
            )
            res <- markers()
            expect_s3_class(res, "data.frame")
            expect_gt(nrow(res), 0)
            expect_identical(colnames(res)[1:2], c("gene", "cluster"))
            expect_identical(rownames(res), as.character(seq_len(nrow(res))))
        }
    )
})

test_that("DE module compares metadata conditions inside a selected cluster", {
    obj <- clustered_fixture()
    cluster_counts <- table(as.character(Seurat::Idents(obj)))
    cluster <- names(cluster_counts)[cluster_counts >= 6][[1]]
    cluster_cells <- names(Seurat::Idents(obj))[as.character(Seurat::Idents(obj)) == cluster]

    obj$test_condition <- "other"
    obj$test_condition[cluster_cells] <- rep(c("A", "B"), length.out = length(cluster_cells))

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(seurat_obj = shiny::reactive(obj), allow_conserved = FALSE),
        {
            session$setInputs(
                de_type = "conditions",
                cluster_1 = cluster,
                condition_column = "test_condition",
                condition_value_1 = "A",
                condition_value_2 = "B",
                test_use = "wilcox",
                pval_cutoff = 0.999999,
                logfc_threshold = 0,
                min_pct = 0,
                only_pos = FALSE,
                run_de = 1
            )
            res <- markers()
            expect_s3_class(res, "data.frame")
            expect_identical(colnames(res)[1:2], c("gene", "cluster"))
            expect_false(anyNA(res$cluster))
            expect_false(anyNA(res$condition_column))
            expect_true(all(res$cluster == cluster))
            expect_true(all(res$condition_column == "test_condition"))
            expect_identical(rownames(res), as.character(seq_len(nrow(res))))
        }
    )
})

# Records showNotification() messages raised by the DE module.
local_de_notifications <- function(env = parent.frame()) {
    messages <- character()
    testthat::local_mocked_bindings(
        showNotification = function(ui, ...) {
            messages <<- c(messages, as.character(ui))
            invisible(NULL)
        },
        .package = "ascseurat",
        .env = env
    )
    function() messages
}

largest_cluster <- function(obj) {
    sizes <- table(as.character(Seurat::Idents(obj)))
    names(sizes)[which.max(sizes)]
}

test_that("blank cutoffs keep the DE session open and explain why nothing ran", {
    obj <- clustered_fixture()
    notifications <- local_de_notifications()

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(seurat_obj = shiny::reactive(obj), allow_conserved = FALSE),
        {
            session$setInputs(
                de_type = "all",
                test_use = "wilcox",
                pval_cutoff = NA,
                logfc_threshold = 0.25,
                min_pct = 0.1,
                only_pos = TRUE
            )
            expect_false(session$isClosed())
            session$setInputs(run_de = 1)
            expect_false(session$isClosed())
            expect_error(markers(), class = "shiny.silent.error")
            expect_match(notifications(), "adjusted p-value cutoff", all = FALSE)

            session$setInputs(pval_cutoff = 0.05, min_pct = NA)
            expect_false(session$isClosed())
            session$setInputs(run_de = 2)
            expect_false(session$isClosed())
            expect_error(markers(), class = "shiny.silent.error")
            expect_match(notifications(), "minimum fraction of cells", all = FALSE)

            session$setInputs(min_pct = 0.1, logfc_threshold = NA)
            session$setInputs(run_de = 3)
            expect_false(session$isClosed())
            expect_error(markers(), class = "shiny.silent.error")
            expect_match(notifications(), "log2 fold-change threshold", all = FALSE)
        }
    )
})

test_that("conserved markers use the chosen test and the largest adjusted p-value", {
    obj <- clustered_fixture()
    obj$sample_group <- rep(c("s1", "s2"), length.out = ncol(obj))
    cluster <- largest_cluster(obj)
    notifications <- local_de_notifications()

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(seurat_obj = shiny::reactive(obj), allow_conserved = TRUE),
        {
            session$setInputs(
                de_type = "conserved",
                cluster_1 = cluster,
                conserved_grouping = "sample_group",
                test_use = "wilcox",
                pval_cutoff = 0.05,
                logfc_threshold = 0.25,
                min_pct = 0.1,
                only_pos = TRUE,
                run_de = 1
            )
            res_wilcox <- markers()
            expect_gt(nrow(res_wilcox), 0)
            expect_equal(
                res_wilcox$max_p_val_adj,
                pmax(res_wilcox$s1_p_val_adj, res_wilcox$s2_p_val_adj)
            )
            expect_true(all(res_wilcox$max_p_val_adj < 0.05))

            session$setInputs(test_use = "t", run_de = 2)
            res_t <- markers()
            expect_false(identical(res_wilcox$s1_p_val, res_t$s1_p_val))

            session$setInputs(test_use = "roc", run_de = 3)
            expect_error(markers(), class = "shiny.silent.error")
            expect_match(notifications(), "'roc' test", all = FALSE)
        }
    )
})

test_that("conserved markers refuse clusters testable in fewer than two groups", {
    obj <- clustered_fixture()
    idents <- as.character(Seurat::Idents(obj))
    cluster <- largest_cluster(obj)
    cluster_cells <- colnames(obj)[idents == cluster]
    small_cluster <- names(sort(table(idents)))[[1]]
    small_cells <- colnames(obj)[idents == small_cluster]

    obj$one_group <- rep(c("s1", "s2"), length.out = ncol(obj))
    obj$one_group[cluster_cells] <- "s1"
    obj$one_group[cluster_cells[1:2]] <- "s2"
    obj$no_group <- rep(c("g1", "g2"), length.out = ncol(obj))
    obj$no_group[small_cells] <- paste0("g", ceiling(seq_along(small_cells) / 2))
    notifications <- local_de_notifications()

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(seurat_obj = shiny::reactive(obj), allow_conserved = TRUE),
        {
            session$setInputs(
                de_type = "conserved",
                cluster_1 = cluster,
                conserved_grouping = "one_group",
                test_use = "wilcox",
                pval_cutoff = 0.05,
                logfc_threshold = 0.25,
                min_pct = 0.1,
                only_pos = TRUE,
                run_de = 1
            )
            expect_error(markers(), class = "shiny.silent.error")
            expect_match(notifications(), "only one group ('s1')", fixed = TRUE, all = FALSE)
            expect_match(notifications(), "Find markers for one cluster", all = FALSE)

            session$setInputs(
                cluster_1 = small_cluster,
                conserved_grouping = "no_group",
                run_de = 2
            )
            expect_error(markers(), class = "shiny.silent.error")
            expect_match(notifications(), "cells in no group", all = FALSE)
        }
    )
})

test_that("DE between conditions offers user-added metadata columns", {
    obj <- clustered_fixture()

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(seurat_obj = shiny::reactive(obj), allow_conserved = FALSE),
        {
            session$setInputs(de_type = "conditions")
            expect_length(condition_columns(), 0)
            expect_match(output$condition_column_ui$html,
                         "No metadata column", fixed = TRUE)
        }
    )

    obj$treatment <- rep(c("ctrl", "drug"), length.out = ncol(obj))
    obj$time_point <- rep(c(0L, 24L), length.out = ncol(obj))
    obj$samples <- rep(c("WT", "rhd6"), length.out = ncol(obj))
    cluster <- largest_cluster(obj)

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(seurat_obj = shiny::reactive(obj), allow_conserved = FALSE),
        {
            session$setInputs(de_type = "conditions", cluster_1 = cluster)
            cols <- condition_columns()
            expect_true(all(c("samples", "treatment", "time_point") %in% cols))
            expect_false(any(c(
                "orig.ident", "seurat_clusters", "nCount_RNA",
                "nFeature_RNA", "percent.mt"
            ) %in% cols))
            expect_false(any(grepl("_snn_res\\.", cols)))
            expect_match(output$condition_column_ui$html,
                         '<option value="samples" selected>', fixed = TRUE)

            session$setInputs(condition_column = "treatment")
            expect_identical(condition_values(), c("ctrl", "drug"))
        }
    )
})

test_that("DE results are cleared when the clustered object changes", {
    obj <- clustered_fixture()
    cluster <- largest_cluster(obj)
    current_obj <- shiny::reactiveVal(obj)

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(
            seurat_obj = shiny::reactive({
                shiny::req(current_obj())
                current_obj()
            }),
            allow_conserved = FALSE
        ),
        {
            session$setInputs(
                de_type = "one",
                cluster_1 = cluster,
                test_use = "wilcox",
                pval_cutoff = 0.05,
                logfc_threshold = 0.25,
                min_pct = 0.1,
                only_pos = TRUE,
                run_de = 1
            )
            expect_true(all(markers()$cluster == cluster))

            renamed <- obj
            Seurat::Idents(renamed) <- factor(paste0("C", as.character(Seurat::Idents(obj))))
            current_obj(renamed)
            session$flushReact()
            expect_error(markers(), class = "shiny.silent.error")

            session$setInputs(cluster_1 = paste0("C", cluster), run_de = 2)
            expect_true(all(markers()$cluster == paste0("C", cluster)))

            current_obj(NULL)
            session$flushReact()
            expect_error(markers(), class = "shiny.silent.error")
        }
    )
})

test_that("DE results survive upstream re-evaluations that return the same object", {
    obj <- clustered_fixture()
    cluster <- largest_cluster(obj)
    tick <- shiny::reactiveVal(0)

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(
            seurat_obj = shiny::reactive({
                tick()
                obj
            }),
            allow_conserved = FALSE
        ),
        {
            session$setInputs(
                de_type = "one",
                cluster_1 = cluster,
                test_use = "wilcox",
                pval_cutoff = 0.05,
                logfc_threshold = 0.25,
                min_pct = 0.1,
                only_pos = TRUE,
                run_de = 1
            )
            res <- markers()
            expect_gt(nrow(res), 0)

            tick(1)
            session$flushReact()
            expect_identical(markers(), res)
        }
    )
})

test_that("cleared DE results report not ready and the download explains why", {
    obj <- clustered_fixture()
    cluster <- largest_cluster(obj)
    current_obj <- shiny::reactiveVal(obj)
    notifications <- local_de_notifications()

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(seurat_obj = shiny::reactive(current_obj()), allow_conserved = FALSE),
        {
            expect_identical(output$de_ready, "false")
            session$setInputs(
                de_type = "one",
                cluster_1 = cluster,
                test_use = "wilcox",
                pval_cutoff = 0.05,
                logfc_threshold = 0.25,
                min_pct = 0.1,
                only_pos = TRUE,
                run_de = 1
            )
            expect_identical(output$de_ready, "true")
            expect_gt(nrow(utils::read.csv(output$download_markers)), 0)

            renamed <- obj
            Seurat::Idents(renamed) <- factor(paste0("C", as.character(Seurat::Idents(obj))))
            current_obj(renamed)
            session$flushReact()
            expect_identical(output$de_ready, "false")
            # The mock session errors when no file is written; a real session
            # answers 404 instead of a server error (see the browser test).
            expect_error(output$download_markers)
            expect_match(notifications(), "no marker results to download", all = FALSE)
            expect_false(session$isClosed())
        }
    )
})

test_that("the browser hides DE results and their download once they are cleared", {
    skip_on_cran()
    skip_if_not_installed("shinytest2")
    skip_if(
        is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL)),
        "Chrome is not available"
    )

    app_dir <- withr::local_tempdir()
    writeLines(c(
        "library(shiny)",
        "library(ascseurat)",
        sprintf("obj <- readRDS(%s)", deparse(clustered_fixture_path())),
        "ui <- bslib::page_fluid(",
        "    shinyFeedback::useShinyFeedback(),",
        "    actionButton('rename', 'Rename clusters'),",
        "    mod_de_ui('de')",
        ")",
        "server <- function(input, output, session) {",
        "    current <- reactiveVal(obj)",
        "    observeEvent(input$rename, {",
        "        renamed <- obj",
        "        Seurat::Idents(renamed) <- factor(paste0('C', as.character(Seurat::Idents(obj))))",
        "        current(renamed)",
        "    })",
        "    mod_de_server('de', seurat_obj = reactive(current()))",
        "}",
        "shinyApp(ui, server)"
    ), file.path(app_dir, "app.R"))

    app <- shinytest2::AppDriver$new(
        app_dir,
        name = "mod-de-ready",
        shiny_args = list(port = httpuv::randomPort(5000, 5900)),
        load_timeout = 180000,
        timeout = 120000
    )
    withr::defer(app$stop())

    download_visible <- "$('#de-download_markers').is(':visible')"
    card_visible <- paste0(
        "$('.card-header').filter(function() {",
        " return $(this).text().trim() === 'Marker Genes'; }).is(':visible')"
    )
    download_status <- paste0(
        "fetch($('#de-download_markers').attr('href'))",
        ".then(function(response) { return response.status; })"
    )

    expect_false(app$get_js(download_visible))
    app$click("de-run_de")
    app$wait_for_idle(duration = 1000, timeout = 120000)
    expect_true(app$get_js(download_visible))
    expect_true(app$get_js(card_visible))
    expect_equal(app$get_js(download_status), 200)

    app$click("rename")
    app$wait_for_idle(duration = 1000)
    expect_false(app$get_js(download_visible))
    expect_false(app$get_js(card_visible))
    expect_equal(app$get_js(download_status), 404)
    app$wait_for_js(
        "$('.shiny-notification').text().indexOf('no marker results to download') >= 0"
    )
    expect_true(app$get_js("Shiny.shinyapp.isConnected()"))
})
