test_that("SingleR module runs when cached references are available", {
    skip_if_not_installed("SingleR")
    skip_if_not_installed("celldex")

    cache_file <- file.path(tools::R_user_dir("ascseurat", "cache"), "celldex", "hpca.rds")
    skip_if_not(file.exists(cache_file), "SingleR reference cache is not available")

    shiny::testServer(
        ascseurat::mod_celltype_annotation_server,
        args = list(
            seurat_obj_single = shiny::reactive(clustered_fixture()),
            seurat_obj_multi = shiny::reactive(clustered_fixture())
        ),
        {
            session$setInputs(
                data_source = "single",
                reference = "hpca",
                granularity = "label.main",
                run_singler = 1
            )
            expect_true("singler_label" %in% colnames(annotated_obj()@meta.data))
        }
    )
})

test_that("SingleR expression prep joins split assay layers", {
    skip_if_not("JoinLayers" %in% getNamespaceExports("SeuratObject"))

    counts <- matrix(rpois(10 * 8, 5), nrow = 10)
    rownames(counts) <- paste0("Gene", seq_len(nrow(counts)))
    colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))

    obj <- Seurat::CreateSeuratObject(counts)
    obj$sample <- rep(c("A", "B"), each = 4)
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)

    expect_error(
        Seurat::GetAssayData(obj, assay = "RNA", layer = "data"),
        "multiple layers"
    )

    prepared <- ascseurat:::prepare_singler_input(obj)
    expect_equal(dim(prepared$expr_matrix), c(nrow(obj), ncol(obj)))
    expect_true("data" %in% SeuratObject::Layers(prepared$obj[["RNA"]]))
    expect_false(any(grepl("^data[.]", SeuratObject::Layers(prepared$obj[["RNA"]]))))
})

# A small reference built from the fixture, cached as "hpca" in a temporary
# cache directory, so SingleR runs in seconds without celldex downloads.
local_singler_reference <- function(obj, env = parent.frame()) {
    cache_dir <- withr::local_tempdir(.local_envir = env)
    withr::local_envvar(R_USER_CACHE_DIR = cache_dir, .local_envir = env)
    ref_file <- file.path(tools::R_user_dir("ascseurat", "cache"), "celldex", "hpca.rds")
    dir.create(dirname(ref_file), recursive = TRUE)
    cells <- colnames(obj)[seq(1, ncol(obj), by = 10)]
    ref <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            logcounts = as.matrix(SeuratObject::LayerData(obj, layer = "data")[, cells])
        ),
        colData = S4Vectors::DataFrame(
            label.main = paste0("type_", as.character(SeuratObject::Idents(obj)[cells]))
        )
    )
    saveRDS(ref, ref_file)
    cache_dir
}

renamed_fixture <- function() {
    obj <- suppressWarnings(clustered_fixture())
    SeuratObject::Idents(obj) <- paste0("c", as.character(SeuratObject::Idents(obj)))
    obj
}

test_that("SingleR results are cleared when the source tab replaces or clears its object", {
    skip_if_not_installed("SingleR")
    skip_if_not_installed("SummarizedExperiment")
    local_mocked_bindings(
        showNotification = function(ui, ...) invisible("id"),
        .package = "ascseurat"
    )
    obj <- suppressWarnings(clustered_fixture())
    renamed <- renamed_fixture()
    local_singler_reference(obj)
    current_obj <- shiny::reactiveVal(obj)
    # Re-runs the upstream reactive without changing its object.
    refresh <- shiny::reactiveVal(0)

    shiny::testServer(
        ascseurat::mod_celltype_annotation_server,
        args = list(seurat_obj_single = shiny::reactive({
            refresh()
            shiny::req(current_obj())
            current_obj()
        })),
        {
            session$setInputs(
                data_source = "single",
                reference = "hpca",
                granularity = "label.main",
                run_singler = 1
            )
            expect_true("singler_label" %in% colnames(annotated_obj()@meta.data))

            # The same object again keeps the results.
            refresh(refresh() + 1)
            session$flushReact()
            expect_false(is.null(annotated_obj()))
            expect_identical(output$singler_state, "ready")

            # Renamed or re-clustered object in the Single Sample tab.
            current_obj(renamed)
            session$flushReact()
            expect_false(session$isClosed())
            expect_null(annotated_obj())
            expect_null(singler_result())
            expect_identical(output$singler_state, "outdated")

            # Annotating again uses the current object.
            session$setInputs(run_singler = 2)
            expect_identical(output$singler_state, "ready")
            expect_identical(
                levels(SeuratObject::Idents(annotated_obj())),
                levels(SeuratObject::Idents(renamed))
            )

            # Choosing another source does not clear results that are current.
            session$setInputs(data_source = "upload")
            expect_identical(output$singler_state, "ready")

            # The Single Sample object is no longer available.
            current_obj(NULL)
            session$flushReact()
            expect_false(session$isClosed())
            expect_null(annotated_obj())
            expect_identical(output$singler_state, "outdated")
        }
    )
})

test_that("SingleR results of an uploaded object ignore the other tabs", {
    skip_if_not_installed("SingleR")
    skip_if_not_installed("SummarizedExperiment")
    local_mocked_bindings(
        showNotification = function(ui, ...) invisible("id"),
        .package = "ascseurat"
    )
    obj <- suppressWarnings(clustered_fixture())
    local_singler_reference(obj)
    current_obj <- shiny::reactiveVal(obj)

    shiny::testServer(
        ascseurat::mod_celltype_annotation_server,
        args = list(seurat_obj_single = shiny::reactive({
            shiny::req(current_obj())
            current_obj()
        })),
        {
            session$setInputs(
                data_source = "upload",
                rds_file = list(datapath = clustered_fixture_path()),
                reference = "hpca",
                granularity = "label.main",
                run_singler = 1
            )
            expect_identical(output$singler_state, "ready")
            current_obj(NULL)
            session$flushReact()
            expect_identical(output$singler_state, "ready")
            expect_false(is.null(annotated_obj()))
        }
    )
})

test_that("SingleR downloads are hidden once the results are cleared", {
    skip_on_cran()
    skip_if_not_installed("SingleR")
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("shinytest2")
    skip_if_not_installed("chromote")
    skip_if(is.null(suppressMessages(chromote::find_chrome())), "Chrome is not available")
    root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
    skip_if_not(file.exists(file.path(root, "DESCRIPTION")), "Package sources are not available")

    cache_dir <- local_singler_reference(suppressWarnings(clustered_fixture()))
    app_dir <- tempfile("singler_app_")
    dir.create(app_dir)
    writeLines(c(
        sprintf("Sys.setenv(R_USER_CACHE_DIR = %s)", deparse(cache_dir)),
        sprintf("pkgload::load_all(%s, quiet = TRUE)", deparse(root)),
        sprintf("obj <- readRDS(%s)", deparse(clustered_fixture_path())),
        "renamed <- obj",
        "SeuratObject::Idents(renamed) <- paste0('c', as.character(SeuratObject::Idents(obj)))",
        "ui <- bslib::page_fluid(shiny::actionButton('swap', 'Swap'), mod_celltype_annotation_ui('celltype'))",
        "server <- function(input, output, session) {",
        "    current_obj <- shiny::reactiveVal(obj)",
        "    shiny::observeEvent(input$swap, current_obj(renamed))",
        "    mod_celltype_annotation_server('celltype', seurat_obj_single = current_obj)",
        "}",
        "shiny::shinyApp(ui, server)"
    ), file.path(app_dir, "app.R"))

    app <- shinytest2::AppDriver$new(
        app_dir,
        load_timeout = 120000,
        timeout = 60000,
        shiny_args = list(port = httpuv::randomPort(min = 5000, max = 5900))
    )
    on.exit(app$stop(), add = TRUE)
    visible <- "$('#celltype-download_annotated').is(':visible')"
    note_visible <- "$('p:contains(\"so the results were cleared\")').is(':visible')"

    app$set_inputs(`celltype-data_source` = "single")
    app$click("celltype-run_singler")
    app$wait_for_js(visible, timeout = 60000)
    expect_true(app$get_js(visible))
    expect_false(app$get_js(note_visible))

    app$click("swap")
    app$wait_for_js(sprintf("!%s", visible), timeout = 60000)
    expect_false(app$get_js(visible))
    expect_true(app$get_js(note_visible))
})
