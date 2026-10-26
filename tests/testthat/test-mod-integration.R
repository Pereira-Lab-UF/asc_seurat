test_that("QC preview cache key changes when a 10X file is overwritten in place", {
    sample_dir <- withr::local_tempdir()
    matrix_file <- file.path(sample_dir, "matrix.mtx")
    writeLines("old", matrix_file)

    shiny::testServer(ascseurat::mod_integration_server, {
        before <- input_version_key(sample_dir)
        writeLines("new contents", matrix_file)
        expect_false(identical(input_version_key(sample_dir), before))
    })
})

# Writes demo cells as an uncompressed 10X (v2) directory.
write_10x_fixture <- function(dir, cells) {
    counts <- SeuratObject::LayerData(demo_fixture()[, cells], layer = "counts")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    entries <- cbind(
        counts@i + 1L,
        rep(seq_len(ncol(counts)), diff(counts@p)),
        counts@x
    )
    writeLines(
        c("%%MatrixMarket matrix coordinate integer general",
          paste(nrow(counts), ncol(counts), nrow(entries))),
        file.path(dir, "matrix.mtx")
    )
    utils::write.table(entries, file.path(dir, "matrix.mtx"), append = TRUE,
                       row.names = FALSE, col.names = FALSE)
    writeLines(colnames(counts), file.path(dir, "barcodes.tsv"))
    writeLines(paste(rownames(counts), rownames(counts), sep = "\t"),
               file.path(dir, "genes.tsv"))
    dir
}

rendered_input_value <- function(html, id) {
    tag <- regmatches(html, regexpr(paste0('<input id="', id, '"[^>]*>'), html))
    sub('.*value="([^"]*)".*', "\\1", tag)
}

local_numeric_updates <- function(env = parent.frame()) {
    updates <- new.env()
    updates$calls <- list()
    testthat::local_mocked_bindings(
        updateNumericInput = function(session, inputId, ...) {
            updates$calls[[length(updates$calls) + 1]] <- c(list(id = inputId), list(...))
        },
        .package = "ascseurat",
        .env = env
    )
    updates
}

test_that("QC controls keep typed values when a sample name or path changes", {
    shiny::testServer(ascseurat::mod_integration_server, {
        session$setInputs(
            n_samples = 2,
            sample_name_1 = "WT", sample_path_1 = "",
            sample_name_2 = "rhd6", sample_path_2 = "",
            qc_mito_pattern_1 = "^ATM", qc_min_genes_1 = 300,
            qc_max_genes_1 = 2500, qc_max_mito_1 = 10
        )
        session$setInputs(sample_name_1 = "WT_renamed", sample_path_2 = "data/rhd6")
        expect_equal(output$qc_title_1, "WT_renamed")
        session$setInputs(n_samples = 3)

        html <- output$qc_params_ui$html
        expect_equal(rendered_input_value(html, session$ns("qc_mito_pattern_1")), "^ATM")
        expect_equal(rendered_input_value(html, session$ns("qc_min_genes_1")), "300")
        expect_equal(rendered_input_value(html, session$ns("qc_max_genes_1")), "2500")
        expect_equal(rendered_input_value(html, session$ns("qc_max_mito_1")), "10")
        expect_equal(rendered_input_value(html, session$ns("qc_min_genes_2")), "200")
        expect_equal(rendered_input_value(html, session$ns("qc_min_genes_3")), "200")
    })
})

test_that("editing a sample name or path does not re-render the QC controls", {
    shiny::testServer(ascseurat::mod_integration_server, {
        session$setInputs(
            n_samples = 2,
            sample_name_1 = "WT", sample_path_1 = "",
            sample_name_2 = "rhd6", sample_path_2 = ""
        )
        before <- output$qc_params_ui$html

        # A re-render would replace what the browser holds with the server's
        # copy, which can lag behind (a suggestion or typing not yet sent back).
        session$setInputs(qc_mito_pattern_1 = "stale", qc_min_genes_1 = 1)
        session$setInputs(sample_name_2 = "rhd6_renamed", sample_path_2 = "data/rhd6")

        expect_identical(output$qc_params_ui$html, before)
        expect_equal(output$qc_title_2, "rhd6_renamed")
    })
})

test_that("changing the number of samples keeps typed sample names and paths", {
    shiny::testServer(ascseurat::mod_integration_server, {
        session$setInputs(
            n_samples = 2,
            sample_name_1 = "WT", sample_path_1 = "data/wt",
            sample_name_2 = "rhd6", sample_path_2 = "data/rhd6"
        )
        session$setInputs(n_samples = 3)

        html <- output$sample_table_ui$html
        expect_equal(rendered_input_value(html, session$ns("sample_name_1")), "WT")
        expect_equal(rendered_input_value(html, session$ns("sample_path_1")), "data/wt")
        expect_equal(rendered_input_value(html, session$ns("sample_name_2")), "rhd6")
        expect_equal(rendered_input_value(html, session$ns("sample_path_2")), "data/rhd6")
        expect_equal(rendered_input_value(html, session$ns("sample_name_3")), "sample_3")
        expect_equal(rendered_input_value(html, session$ns("sample_path_3")), "")
    })
})

test_that("QC suggestions only update samples whose data path changed", {
    root <- withr::local_tempdir()
    path_1 <- write_10x_fixture(file.path(root, "s1"), 1:300)
    path_2 <- write_10x_fixture(file.path(root, "s2"), 301:600)
    updates <- local_numeric_updates()
    updated_ids <- function() vapply(updates$calls, `[[`, character(1), "id")

    shiny::testServer(ascseurat::mod_integration_server, {
        session$setInputs(
            n_samples = 2,
            sample_name_1 = "s1", sample_path_1 = "",
            sample_name_2 = "s2", sample_path_2 = ""
        )
        session$elapse(1500)

        suppressMessages({
            session$setInputs(sample_path_1 = path_1)
            session$elapse(1500)
        })
        expect_setequal(
            updated_ids(),
            c("qc_min_genes_1", "qc_max_genes_1", "qc_max_mito_1")
        )

        updates$calls <- list()
        session$setInputs(qc_min_genes_1 = 333, sample_name_2 = "s2_renamed", n_samples = 3)
        session$elapse(1500)
        expect_length(updates$calls, 0)

        suppressMessages({
            session$setInputs(sample_path_2 = path_2)
            session$elapse(1500)
        })
        expect_setequal(
            updated_ids(),
            c("qc_min_genes_2", "qc_max_genes_2", "qc_max_mito_2")
        )
    })
})

test_that("a restored bookmark keeps its QC thresholds", {
    root <- withr::local_tempdir()
    path_1 <- write_10x_fixture(file.path(root, "s1"), 1:300)
    path_2 <- write_10x_fixture(file.path(root, "s2"), 301:600)
    updates <- local_numeric_updates()
    restore <- new.env()
    # MockShinySession ignores onRestore, so capture the callback instead.
    testthat::local_mocked_bindings(
        onRestore = function(fun, session = NULL) restore$fun <- fun,
        .package = "shiny"
    )

    shiny::testServer(ascseurat::mod_integration_server, {
        restored <- list(
            n_samples = 1, sample_name_1 = "s1", sample_path_1 = path_1,
            qc_mito_pattern_1 = "", qc_min_genes_1 = 333,
            qc_max_genes_1 = 2222, qc_max_mito_1 = 7
        )
        restore$fun(list(input = restored))
        # The sample and QC inputs arrive once the rendered controls bind.
        session$setInputs(n_samples = 1)
        session$elapse(1500)
        do.call(session$setInputs, restored)
        suppressMessages(session$elapse(1500))
        expect_length(updates$calls, 0)

        suppressMessages({
            session$setInputs(sample_path_1 = path_2)
            session$elapse(1500)
        })
        expect_length(updates$calls, 3)
    })
})

local_mocked_integration <- function(n_pcs, env = parent.frame()) {
    testthat::local_mocked_bindings(
        run_integration = function(obj, ...) {
            embeddings <- matrix(
                stats::rnorm(ncol(obj) * n_pcs),
                nrow = ncol(obj),
                dimnames = list(colnames(obj), paste0("PC_", seq_len(n_pcs)))
            )
            obj[["pca"]] <- SeuratObject::CreateDimReducObject(
                embeddings = embeddings,
                key = "PC_",
                assay = "RNA"
            )
            obj
        },
        .package = "ascseurat",
        .env = env
    )
}

test_that("the first integration updates the number of PCs for clustering", {
    root <- withr::local_tempdir()
    sample_path <- write_10x_fixture(file.path(root, "s1"), 1:300)
    updates <- local_numeric_updates()
    local_mocked_integration(n_pcs = 12)

    shiny::testServer(ascseurat::mod_integration_server, {
        session$setInputs(
            n_samples = 1,
            sample_name_1 = "s1", sample_path_1 = sample_path,
            qc_mito_pattern_1 = "", qc_min_genes_1 = 0,
            qc_max_genes_1 = 1e6, qc_max_mito_1 = 100,
            int_method = "harmony", norm_method = "lognorm", n_dims = 30
        )
        updates$calls <- list()
        suppressMessages(session$setInputs(run_integration = 1))
        expect_s4_class(integrated_obj(), "Seurat")

        expect_equal(
            Filter(function(call) identical(call$id, "cluster_n_dims"), updates$calls),
            list(list(id = "cluster_n_dims", value = 12, max = 12))
        )
    })
})

test_that("integration reuses the counts read for the QC preview", {
    root <- withr::local_tempdir()
    sample_path <- write_10x_fixture(file.path(root, "s1"), 1:300)
    local_mocked_integration(n_pcs = 12)
    load_calls <- 0
    real_load_input_data <- ascseurat:::load_input_data
    testthat::local_mocked_bindings(
        load_input_data = function(...) {
            load_calls <<- load_calls + 1
            real_load_input_data(...)
        },
        .package = "ascseurat"
    )

    shiny::testServer(ascseurat::mod_integration_server, {
        session$setInputs(n_samples = 1, sample_name_1 = "s1", sample_path_1 = "")
        session$elapse(1500)
        suppressMessages({
            session$setInputs(sample_path_1 = sample_path)
            session$elapse(1500)
        })
        expect_equal(load_calls, 1)
        expect_length(ls(preview_counts, all.names = TRUE), 1)

        suppressMessages(session$setInputs(
            qc_mito_pattern_1 = "", qc_min_genes_1 = 0,
            qc_max_genes_1 = 1e6, qc_max_mito_1 = 100,
            int_method = "harmony", norm_method = "lognorm", n_dims = 30,
            run_integration = 1
        ))
        expect_s4_class(integrated_obj(), "Seurat")
        expect_equal(load_calls, 1)
        expect_length(ls(preview_counts, all.names = TRUE), 0)
    })
})
