loading_rds_path <- local({
    cached_path <- NULL
    function() {
        if (is.null(cached_path) || !file.exists(cached_path)) {
            obj <- demo_fixture()
            # Drop percent.mt so loading has to compute it from the regex
            obj$percent.mt <- NULL
            cached_path <<- tempfile("ascseurat_loading_", fileext = ".rds")
            saveRDS(obj, cached_path)
        }
        cached_path
    }
})

upload_input <- function(path) {
    list(
        name = basename(path),
        datapath = path,
        size = file.size(path),
        type = ""
    )
}

local_notifications <- function(env = parent.frame()) {
    notes <- new.env()
    notes$messages <- character()
    testthat::local_mocked_bindings(
        showNotification = function(ui, ...) {
            notes$messages <- c(notes$messages, as.character(ui))
            invisible("notification")
        },
        .env = env
    )
    notes
}

test_that("an invalid mitochondrial regex is reported without closing the session", {
    notes <- local_notifications()

    shiny::testServer(ascseurat::mod_data_loading_server, {
        session$setInputs(
            input_type = "rds",
            data_file = upload_input(loading_rds_path()),
            mito_pattern = "^MT-(",
            load_data = 1
        )

        expect_false(session$isClosed())
        expect_null(state_obj())
        expect_true(any(grepl("Invalid mitochondrial gene pattern", notes$messages)))

        session$setInputs(mito_pattern = "^MT-", load_data = 2)

        expect_false(session$isClosed())
        expect_s4_class(state_obj(), "Seurat")
        expect_true(any(state_obj()$percent.mt > 0))
    })
})

test_that("metadata with duplicated, NA or blank barcodes is rejected", {
    notes <- local_notifications()
    cells <- colnames(demo_fixture())[1:20]
    good <- data.frame(barcode = cells, group = rep(c("a", "b"), 10))
    tables <- list(
        duplicated = rbind(good, data.frame(barcode = cells[[1]], group = "c")),
        na = rbind(good, data.frame(barcode = NA, group = "c")),
        blank = rbind(good, data.frame(barcode = c("", ""), group = c("c", "d")))
    )

    shiny::testServer(ascseurat::mod_data_loading_server, {
        session$setInputs(
            input_type = "rds",
            data_file = upload_input(loading_rds_path()),
            mito_pattern = "",
            load_data = 1
        )
        expect_s4_class(state_obj(), "Seurat")

        clicks <- 0
        for (case in names(tables)) {
            csv_path <- tempfile(paste0("metadata_", case, "_"), fileext = ".csv")
            utils::write.csv(tables[[case]], csv_path, row.names = FALSE)
            clicks <- clicks + 1
            notes$messages <- character()
            session$setInputs(metadata_file = upload_input(csv_path),
                              apply_metadata = clicks)

            expect_false(session$isClosed())
            expect_false("group" %in% colnames(state_obj()@meta.data))
            expect_true(any(grepl("unique cell barcode", notes$messages)), info = case)
        }

        csv_path <- tempfile("metadata_good_", fileext = ".csv")
        utils::write.csv(good, csv_path, row.names = FALSE)
        session$setInputs(metadata_file = upload_input(csv_path),
                          apply_metadata = clicks + 1)

        expect_false(session$isClosed())
        expect_equal(state_obj()@meta.data[cells, "group"], good$group)
    })
})

test_that("unexpected metadata errors are reported without closing the session", {
    notes <- local_notifications()
    empty_path <- tempfile("metadata_empty_", fileext = ".csv")
    file.create(empty_path)

    shiny::testServer(ascseurat::mod_data_loading_server, {
        session$setInputs(
            input_type = "rds",
            data_file = upload_input(loading_rds_path()),
            mito_pattern = "",
            load_data = 1
        )
        session$setInputs(metadata_file = upload_input(empty_path),
                          apply_metadata = 1)

        expect_false(session$isClosed())
        expect_s4_class(state_obj(), "Seurat")
        expect_true(any(grepl("Adding metadata failed", notes$messages)))
    })
})
