test_that("QC module can pass the loaded object through with no filtering", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(session$setInputs(filter_mode = "none"))
            expect_s4_class(qc_result(), "Seurat")
            expect_equal(ncol(qc_result()), ncol(obj))
            expect_identical(
                attr(qc_result(), "ascseurat_filter_mode", exact = TRUE),
                "none"
            )
        }
    )
})

test_that("QC module applies valid filtering thresholds", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(session$setInputs(
                filter_mode = "filter",
                min_genes = 0,
                max_genes = max(obj$nFeature_RNA, na.rm = TRUE) + 1,
                max_mito = 100,
                detect_doublets = FALSE,
                apply_filter = 1
            ))

            filtered <- suppressWarnings(qc_result())
            expect_s4_class(filtered, "Seurat")
            expect_equal(ncol(filtered), ncol(obj))
            expect_identical(
                attr(filtered, "ascseurat_filter_mode", exact = TRUE),
                "filter"
            )
        }
    )
})

test_that("QC module can apply no metric filtering when requested", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(session$setInputs(
                filter_mode = "none",
                detect_doublets = FALSE,
                apply_filter = 1
            ))

            filtered <- suppressWarnings(qc_result())
            expect_s4_class(filtered, "Seurat")
            expect_equal(ncol(filtered), ncol(obj))
            expect_identical(
                attr(filtered, "ascseurat_filter_mode", exact = TRUE),
                "none"
            )
        }
    )
})

test_that("QC module blocks thresholds that remove every cell", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(session$setInputs(
                filter_mode = "filter",
                min_genes = max(obj$nFeature_RNA, na.rm = TRUE) + 1,
                max_genes = max(obj$nFeature_RNA, na.rm = TRUE) + 2,
                max_mito = 100,
                detect_doublets = FALSE,
                apply_filter = 1
            ))

            expect_null(filtered_obj())
            expect_error(suppressWarnings(qc_result()), class = "shiny.silent.error")
        }
    )
})

capture_session_messages <- function(session) {
    captured <- new.env()
    captured$notifications <- character()
    captured$custom <- list()
    root <- session$rootScope()
    root$sendNotification <- function(type, message) {
        if (identical(type, "show")) {
            captured$notifications <- c(captured$notifications, as.character(message$html))
        }
    }
    root$sendCustomMessage <- function(type, message) {
        captured$custom[[length(captured$custom) + 1]] <- list(type = type, message = message)
    }
    captured
}

test_that("QC module flags a cleared max mitochondrial % without closing the session", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            captured <- capture_session_messages(session)
            suppressWarnings(session$setInputs(filter_mode = "filter", max_mito = 5))
            suppressWarnings(session$setInputs(max_mito = NA))

            expect_false(session$isClosed())
            feedback <- Filter(function(msg) identical(msg$type, "showFeedback"), captured$custom)
            expect_length(feedback, 1)
            expect_match(feedback[[1]]$message$text, "between 0 and 100")
        }
    )
})

test_that("QC module refuses to apply filters with a blank threshold", {
    obj <- demo_fixture()

    for (blank_input in c("min_genes", "max_genes", "max_mito")) {
        thresholds <- list(min_genes = 0, max_genes = 1e6, max_mito = 100)
        thresholds[[blank_input]] <- NA

        shiny::testServer(
            ascseurat::mod_qc_server,
            args = list(seurat_obj = shiny::reactive(obj)),
            {
                captured <- capture_session_messages(session)
                suppressWarnings(do.call(session$setInputs, c(
                    list(filter_mode = "filter", detect_doublets = FALSE),
                    thresholds
                )))
                suppressWarnings(session$setInputs(apply_filter = 1))

                expect_false(session$isClosed())
                expect_null(filtered_obj())
                expect_length(captured$notifications, 1)
                expect_match(captured$notifications, "Enter a value for every filtering threshold")
            }
        )
    }
})

test_that("QC result follows filter mode changes made without Apply", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(session$setInputs(
                filter_mode = "filter",
                min_genes = 0,
                max_genes = max(obj$nFeature_RNA, na.rm = TRUE) + 1,
                max_mito = 100,
                detect_doublets = FALSE,
                apply_filter = 1
            ))
            expect_identical(attr(qc_result(), "ascseurat_filter_mode", exact = TRUE), "filter")

            suppressWarnings(session$setInputs(filter_mode = "none"))
            expect_identical(attr(qc_result(), "ascseurat_filter_mode", exact = TRUE), "none")

            suppressWarnings(session$setInputs(filter_mode = "filter"))
            expect_null(filtered_obj())
            expect_error(qc_result(), class = "shiny.silent.error")

            suppressWarnings(session$setInputs(filter_mode = "none", detect_doublets = TRUE))
            expect_error(qc_result(), class = "shiny.silent.error")
        }
    )
})

test_that("QC module drops its filtered object when the loaded data becomes unavailable", {
    obj <- demo_fixture()
    loaded <- shiny::reactiveVal(obj)

    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive({
            shiny::req(loaded())
            loaded()
        })),
        {
            suppressWarnings(session$setInputs(
                filter_mode = "filter",
                min_genes = 0,
                max_genes = max(obj$nFeature_RNA, na.rm = TRUE) + 1,
                max_mito = 100,
                detect_doublets = FALSE,
                apply_filter = 1
            ))
            expect_s4_class(filtered_obj(), "Seurat")

            loaded(NULL)
            suppressWarnings(session$flushReact())
            expect_null(filtered_obj())

            loaded(obj)
            suppressWarnings(session$flushReact())
            expect_error(qc_result(), class = "shiny.silent.error")
        }
    )
})
