capture_notifications <- function(session) {
    captured <- new.env()
    captured$messages <- character()
    root <- session$rootScope()
    root$sendNotification <- function(type, message) {
        if (identical(type, "show")) {
            captured$messages <- c(captured$messages, as.character(message$html))
        }
    }
    captured
}

capture_input_updates <- function(session) {
    captured <- new.env()
    captured$updates <- list()
    root <- session$rootScope()
    root$sendInputMessage <- function(inputId, message) {
        captured$updates[[length(captured$updates) + 1]] <- list(
            id = sub("^.*-", "", inputId),
            message = message
        )
    }
    captured
}

updates_for <- function(captured, id) {
    Filter(function(update) identical(update$id, id), captured$updates)
}

lognorm_inputs <- list(
    norm_method = "lognorm",
    n_variable_features = 500,
    scale_factor = 10000,
    score_cc = FALSE,
    regress_cc = FALSE
)

test_that("Normalization ignores cell-cycle regression once scoring is unticked", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_normalization_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(do.call(session$setInputs, lognorm_inputs))
            suppressWarnings(session$setInputs(score_cc = TRUE, regress_cc = TRUE))
            suppressWarnings(session$setInputs(score_cc = FALSE))
            suppressMessages(suppressWarnings(session$setInputs(run_norm = 1)))

            expect_false(session$isClosed())
            expect_s4_class(normalized_obj(), "Seurat")
            expect_null(normalized_obj()@commands[["ScaleData.RNA"]]$vars.to.regress)
        }
    )
})

test_that("LogNormalize ignores regression variables chosen for SCTransform", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_normalization_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(do.call(session$setInputs, lognorm_inputs))
            suppressWarnings(session$setInputs(
                norm_method = "sctransform",
                vars_to_regress = "percent.mt"
            ))
            suppressWarnings(session$setInputs(norm_method = "lognorm"))
            suppressMessages(suppressWarnings(session$setInputs(run_norm = 1)))

            expect_s4_class(normalized_obj(), "Seurat")
            expect_null(normalized_obj()@commands[["ScaleData.RNA"]]$vars.to.regress)
        }
    )
})

test_that("Visible regression choices are still passed to run_normalization", {
    obj <- demo_fixture()
    calls <- list()
    local_mocked_bindings(run_normalization = function(obj, ...) {
        calls[[length(calls) + 1]] <<- list(...)
        obj
    })

    shiny::testServer(
        ascseurat::mod_normalization_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(do.call(session$setInputs, lognorm_inputs))
            suppressWarnings(session$setInputs(
                norm_method = "sctransform",
                vars_to_regress = "percent.mt",
                score_cc = TRUE,
                regress_cc = TRUE
            ))
            suppressWarnings(session$setInputs(run_norm = 1))
        }
    )

    expect_length(calls, 1)
    expect_identical(calls[[1]]$method, "sctransform")
    expect_identical(calls[[1]]$vars_to_regress, "percent.mt")
    expect_true(calls[[1]]$score_cell_cycle)
    expect_true(calls[[1]]$regress_cc)
})

test_that("Normalization errors show a notification instead of closing the session", {
    obj <- demo_fixture()
    local_mocked_bindings(run_normalization = function(obj, ...) {
        stop("None of the requested variables to regress are present")
    })

    shiny::testServer(
        ascseurat::mod_normalization_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            captured <- capture_notifications(session)
            suppressWarnings(do.call(session$setInputs, lognorm_inputs))
            suppressWarnings(session$setInputs(run_norm = 1))

            expect_false(session$isClosed())
            expect_null(normalized_obj())
            expect_length(captured$messages, 1)
            expect_match(captured$messages, "Normalization failed: None of the requested variables")
        }
    )
})

test_that("SCTransform ignores the hidden LogNormalize scale factor", {
    obj <- demo_fixture()
    calls <- list()
    local_mocked_bindings(run_normalization = function(obj, ...) {
        calls[[length(calls) + 1]] <<- list(...)
        obj
    })

    shiny::testServer(
        ascseurat::mod_normalization_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(do.call(session$setInputs, lognorm_inputs))
            # Cleared under LogNormalize, then hidden by choosing SCTransform.
            suppressWarnings(session$setInputs(scale_factor = NA))
            suppressWarnings(session$setInputs(norm_method = "sctransform"))
            suppressWarnings(session$setInputs(run_norm = 1))
            suppressWarnings(session$setInputs(norm_method = "lognorm", scale_factor = 1e6))
            suppressWarnings(session$setInputs(norm_method = "sctransform"))
            suppressWarnings(session$setInputs(run_norm = 2))
            # A visible LogNormalize scale factor is still used.
            suppressWarnings(session$setInputs(norm_method = "lognorm", scale_factor = 5e4))
            suppressWarnings(session$setInputs(run_norm = 3))
        }
    )

    expect_length(calls, 3)
    expect_identical(calls[[1]]$method, "sctransform")
    expect_equal(calls[[1]]$scale_factor, 10000)
    expect_equal(calls[[2]]$scale_factor, 10000)
    expect_identical(calls[[3]]$method, "lognorm")
    expect_equal(calls[[3]]$scale_factor, 5e4)
})

test_that("Blank or non-positive normalization settings are refused with a notification", {
    obj <- demo_fixture()
    calls <- 0L
    local_mocked_bindings(run_normalization = function(obj, ...) {
        calls <<- calls + 1L
        obj
    })

    shiny::testServer(
        ascseurat::mod_normalization_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            captured <- capture_notifications(session)
            suppressWarnings(do.call(session$setInputs, lognorm_inputs))
            suppressWarnings(session$setInputs(scale_factor = NA, run_norm = 1))
            suppressWarnings(session$setInputs(scale_factor = 0, run_norm = 2))
            suppressWarnings(session$setInputs(
                scale_factor = 10000,
                n_variable_features = NA,
                run_norm = 3
            ))
            suppressWarnings(session$setInputs(norm_method = "sctransform", run_norm = 4))

            expect_false(session$isClosed())
            expect_null(normalized_obj())
            expect_length(captured$messages, 4)
            expect_match(captured$messages[1:2], "positive value for Scale factor")
            expect_match(captured$messages[3:4], "positive value for Number of variable features")
            expect_no_match(captured$messages[4], "Scale factor")
        }
    )

    expect_identical(calls, 0L)
})

test_that("Hidden normalization controls are reset so logged inputs match the run", {
    obj <- demo_fixture()

    shiny::testServer(
        ascseurat::mod_normalization_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            captured <- capture_input_updates(session)
            suppressWarnings(do.call(session$setInputs, lognorm_inputs))
            captured$updates <- list()
            suppressWarnings(session$setInputs(score_cc = TRUE, regress_cc = TRUE))
            expect_length(updates_for(captured, "regress_cc"), 0)

            suppressWarnings(session$setInputs(score_cc = FALSE))
            regress_updates <- updates_for(captured, "regress_cc")
            expect_length(regress_updates, 1)
            expect_false(regress_updates[[1]]$message$value)

            captured$updates <- list()
            suppressWarnings(session$setInputs(
                norm_method = "sctransform",
                vars_to_regress = "percent.mt"
            ))
            expect_length(updates_for(captured, "vars_to_regress"), 0)
            scale_updates <- updates_for(captured, "scale_factor")
            expect_length(scale_updates, 1)
            expect_equal(as.numeric(scale_updates[[1]]$message$value), 10000)

            suppressWarnings(session$setInputs(norm_method = "lognorm"))
            vars_updates <- updates_for(captured, "vars_to_regress")
            expect_length(vars_updates, 1)
            expect_length(vars_updates[[1]]$message$value, 0)
        }
    )
})

test_that("Normalized object is dropped when the QC result becomes unavailable", {
    obj <- demo_fixture()
    local_mocked_bindings(run_normalization = function(obj, ...) obj)
    qc_obj <- shiny::reactiveVal(obj)

    shiny::testServer(
        ascseurat::mod_normalization_server,
        args = list(seurat_obj = shiny::reactive({
            shiny::req(qc_obj())
            qc_obj()
        })),
        {
            suppressWarnings(do.call(session$setInputs, lognorm_inputs))
            suppressWarnings(session$setInputs(run_norm = 1))
            expect_s4_class(normalized_obj(), "Seurat")

            qc_obj(NULL)
            session$flushReact()
            expect_null(normalized_obj())
            expect_error(session$getReturned()$obj(), class = "shiny.silent.error")
        }
    )
})

test_that("Changing QC without Apply hides the old normalized object", {
    obj <- demo_fixture()
    local_mocked_bindings(run_normalization = function(obj, ...) obj)
    loaded <- shiny::reactiveVal(obj)

    chained_server <- function(input, output, session) {
        qc <- ascseurat::mod_qc_server("qc", seurat_obj = shiny::reactive({
            shiny::req(loaded())
            loaded()
        }))
        norm <- ascseurat::mod_normalization_server("norm", seurat_obj = qc)
    }
    norm_available <- function(norm) {
        !is.null(tryCatch(norm$obj(), error = function(e) NULL))
    }

    shiny::testServer(chained_server, {
        suppressWarnings(session$setInputs(
            `qc-filter_mode` = "none",
            `qc-detect_doublets` = FALSE,
            `norm-norm_method` = "lognorm"
        ))
        suppressWarnings(session$setInputs(`norm-run_norm` = 1))
        expect_true(norm_available(norm))

        # Switching to filtering thresholds without pressing Apply.
        suppressWarnings(session$setInputs(`qc-filter_mode` = "filter"))
        expect_false(norm_available(norm))

        suppressWarnings(session$setInputs(`qc-filter_mode` = "none"))
        suppressWarnings(session$setInputs(`norm-run_norm` = 2))
        expect_true(norm_available(norm))

        # Ticking doublet removal without pressing Apply.
        suppressWarnings(session$setInputs(`qc-detect_doublets` = TRUE))
        expect_false(norm_available(norm))

        # Loading a new dataset while filtering thresholds are in use.
        suppressWarnings(session$setInputs(
            `qc-detect_doublets` = FALSE,
            `qc-filter_mode` = "filter",
            `qc-min_genes` = 0,
            `qc-max_genes` = 1e6,
            `qc-max_mito` = 100
        ))
        suppressWarnings(session$setInputs(`qc-apply_filter` = 1))
        suppressWarnings(session$setInputs(`norm-run_norm` = 3))
        expect_true(norm_available(norm))

        loaded(obj[, seq_len(ncol(obj) - 100)])
        suppressWarnings(session$flushReact())
        expect_false(norm_available(norm))
    })
})
