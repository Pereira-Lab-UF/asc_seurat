#' Data Loading Module — UI
#'
#' @param id Module namespace ID.
#' @param demo_mode Logical. When TRUE, hides input controls (used by Demo tab).
#' @export
#' @keywords internal
mod_data_loading_ui <- function(id, demo_mode = FALSE) {
    ns <- NS(id)
    current_workdir <- ascseurat_workdir()

    input_controls <- if (demo_mode) {
        NULL
    } else {
        list(
            layout_columns(
                col_widths = c(4, 4, 4),
                card(
                    card_body(
                        radioButtons(ns("input_type"),
                                     "Select input type:",
                                     choices = list(
                                         "10X Genomics directory" = "10x_dir",
                                         "10X h5 file" = "10x_h5",
                                         "CSV count matrix" = "csv",
                                         "TSV count matrix" = "tsv",
                                         "h5ad (AnnData)" = "h5ad",
                                         "Existing Seurat RDS" = "rds"
                                     ),
                                     selected = "10x_dir"),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == '10x_dir'", ns("input_type")),
                            textInput(ns("data_path"),
                                      "Path to 10X data directory:",
                                      value = "data/",
                                      placeholder = "e.g., /path/to/sample/filtered_feature_bc_matrix"),
                            data_path_help(current_workdir)
                        )
                    )
                ),
                card(
                    card_body(
                        textInput(ns("project_name"),
                                  label = "Project name",
                                  value = "",
                                  placeholder = "Optional, e.g. AscSeurat"),
                        input_mito_regex(ns("mito_pattern")),
                        conditionalPanel(
                            condition = sprintf("input['%s'] != 'rds'", ns("input_type")),
                            numericInput(ns("min_cells"),
                                         "Min cells per gene", value = 3, min = 0),
                            numericInput(ns("min_features"),
                                         "Min genes per cell", value = 200, min = 0)
                        )
                    )
                ),
                card(
                    card_body(
                        conditionalPanel(
                            condition = sprintf("input['%s'] != '10x_dir'", ns("input_type")),
                            fileInput(ns("data_file"),
                                      "Upload file:",
                                      accept = c(".h5", ".h5ad", ".csv", ".tsv", ".rds"))
                        ),
                        tags$br(),
                        action_btn(ns("load_data"), "Load Data", icon = icon("upload"))
                    )
                )
            ),
            tags$details(
                class = "asc-collapse-panel",
                tags$summary(
                    class = "asc-collapse-summary",
                    icon("table", class = "me-2"),
                    "Add cell-level metadata (optional)"
                ),
                tags$div(
                    class = "asc-collapse-body",
                    fileInput(ns("metadata_file"),
                              "Upload metadata CSV/TSV:",
                              accept = c(".csv", ".tsv")),
                    tags$p(class = "text-muted small",
                           "First column must be cell barcode. Other columns will be added ",
                           "to the Seurat object metadata."),
                    action_btn(ns("apply_metadata"),
                               "Add Metadata",
                               icon = icon("table"))
                )
            )
        )
    }

    summary_condition <- if (demo_mode) {
        "true"
    } else {
        sprintf(
            "input['%s'] > 0 || input['%s'] > 0",
            ns("load_data"), ns("apply_metadata")
        )
    }

    step_card(
        title = "Step 1: Load Data",
        next_step = "After loading, proceed to Step 2: Quality Control & Filtering",
        input_controls,
        conditionalPanel(
            condition = summary_condition,
            tags$hr(),
            layout_columns(
                col_widths = c(6, 6),
                card(
                    card_header("Data Summary"),
                    card_body(
                        with_spinner(verbatimTextOutput(ns("data_summary")))
                    )
                ),
                card(
                    card_header("Mitochondrial Gene Targets"),
                    card_body(
                        with_spinner(verbatimTextOutput(ns("mito_genes")))
                    )
                )
            )
        )
    )
}


#' Data Loading Module — Server
#'
#' @param id Module namespace ID.
#' @param trigger_demo Optional reactive value. When incremented above 0, demo
#'   data loads automatically (used by the Demo tab).
#' @return Reactive Seurat object.
#' @export
#' @keywords internal
mod_data_loading_server <- function(id, trigger_demo = NULL) {
    moduleServer(id, function(input, output, session) {
        setBookmarkExclude(c("data_file", "metadata_file"))
        state_obj <- reactiveVal(NULL)

        project_name <- function() {
            project <- trimws(as.character(input$project_name %||% ""))
            if (nzchar(project)) project else "AscSeurat"
        }

        finalize_loaded_object <- function(obj, mito_pattern) {
            mito_pattern <- trimws(as.character(mito_pattern %||% ""))
            if (!"percent.mt" %in% colnames(obj@meta.data)) {
                obj[["percent.mt"]] <- if (nzchar(mito_pattern)) {
                    PercentageFeatureSet(obj, pattern = mito_pattern)
                } else {
                    0
                }
            }
            attr(obj, "ascseurat_mito_pattern") <- mito_pattern
            obj
        }

        output$data_ready <- shiny::renderText({
            if (is.null(state_obj())) "false" else "true"
        })
        shiny::outputOptions(output, "data_ready", suspendWhenHidden = FALSE)

        do_load_demo <- function() {
            withProgress(message = "Loading demo data...", value = 0.4, {
                demo_path <- system.file("extdata", "pbmc_demo.rds", package = "ascseurat")
                if (demo_path == "") {
                    demo_path <- file.path(getwd(), "inst", "extdata", "pbmc_demo.rds")
                }
                if (!file.exists(demo_path)) {
                    showNotification(
                        "Demo dataset not found in the package installation.",
                        type = "error",
                        duration = 10
                    )
                    return(NULL)
                }
                obj <- load_seurat_rds(demo_path)
                obj <- finalize_loaded_object(obj, isolate(input$mito_pattern))
                state_obj(obj)
                setProgress(1, message = "Demo data ready.")
            })
        }

        if (!is.null(trigger_demo)) {
            observeEvent(trigger_demo(), {
                req(trigger_demo() > 0)
                do_load_demo()
            }, ignoreInit = TRUE)
        }

        observeEvent(input$load_data, {
            req(input$input_type)

            withProgress(message = "Loading data...", value = 0.3, {
                path <- if (input$input_type == "10x_dir") {
                    req(input$data_path)
                    path <- resolve_input_path(input$data_path)
                    if (!dir.exists(path)) {
                        showNotification(
                            paste0(
                                "10X data directory not found: ", path,
                                ". Relative paths start from: ", ascseurat_workdir(),
                                ". In Docker, confirm that the host data folder was mounted and use a container-visible path such as data/sample_name."
                            ),
                            type = "error",
                            duration = 14
                        )
                        return(NULL)
                    }
                    if (!identical(detect_input_format(path), "10x_dir")) {
                        showNotification(
                            paste0(
                                "Directory exists but does not look like a 10X matrix folder: ",
                                path,
                                ". Expected matrix.mtx(.gz), barcodes.tsv(.gz), and genes.tsv/features.tsv(.gz)."
                            ),
                            type = "error",
                            duration = 14
                        )
                        return(NULL)
                    }
                    path
                } else {
                    req(input$data_file)
                    file_size <- file.info(input$data_file$datapath)$size
                    if (is.finite(file_size) && file_size > 500 * 1024^2) {
                        showNotification(
                            "Large uploads can take several minutes and may exceed available memory.",
                            type = "warning",
                            duration = 10
                        )
                    }
                    input$data_file$datapath
                }

                obj <- tryCatch(
                    load_input_data(
                        path = path,
                        project = project_name(),
                        min.cells = if (input$input_type != "rds") input$min_cells else 3,
                        min.features = if (input$input_type != "rds") input$min_features else 200
                    ),
                    error = function(e) {
                        showNotification(
                            paste0("Data loading failed: ", conditionMessage(e)),
                            type = "error",
                            duration = 12
                        )
                        NULL
                    }
                )
                req(!is.null(obj))

                setProgress(0.7, message = "Computing mitochondrial content...")
                obj <- finalize_loaded_object(obj, isolate(input$mito_pattern))
                state_obj(obj)
                setProgress(1, message = "Done!")
            })
        })

        observeEvent(input$load_demo, { do_load_demo() })

        observeEvent(input$apply_metadata, {
            req(state_obj(), input$metadata_file)

            metadata_df <- read_uploaded_table(
                input$metadata_file$datapath,
                stringsAsFactors = FALSE
            )

            if (ncol(metadata_df) < 2) {
                showNotification(
                    "Metadata file must have a barcode column and at least one annotation column.",
                    type = "error",
                    duration = 10
                )
                return(NULL)
            }

            barcodes <- as.character(metadata_df[[1]])
            metadata_df <- metadata_df[, -1, drop = FALSE]
            rownames(metadata_df) <- barcodes

            obj <- state_obj()
            overlap <- intersect(colnames(obj), rownames(metadata_df))
            if (!length(overlap)) {
                showNotification(
                    "No barcodes from the metadata file matched the loaded Seurat object.",
                    type = "error",
                    duration = 10
                )
                return(NULL)
            }

            for (column in colnames(metadata_df)) {
                values <- metadata_df[colnames(obj), column]
                obj[[column]] <- values
            }

            state_obj(obj)
            showNotification(
                paste0(
                    "Added ", ncol(metadata_df), " metadata columns to ",
                    length(overlap), " matching cells. ",
                    nrow(metadata_df) - length(overlap),
                    " metadata rows did not match."
                ),
                type = "message",
                duration = 6
            )
        })

        output$data_summary <- renderPrint({
            req(state_obj())
            obj <- state_obj()
            cat("Cells:", ncol(obj), "\n")
            cat("Features:", nrow(obj), "\n")
            cat("Samples:", length(unique(obj$orig.ident)), "\n")
            cat("Assay:", SeuratObject::DefaultAssay(obj), "\n")
        })

        output$mito_genes <- renderPrint({
            req(state_obj())
            obj <- state_obj()
            pattern <- attr(obj, "ascseurat_mito_pattern", exact = TRUE) %||%
                input$mito_pattern
            pattern <- trimws(as.character(pattern %||% ""))
            if (!nzchar(pattern)) {
                cat("No mitochondrial gene pattern set.")
                return(invisible(NULL))
            }
            mito_genes <- grep(pattern, rownames(obj), value = TRUE)
            if (length(mito_genes) == 0) {
                cat("No mitochondrial genes found with pattern:", pattern)
            } else {
                cat("Found", length(mito_genes), "mitochondrial genes:\n")
                cat(paste(mito_genes, collapse = ", "))
            }
        })

        reactive({
            req(state_obj())
            state_obj()
        })
    })
}
