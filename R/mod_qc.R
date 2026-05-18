#' QC Module — UI
#'
#' @param id Module namespace ID.
#' @export
#' @keywords internal
mod_qc_ui <- function(id) {
    ns <- NS(id)

    step_card(
        title = "Step 2: Quality Control & Filtering",
        next_step = "After QC filtering, or choosing No filtering, proceed to Step 3: Normalization & PCA",

        layout_columns(
            col_widths = c(9, 3),
            card(
                card_header("QC Violin Plots"),
                card_body(
                    class = "plot-container plot-container-wide",
                    with_spinner(plotOutput(ns("vln_plot"), height = "420px", width = "100%"))
                )
            ),
            card(
                card_header("Filtering Parameters"),
                card_body(
                    radioButtons(
                        ns("filter_mode"),
                        "Filtering mode",
                        choices = list(
                            "No filtering" = "none",
                            "Use filtering thresholds" = "filter"
                        ),
                        selected = "none"
                    ),
                    checkboxInput(ns("detect_doublets"),
                                  "Detect and remove doublets (scDblFinder)",
                                  value = FALSE),
                    conditionalPanel(
                        condition = sprintf("input['%s'] == 'none'", ns("filter_mode")),
                        tags$p(
                            class = "text-muted small mb-0",
                            "Metric filtering is off. Normalization will use all loaded cells unless doublet removal is selected."
                        )
                    ),
                    conditionalPanel(
                        condition = sprintf("input['%s'] == 'filter'", ns("filter_mode")),
                        numericInput(ns("min_genes"),
                                     "Min genes per cell",
                                     value = 200, min = 0),
                        numericInput(ns("max_genes"),
                                     "Max genes per cell",
                                     value = 5000, min = 0),
                        numericInput(ns("max_mito"),
                                     "Max mitochondrial %",
                                     value = 5, min = 0, max = 100, step = 0.5),
                        tags$p(
                            class = "text-muted small",
                            "Suggested thresholds are calculated as median +/- 3 MADs from this dataset, capped to valid ranges."
                        )
                    ),
                    conditionalPanel(
                        condition = sprintf(
                            "input['%s'] == 'filter' || input['%s']",
                            ns("filter_mode"),
                            ns("detect_doublets")
                        ),
                        tags$hr(),
                        action_btn(ns("apply_filter"),
                                   "Apply Filters",
                                   icon = icon("filter"))
                    )
                )
            )
        ),

        conditionalPanel(
            condition = sprintf(
                "input['%s'] > 0 && (input['%s'] == 'filter' || input['%s'])",
                ns("apply_filter"),
                ns("filter_mode"),
                ns("detect_doublets")
            ),
            tags$hr(),
            layout_columns(
                col_widths = c(9, 3),
                card(
                    card_header("After Filtering"),
                    card_body(
                        class = "plot-container plot-container-wide",
                        with_spinner(plotOutput(ns("vln_plot_filtered"), height = "420px", width = "100%"))
                    )
                ),
                card(
                    card_header("Summary"),
                    card_body(
                        with_spinner(verbatimTextOutput(ns("filter_summary"))),
                        card_plot_download("qc", ns("download_qc_plot"), ns = ns)
                    )
                )
            )
        )
    )
}


#' QC Module — Server
#'
#' @param id Module namespace ID.
#' @param seurat_obj Reactive Seurat object from data loading module.
#' @return Reactive filtered Seurat object.
#' @export
#' @keywords internal
mod_qc_server <- function(id, seurat_obj) {
    moduleServer(id, function(input, output, session) {
        plot_width <- function(output_id, fallback = 1100) {
            width <- session$clientData[[sprintf(
                "output_%s_width",
                session$ns(output_id)
            )]]
            if (is.null(width) || !is.finite(width) || width < 200) {
                fallback
            } else {
                width
            }
        }

        filtered_obj <- reactiveVal(NULL)

        observeEvent(seurat_obj(), {
            defaults <- qc_filter_defaults(seurat_obj()@meta.data)
            updateNumericInput(session, "min_genes", value = defaults$min_genes)
            updateNumericInput(session, "max_genes", value = defaults$max_genes)
            updateNumericInput(session, "max_mito", value = defaults$max_mito)
            filtered_obj(NULL)
        }, ignoreNULL = TRUE)

        observeEvent(input$filter_mode, {
            filtered_obj(NULL)
        }, ignoreInit = TRUE)

        observeEvent(input$detect_doublets, {
            filtered_obj(NULL)
        }, ignoreInit = TRUE)

        observe({
            hideFeedback("max_mito")
            if (!is.null(input$max_mito) && (input$max_mito < 0 || input$max_mito > 100)) {
                showFeedbackDanger("max_mito", "Mitochondrial percentage must be between 0 and 100")
            }
        })

        output$vln_plot <- renderPlot({
            req(seurat_obj())
            scCustomize::VlnPlot_scCustom(
                seurat_object = seurat_obj(),
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                pt.size = 0.1
            )
        }, height = function() 420, width = function() plot_width("vln_plot"))

        observeEvent(input$apply_filter, {
            req(seurat_obj())

            withProgress(message = "Filtering cells...", value = 0.4, {
                obj <- seurat_obj()
                meta <- obj@meta.data
                keep <- rep(TRUE, nrow(meta))
                names(keep) <- rownames(meta)

                filter_mode <- input$filter_mode %||% "none"
                if (identical(filter_mode, "filter")) {
                    if ("nFeature_RNA" %in% colnames(meta)) {
                        keep <- keep &
                            meta$nFeature_RNA > input$min_genes &
                            meta$nFeature_RNA < input$max_genes
                    }
                    if ("percent.mt" %in% colnames(meta)) {
                        keep <- keep & meta$percent.mt < input$max_mito
                    }
                }
                keep[is.na(keep)] <- FALSE
                cells_to_keep <- names(keep)[keep]
                if (!length(cells_to_keep)) {
                    showNotification(
                        "No cells remain after applying these QC thresholds. Adjust the filters or choose No filtering.",
                        type = "error",
                        duration = 10
                    )
                    filtered_obj(NULL)
                    return(invisible(NULL))
                }
                obj <- subset(obj, cells = cells_to_keep)

                doublets_removed <- 0L
                if (isTRUE(input$detect_doublets)) {
                    setProgress(0.7, message = "Detecting doublets...")
                    obj <- tryCatch(
                        detect_doublets(
                            obj,
                            sample_col = if ("orig.ident" %in% colnames(obj@meta.data)) {
                                "orig.ident"
                            } else {
                                NULL
                            }
                        ),
                        error = function(e) {
                            showNotification(
                                paste0("Doublet detection failed: ", conditionMessage(e)),
                                type = "error",
                                duration = 12
                            )
                            NULL
                        }
                    )
                    req(!is.null(obj))

                    doublets_removed <- sum(obj$scDblFinder_class == "doublet", na.rm = TRUE)
                    singlets <- colnames(obj)[obj$scDblFinder_class == "singlet"]
                    if (!length(singlets)) {
                        showNotification(
                            "Doublet detection marked every remaining cell as a doublet. Adjust filters or disable doublet removal.",
                            type = "error",
                            duration = 10
                        )
                        filtered_obj(NULL)
                        return(invisible(NULL))
                    }
                    obj <- subset(obj, cells = singlets)
                    showNotification(
                        paste0("Removed ", doublets_removed, " doublets."),
                        type = "message",
                        duration = 6
                    )
                }

                attr(obj, "ascseurat_doublets_removed") <- doublets_removed
                attr(obj, "ascseurat_filter_mode") <- if (identical(filter_mode, "filter")) {
                    "filter"
                } else if (isTRUE(input$detect_doublets)) {
                    "doublet_only"
                } else {
                    "none"
                }
                setProgress(1, message = paste0(ncol(obj), " cells remaining"))
                filtered_obj(obj)
            })
        })

        qc_result <- reactive({
            req(seurat_obj())
            if (identical(input$filter_mode, "none")) {
                if (isTRUE(input$detect_doublets)) {
                    req(filtered_obj())
                    return(filtered_obj())
                }
                obj <- seurat_obj()
                attr(obj, "ascseurat_filter_mode") <- "none"
                attr(obj, "ascseurat_doublets_removed") <- 0L
                return(obj)
            }
            req(filtered_obj())
            filtered_obj()
        })

        output$vln_plot_filtered <- renderPlot({
            req(filtered_obj())
            scCustomize::VlnPlot_scCustom(
                seurat_object = filtered_obj(),
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                pt.size = 0.1
            )
        }, height = function() 420, width = function() plot_width("vln_plot_filtered"))

        output$filter_summary <- renderPrint({
            req(filtered_obj())
            before <- ncol(seurat_obj())
            after <- ncol(filtered_obj())
            removed <- before - after
            doublets_removed <- attr(filtered_obj(), "ascseurat_doublets_removed", exact = TRUE)
            if (is.null(doublets_removed)) {
                doublets_removed <- 0L
            }

            cat("Cells before filtering:", before, "\n")
            cat("Cells after filtering:", after, "\n")
            cat("Cells removed:", removed,
                sprintf("(%.1f%%)", (removed / before) * 100), "\n")
            if (isTRUE(input$detect_doublets)) {
                cat("Doublets removed:", doublets_removed, "\n")
            }
        })

        output$download_qc_plot <- plot_download_handler(
            plot_reactive = reactive({
                scCustomize::VlnPlot_scCustom(
                    seurat_object = filtered_obj(),
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    pt.size = 0.1
                )
            }),
            filename_base = "qc_violin",
            input = input,
            prefix = "qc"
        )

        qc_result
    })
}
