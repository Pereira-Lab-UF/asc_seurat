#' Automated Cell-Type Annotation Module - UI
#'
#' @param id Module namespace ID.
#' @export
#' @keywords internal
mod_celltype_annotation_ui <- function(id) {
    ns <- NS(id)

    step_card(
        title = "Automated Cell-Type Annotation (SingleR)",
        layout_columns(
            col_widths = c(4, 4, 4),
            card(
                card_header("Data"),
                card_body(
                    radioButtons(
                        ns("data_source"),
                        "Use data from:",
                        choices = list(
                            "Upload processed Seurat RDS" = "upload",
                            "Single Sample tab" = "single",
                            "Integration tab" = "integration"
                        ),
                        selected = "upload"
                    ),
                    conditionalPanel(
                        condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
                        fileInput(ns("rds_file"),
                                  "Load processed Seurat RDS:",
                                  accept = ".rds")
                    ),
                    tags$p(class = "text-muted small",
                           "Choose a processed Seurat object with normalized expression.")
                )
            ),
            card(
                card_header("Reference"),
                card_body(
                    selectInput(
                        ns("reference"),
                        "Reference atlas:",
                        choices = list(
                            "Human Primary Cell Atlas (HPCA)" = "hpca",
                            "Blueprint/ENCODE (human)" = "blueprint",
                            "Monaco Immune (human, immune-only)" = "monaco",
                            "Mouse RNA-seq (Tabula Muris)" = "mousernaseq",
                            "ImmGen (mouse, immune-only)" = "immgen"
                        ),
                        selected = "hpca"
                    ),
                    radioButtons(
                        ns("granularity"),
                        "Annotation granularity:",
                        choices = list(
                            "Main labels" = "label.main",
                            "Fine labels" = "label.fine"
                        ),
                        selected = "label.main"
                    )
                )
            ),
            card(
                card_header("Run"),
                card_body(
                    action_btn(ns("run_singler"),
                               "Annotate Cells",
                               icon = icon("tags"))
                )
            )
        ),

        conditionalPanel(
            condition = sprintf("input['%s'] > 0", ns("run_singler")),
            tags$hr(),
            layout_columns(
                col_widths = c(6, 6),
                card(
                    card_header("UMAP - SingleR labels"),
                    card_body(
                        class = "plot-container",
                        with_spinner(plotOutput(ns("singler_umap"), height = "500px")),
                        card_plot_download(
                            "singler_umap",
                            ns("download_singler_umap"),
                            ns = ns,
                            height_default = 10,
                            width_default = 12
                        )
                    )
                ),
                card(
                    card_header("Score heatmap"),
                    card_body(
                        class = "plot-container",
                        with_spinner(plotOutput(ns("singler_scores"), height = "500px")),
                        card_plot_download(
                            "singler_scores",
                            ns("download_singler_scores"),
                            ns = ns,
                            height_default = 12,
                            width_default = 16
                        )
                    )
                )
            ),
            card(
                card_header("Per-cluster annotations"),
                card_body(
                    with_spinner(reactableOutput(ns("singler_table"))),
                    tags$br(),
                    layout_columns(
                        col_widths = c(6, 6),
                        downloadButton(
                            ns("download_singler_table"),
                            "Download Annotation Table (CSV)",
                            class = "btn btn-outline-secondary w-100"
                        ),
                        downloadButton(
                            ns("download_annotated"),
                            "Download annotated Seurat object",
                            class = "btn btn-outline-secondary w-100"
                        )
                    )
                )
            )
        )
    )
}


#' Prepare a Seurat Object for SingleR
#'
#' Joins split Seurat v5 assay layers and returns an expression matrix.
#'
#' @param obj A Seurat object.
#' @return A list with `obj`, `expr_matrix`, and `assay`.
#' @noRd
prepare_singler_input <- function(obj) {
    if (!inherits(obj, "Seurat")) {
        stop("SingleR requires a Seurat object.", call. = FALSE)
    }

    assay <- if ("RNA" %in% names(obj@assays)) {
        "RNA"
    } else {
        SeuratObject::DefaultAssay(obj)
    }
    SeuratObject::DefaultAssay(obj) <- assay

    layers <- tryCatch(
        SeuratObject::Layers(obj[[assay]]),
        error = function(e) character()
    )
    split_layers <- grep("^(counts|data)($|[.])", layers, value = TRUE)
    if (length(split_layers) > 2L || any(grepl("^(counts|data)[.]", split_layers))) {
        obj <- tryCatch(
            SeuratObject::JoinLayers(obj, assay = assay),
            error = function(e) {
                stop(
                    "This integrated Seurat object has split assay layers that could not be joined for SingleR: ",
                    conditionMessage(e),
                    call. = FALSE
                )
            }
        )
    }

    expr_matrix <- tryCatch(
        Seurat::GetAssayData(obj, assay = assay, layer = "data"),
        error = function(e) NULL
    )

    if (is.null(expr_matrix) || !nrow(expr_matrix) || !ncol(expr_matrix)) {
        obj <- tryCatch(
            Seurat::NormalizeData(obj, assay = assay, verbose = FALSE),
            error = function(e) obj
        )
        expr_matrix <- tryCatch(
            Seurat::GetAssayData(obj, assay = assay, layer = "data"),
            error = function(e) NULL
        )
    }

    if (is.null(expr_matrix) || !nrow(expr_matrix) || !ncol(expr_matrix)) {
        expr_matrix <- tryCatch(
            Seurat::GetAssayData(obj, assay = assay, layer = "counts"),
            error = function(e) {
                stop(
                    "Could not retrieve expression values for SingleR from assay '",
                    assay,
                    "'.",
                    call. = FALSE
                )
            }
        )
    }

    if (!nrow(expr_matrix) || !ncol(expr_matrix)) {
        stop("SingleR expression matrix is empty.", call. = FALSE)
    }

    list(obj = obj, expr_matrix = expr_matrix, assay = assay)
}

#' Automated Cell-Type Annotation Module - Server
#'
#' @param id Module namespace ID.
#' @param seurat_obj_single Optional reactive Seurat object from the single-sample workflow.
#' @param seurat_obj_multi Optional reactive Seurat object from the integration workflow.
#' @export
#' @keywords internal
mod_celltype_annotation_server <- function(id,
                                           seurat_obj_single = NULL,
                                           seurat_obj_multi = NULL) {
    moduleServer(id, function(input, output, session) {
        setBookmarkExclude(c("rds_file"))
        annotated_obj <- reactiveVal(NULL)
        singler_result <- reactiveVal(NULL)

        get_reference <- function(reference_name) {
            if (!requireNamespace("celldex", quietly = TRUE)) {
                stop(
                    "Package 'celldex' is not installed. Install with: BiocManager::install('celldex')",
                    call. = FALSE
                )
            }

            cache_dir <- file.path(tools::R_user_dir("ascseurat", "cache"), "celldex")
            dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
            cache_file <- file.path(cache_dir, paste0(reference_name, ".rds"))
            if (file.exists(cache_file)) {
                return(readRDS(cache_file))
            }

            ref <- switch(
                reference_name,
                hpca = celldex::HumanPrimaryCellAtlasData(),
                blueprint = celldex::BlueprintEncodeData(),
                monaco = celldex::MonacoImmuneData(),
                mousernaseq = celldex::MouseRNAseqData(),
                immgen = celldex::ImmGenData()
            )
            saveRDS(ref, cache_file)
            ref
        }

        get_input_object <- function() {
            switch(
                input$data_source,
                upload = {
                    req(input$rds_file)
                    load_seurat_rds(input$rds_file$datapath)
                },
                single = {
                    req(seurat_obj_single, seurat_obj_single())
                    seurat_obj_single()
                },
                integration = {
                    req(seurat_obj_multi, seurat_obj_multi())
                    seurat_obj_multi()
                }
            )
        }

        observeEvent(input$run_singler, {
            if (!requireNamespace("SingleR", quietly = TRUE)) {
                showNotification(
                    "Package 'SingleR' is not installed. Install with: BiocManager::install('SingleR')",
                    type = "error",
                    duration = 12
                )
                return(NULL)
            }

            withProgress(message = "Running SingleR...", value = 0.2, {
                tryCatch({
                    obj <- get_input_object()
                    ref <- get_reference(input$reference)
                    labels <- ref[[input$granularity]]
                    if (all(is.na(labels))) {
                        labels <- ref$label.main
                    }

                    setProgress(0.45, message = "Preparing expression matrix...")
                    singler_input <- prepare_singler_input(obj)
                    obj <- singler_input$obj
                    expr_matrix <- singler_input$expr_matrix

                    setProgress(0.65, message = "Annotating cells...")
                    pred <- SingleR::SingleR(
                        test = expr_matrix,
                        ref = ref,
                        labels = labels
                    )

                    obj$singler_label <- pred$labels
                    annotated_obj(obj)
                    singler_result(pred)
                    setProgress(1, message = "SingleR annotation complete.")
                }, error = function(e) {
                    showNotification(
                        paste0("SingleR annotation failed: ", conditionMessage(e)),
                        type = "error",
                        duration = 15
                    )
                })
            })
        })

        singler_umap_plot <- reactive({
            req(annotated_obj())
            scCustomize::DimPlot_scCustom(
                seurat_object = annotated_obj(),
                group.by = "singler_label",
                label = TRUE,
                figure_plot = TRUE
            )
        })

        singler_scores_plot <- reactive({
            req(singler_result())
            function() {
                SingleR::plotScoreHeatmap(singler_result())
            }
        })

        singler_annotation_table <- reactive({
            req(annotated_obj())
            cluster_labels <- as.character(Idents(annotated_obj()))
            singler_labels <- as.character(annotated_obj()$singler_label)
            tab <- table(cluster = cluster_labels, label = singler_labels)
            modal_label <- apply(tab, 1, function(values) {
                names(values)[which.max(values)]
            })
            modal_cells <- apply(tab, 1, max)

            df <- data.frame(
                Cluster = names(modal_label),
                Label = unname(modal_label),
                Cells = as.integer(modal_cells),
                stringsAsFactors = FALSE
            )
            df
        })

        output$singler_umap <- renderPlot({
            singler_umap_plot()
        })

        output$singler_scores <- renderPlot({
            singler_scores_plot()()
        })

        output$download_singler_umap <- plot_download_handler(
            plot_reactive = singler_umap_plot,
            filename_base = "singler_umap",
            input = input,
            prefix = "singler_umap"
        )

        output$download_singler_scores <- plot_download_handler(
            plot_reactive = singler_scores_plot,
            filename_base = "singler_score_heatmap",
            input = input,
            prefix = "singler_scores"
        )

        output$singler_table <- renderReactable({
            reactable(
                singler_annotation_table(),
                striped = TRUE,
                highlight = TRUE,
                defaultPageSize = 20
            )
        })

        output$download_singler_table <- downloadHandler(
            filename = function() paste0("singler_annotations_", Sys.Date(), ".csv"),
            content = function(file) {
                write.csv(singler_annotation_table(), file, row.names = FALSE)
            }
        )

        output$download_annotated <- downloadHandler(
            filename = function() paste0("annotated_singler_", Sys.Date(), ".rds"),
            content = function(file) {
                saveRDS(annotated_obj(), file)
            }
        )

        list(
            annotated_obj = reactive({
                req(annotated_obj())
                annotated_obj()
            }),
            singler_result = reactive({
                req(singler_result())
                singler_result()
            })
        )
    })
}
