#' Advanced Plots Module — UI
#'
#' @param id Module namespace ID.
#' @export
#' @keywords internal
mod_advanced_plots_ui <- function(id) {
    ns <- NS(id)

    step_card(
        title = "Advanced Plots",

        layout_columns(
            col_widths = c(4, 4, 4),
            card(
                card_header("Data Source"),
                card_body(
                    radioButtons(ns("data_source"),
                                 "Use data from:",
                                 choices = list(
                                     "Upload RDS file" = "upload",
                                     "Single Sample tab" = "single",
                                     "Integration tab" = "integration"
                                 ),
                                 selected = "upload"),
                    conditionalPanel(
                        condition = sprintf("input['%s'] == 'upload'", ns("data_source")),
                        fileInput(ns("rds_file"),
                                  "Upload Seurat RDS:",
                                  accept = ".rds")
                    )
                )
            ),
            card(
                card_header("Plot Type"),
                card_body(
                    radioButtons(ns("plot_type"),
                                 "Select plot:",
                                 choices = list(
                                     "Stacked Violin Plot" = "stacked_vln",
                                     "Multi-Gene Dot Plot" = "dotplot"
                                 ),
                                 selected = "stacked_vln"),
                    input_markers_file(ns("gene_list")),
                    radioButtons(
                        ns("gene_count_mode"),
                        "Genes to plot:",
                        choices = list(
                            "All genes" = "all",
                            "First N genes" = "top_n"
                        ),
                        selected = "all"
                    ),
                    conditionalPanel(
                        condition = sprintf("input['%s'] == 'top_n'", ns("gene_count_mode")),
                        numericInput(ns("gene_count"),
                                     "Number of genes",
                                     value = 20,
                                     min = 1,
                                     step = 1)
                    ),
                    action_btn(ns("make_plot"), "Generate Plot", icon = icon("chart-bar"))
                )
            ),
            card(
                card_header("Download"),
                card_body(
                    downloadButton(ns("download_per_gene"),
                                   "Download per-gene bundle (zip)",
                                   class = "btn btn-outline-secondary w-100 mb-3"),
                    card_plot_download("adv", ns("download_adv"),
                                       height_default = 15, width_default = 20,
                                       ns = ns)
                )
            )
        ),

        uiOutput(ns("plot_and_ordering"))
    )
}


#' Advanced Plots Module — Server
#'
#' @param id Module namespace ID.
#' @param seurat_obj_single Reactive Seurat object from single-sample pipeline.
#' @param seurat_obj_multi Reactive Seurat object from integration pipeline.
#' @export
#' @keywords internal
mod_advanced_plots_server <- function(id, seurat_obj_single = NULL,
                                      seurat_obj_multi = NULL) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        setBookmarkExclude(c("rds_file", "gene_list"))

        gene_order <- reactiveVal(character())
        cluster_order <- reactiveVal(character())
        applied_gene_order <- reactiveVal(character())
        applied_cluster_order <- reactiveVal(character())

        load_advanced_object <- function() {
            current_obj <- function(obj_reactive, label) {
                if (is.null(obj_reactive)) {
                    return(NULL)
                }
                obj <- tryCatch(obj_reactive(), error = function(e) NULL)
                if (is.null(obj)) {
                    showNotification(
                        paste("No processed object is available from the", label, "tab yet."),
                        type = "warning",
                        duration = 8
                    )
                }
                obj
            }

            obj <- switch(
                input$data_source,
                upload = {
                    if (is.null(input$rds_file)) {
                        showNotification(
                            "Upload a Seurat RDS file or choose an object from another tab.",
                            type = "warning",
                            duration = 8
                        )
                        return(NULL)
                    }
                    load_seurat_rds(input$rds_file$datapath)
                },
                single = {
                    current_obj(seurat_obj_single, "Single Sample")
                },
                integration = {
                    current_obj(seurat_obj_multi, "Integration")
                }
            )

            if (is.null(obj)) {
                return(NULL)
            }
            if ("RNA" %in% names(obj@assays)) {
                SeuratObject::DefaultAssay(obj) <- "RNA"
            }
            obj
        }

        read_plot_genes <- function(obj) {
            if (is.null(input$gene_list)) {
                showNotification(
                    "Upload a marker gene list before generating an advanced plot.",
                    type = "warning",
                    duration = 8
                )
                return(character())
            }

            genes_df <- read_uploaded_table(
                input$gene_list$datapath,
                stringsAsFactors = FALSE
            )
            genes <- trimws(as.character(genes_df[[1]]))
            genes <- unique(genes[nzchar(genes)])
            genes <- genes[genes %in% rownames(obj)]

            if (identical(input$gene_count_mode %||% "all", "top_n")) {
                n_genes <- suppressWarnings(as.integer(input$gene_count))
                if (!is.finite(n_genes) || n_genes < 1) {
                    n_genes <- 1
                }
                genes <- head(genes, n_genes)
            }

            genes
        }

        current_identity_levels <- function(obj) {
            current_levels <- levels(SeuratObject::Idents(obj))
            if (is.null(current_levels)) {
                current_levels <- unique(as.character(SeuratObject::Idents(obj)))
            }
            as.character(current_levels)
        }

        plot_state <- eventReactive(input$make_plot, {
            obj <- load_advanced_object()
            if (is.null(obj)) {
                return(NULL)
            }

            genes <- tryCatch(
                read_plot_genes(obj),
                error = function(e) {
                    showNotification(
                        paste("Could not read the marker gene list:", conditionMessage(e)),
                        type = "error",
                        duration = 10
                    )
                    character()
                }
            )
            if (length(genes) == 0) {
                showNotification("No matching genes found in the dataset.",
                                 type = "warning")
                return(NULL)
            }

            gene_order(genes)
            identity_levels <- current_identity_levels(obj)
            cluster_order(identity_levels)
            applied_gene_order(genes)
            applied_cluster_order(identity_levels)

            list(
                obj = obj,
                genes = genes,
                plot_type = input$plot_type %||% "stacked_vln"
            )
        })

        selected_genes <- reactive({
            state <- plot_state()
            req(state)
            requested_order <- applied_gene_order()
            ordered_genes <- requested_order[requested_order %in% state$genes]
            c(ordered_genes, setdiff(state$genes, ordered_genes))
        })

        ordered_obj <- reactive({
            state <- plot_state()
            req(state)
            obj <- state$obj
            requested_order <- applied_cluster_order()

            if (length(requested_order)) {
                current_levels <- current_identity_levels(obj)
                ordered_levels <- requested_order[requested_order %in% current_levels]
                if (length(ordered_levels)) {
                    new_levels <- c(ordered_levels, setdiff(current_levels, ordered_levels))
                    SeuratObject::Idents(obj) <- factor(
                        as.character(SeuratObject::Idents(obj)),
                        levels = new_levels
                    )
                }
            }

            obj
        })

        move_order_value <- function(values, index, direction) {
            index <- suppressWarnings(as.integer(index))
            if (!length(values) || !is.finite(index)) {
                return(values)
            }
            target <- if (identical(direction, "up")) index - 1L else index + 1L
            if (index < 1L || index > length(values) ||
                target < 1L || target > length(values)) {
                return(values)
            }
            values[c(index, target)] <- values[c(target, index)]
            values
        }

        observeEvent(input$order_move, {
            move <- input$order_move
            if (is.null(move$type) || is.null(move$index) || is.null(move$direction)) {
                return(NULL)
            }

            if (identical(move$type, "cluster")) {
                cluster_order(move_order_value(cluster_order(), move$index, move$direction))
            } else {
                gene_order(move_order_value(gene_order(), move$index, move$direction))
            }
        })

        observeEvent(input$apply_order, {
            req(plot_state())
            applied_gene_order(gene_order())
            applied_cluster_order(cluster_order())
        })

        order_list_ui <- function(type, values) {
            values <- as.character(values)
            if (!length(values)) {
                return(tags$p(class = "text-muted small", "No items are available yet."))
            }

            tags$div(
                class = "asc-order-list",
                lapply(seq_along(values), function(i) {
                    tags$div(
                        class = "asc-order-row",
                        tags$span(class = "asc-order-label", title = values[[i]], values[[i]]),
                        tags$span(
                            class = "asc-order-controls",
                            tags$button(
                                type = "button",
                                class = "btn btn-sm btn-outline-secondary",
                                disabled = if (i == 1L) "disabled" else NULL,
                                onclick = sprintf(
                                    "Shiny.setInputValue('%s', {type: '%s', direction: 'up', index: %d, nonce: Math.random()}, {priority: 'event'})",
                                    ns("order_move"), type, i
                                ),
                                icon("arrow-up")
                            ),
                            tags$button(
                                type = "button",
                                class = "btn btn-sm btn-outline-secondary",
                                disabled = if (i == length(values)) "disabled" else NULL,
                                onclick = sprintf(
                                    "Shiny.setInputValue('%s', {type: '%s', direction: 'down', index: %d, nonce: Math.random()}, {priority: 'event'})",
                                    ns("order_move"), type, i
                                ),
                                icon("arrow-down")
                            )
                        )
                    )
                })
            )
        }

        output$ordering_ui <- renderUI({
            req(plot_state())
            card(
                card_header("Ordering"),
                card_body(
                    class = "asc-order-card",
                    tags$p(
                        class = "text-muted small",
                        "Use the arrows to adjust the order, then regenerate the plot."
                    ),
                    tags$h6("Genes"),
                    order_list_ui("gene", gene_order()),
                    tags$hr(),
                    tags$h6("Clusters"),
                    order_list_ui("cluster", cluster_order()),
                    tags$hr(),
                    action_btn(ns("apply_order"),
                               "Regenerate Plot",
                               icon = icon("rotate"))
                )
            )
        })

        output$plot_and_ordering <- renderUI({
            req(plot_state())
            tagList(
                tags$hr(),
                layout_columns(
                    col_widths = c(8, 4),
                    card(
                        card_header("Plot"),
                        card_body(
                            class = "plot-container",
                            with_spinner(uiOutput(ns("plot_output")))
                        )
                    ),
                    uiOutput(ns("ordering_ui"))
                )
            )
        })

        adv_plot <- reactive({
            state <- plot_state()
            req(state)
            obj <- tryCatch(
                ordered_obj(),
                error = function(e) {
                    showNotification(
                        "Choose or upload a Seurat object before generating an advanced plot.",
                        type = "warning",
                        duration = 8
                    )
                    NULL
                }
            )
            if (is.null(obj)) {
                return(NULL)
            }

            genes <- selected_genes()
            if (length(genes) == 0) {
                showNotification("No matching genes found in the dataset.",
                                 type = "warning")
                return(NULL)
            }

            if (identical(state$plot_type, "stacked_vln")) {
                scCustomize::Stacked_VlnPlot(
                    seurat_object = obj,
                    features = genes,
                    x_lab_rotate = TRUE
                )
            } else {
                scCustomize::DotPlot_scCustom(
                    seurat_object = obj,
                    features = genes,
                    flip_axes = TRUE
                )
            }
        })

        output$plot_output <- renderUI({
            n_genes <- tryCatch(length(selected_genes()), error = function(e) 10)
            height <- max(400, n_genes * 40)
            plotOutput(ns("adv_plot_render"), height = paste0(height, "px"))
        })

        output$adv_plot_render <- renderPlot({
            req(adv_plot())
            adv_plot()
        })

        output$download_adv <- plot_download_handler(
            plot_reactive = reactive(adv_plot()),
            filename_base = "advanced_plot",
            input = input,
            prefix = "adv"
        )

        output$download_per_gene <- downloadHandler(
            filename = function() paste0("advanced_per_gene_", Sys.Date(), ".zip"),
            content = function(file) {
                obj <- ordered_obj()
                genes <- selected_genes()
                req(length(genes) > 0)

                tmpdir <- tempfile("ascseurat_advanced_per_gene_")
                dir.create(tmpdir)
                on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

                withProgress(message = "Generating per-gene plots...",
                             min = 0, max = length(genes), value = 0, {
                    for (i in seq_along(genes)) {
                        gene <- genes[[i]]
                        dot_plot <- scCustomize::DotPlot_scCustom(
                            seurat_object = obj,
                            features = gene,
                            flip_axes = TRUE
                        )
                        violin_plot <- scCustomize::VlnPlot_scCustom(
                            seurat_object = obj,
                            features = gene,
                            pt.size = 0.1
                        )

                        ggplot2::ggsave(
                            file.path(tmpdir, paste0("dotplot_", gene, ".png")),
                            dot_plot,
                            width = 12,
                            height = 10,
                            units = "cm",
                            dpi = 300
                        )
                        ggplot2::ggsave(
                            file.path(tmpdir, paste0("violin_", gene, ".png")),
                            violin_plot,
                            width = 12,
                            height = 10,
                            units = "cm",
                            dpi = 300
                        )
                        setProgress(i, message = paste0("Plotting ", gene, " (",
                                                        i, "/", length(genes), ")"))
                    }
                })

                utils::zip(
                    zipfile = file,
                    files = list.files(tmpdir, full.names = TRUE),
                    flags = "-j"
                )
            },
            contentType = "application/zip"
        )
    })
}
