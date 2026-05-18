expression_layer_label <- function(layer) {
    switch(
        layer,
        "counts" = "counts",
        "data" = "normalized",
        "scale.data" = "scaled",
        layer
    )
}

ensure_expression_layer <- function(obj, genes, layer = "data") {
    genes <- intersect(unique(genes), rownames(obj))
    if (!identical(layer, "scale.data") || !length(genes)) {
        return(obj)
    }

    assay <- SeuratObject::DefaultAssay(obj)
    scaled_features <- tryCatch(
        rownames(SeuratObject::LayerData(obj, assay = assay, layer = "scale.data")),
        error = function(e) character()
    )
    missing_features <- setdiff(genes, scaled_features)
    if (length(missing_features)) {
        obj <- suppressWarnings(Seurat::ScaleData(
            obj,
            assay = assay,
            features = genes,
            verbose = FALSE
        ))
    }
    obj
}

fetch_expression_layer <- function(obj, genes, layer = "data") {
    genes <- intersect(unique(genes), rownames(obj))
    if (!length(genes)) {
        stop("None of the selected genes are present in the Seurat object.", call. = FALSE)
    }
    obj <- ensure_expression_layer(obj, genes, layer)

    data <- Seurat::FetchData(
        object = obj,
        vars = genes,
        layer = layer,
        clean = FALSE
    )
    genes <- intersect(genes, colnames(data))
    if (!length(genes)) {
        stop(
            "The selected genes are not available in the ",
            expression_layer_label(layer),
            " expression layer.",
            call. = FALSE
        )
    }

    mat <- as.matrix(data[, genes, drop = FALSE])
    storage.mode(mat) <- "numeric"
    mat
}

clustered_expression_dot_plot <- function(obj, genes, layer = "data") {
    obj <- ensure_expression_layer(obj, genes, layer)
    expr_mat <- fetch_expression_layer(obj, genes, layer)
    genes <- colnames(expr_mat)

    assay_layers <- tryCatch(
        names(obj[[SeuratObject::DefaultAssay(obj)]]@layers),
        error = function(e) character()
    )
    pct_layer <- if ("counts" %in% assay_layers) {
        "counts"
    } else {
        "data"
    }
    pct_mat_source <- tryCatch(
        fetch_expression_layer(obj, genes, pct_layer),
        error = function(e) expr_mat
    )
    pct_mat_source <- pct_mat_source[rownames(expr_mat), genes, drop = FALSE]

    idents <- SeuratObject::Idents(obj)
    clusters <- as.character(idents[rownames(expr_mat)])
    clusters[is.na(clusters) | !nzchar(clusters)] <- "Unassigned"
    cluster_levels <- levels(idents)
    cluster_levels <- cluster_levels[cluster_levels %in% clusters]
    if (!length(cluster_levels)) {
        cluster_levels <- sort(unique(clusters))
    }

    avg_mat <- vapply(
        cluster_levels,
        function(cluster) {
            colMeans(expr_mat[clusters == cluster, , drop = FALSE], na.rm = TRUE)
        },
        numeric(length(genes))
    )
    rownames(avg_mat) <- genes

    pct_mat <- vapply(
        cluster_levels,
        function(cluster) {
            colMeans(pct_mat_source[clusters == cluster, , drop = FALSE] > 0,
                     na.rm = TRUE) * 100
        },
        numeric(length(genes))
    )
    rownames(pct_mat) <- genes

    order_mat <- avg_mat
    order_mat[!is.finite(order_mat)] <- 0
    gene_order <- rownames(order_mat)
    cluster_order <- colnames(order_mat)
    if (nrow(order_mat) > 1L) {
        gene_order <- rownames(order_mat)[stats::hclust(stats::dist(order_mat))$order]
    }
    if (ncol(order_mat) > 1L) {
        cluster_order <- colnames(order_mat)[stats::hclust(stats::dist(t(order_mat)))$order]
    }

    plot_df <- data.frame(
        gene = rep(rownames(avg_mat), times = ncol(avg_mat)),
        cluster = rep(colnames(avg_mat), each = nrow(avg_mat)),
        avg_exp = as.vector(avg_mat),
        pct_exp = as.vector(pct_mat),
        stringsAsFactors = FALSE
    )
    plot_df$gene <- factor(plot_df$gene, levels = rev(gene_order))
    plot_df$cluster <- factor(plot_df$cluster, levels = cluster_order)

    plot <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = .data$cluster, y = .data$gene)
    ) +
        ggplot2::geom_point(
            ggplot2::aes(size = .data$pct_exp, color = .data$avg_exp),
            alpha = 0.9
        ) +
        ggplot2::scale_size(
            range = c(1.5, 8),
            limits = c(0, 100),
            name = "Percent expressing"
        ) +
        ggplot2::labs(
            x = "Cluster",
            y = NULL,
            color = paste0("Mean ", expression_layer_label(layer), " expression")
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_line(color = "#e6e0d2", linewidth = 0.3),
            panel.grid.minor = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.position = "right"
        )

    if (identical(layer, "scale.data")) {
        plot + ggplot2::scale_color_gradient2(
            low = "#2b6cb0",
            mid = "#f7f7f7",
            high = "#b83232",
            midpoint = 0,
            na.value = "grey85"
        )
    } else {
        plot + ggplot2::scale_color_gradient(
            low = "#f0efe7",
            high = "#087f5b",
            na.value = "grey85"
        )
    }
}

#' Visualization Module — UI
#'
#' @param id Module namespace ID.
#' @param step_number Numeric step number to show in the workflow title.
#' @export
#' @keywords internal
mod_visualization_ui <- function(id, step_number = 6) {
    ns <- NS(id)
    if (is.null(step_number) || !is.finite(step_number)) {
        step_number <- 6L
    }
    step_number <- as.integer(step_number)
    marker_source_choices <- list("Upload file" = "upload")
    marker_source_choices[[paste0(
        "Use de novo marker genes detected in Step ",
        step_number - 1L
    )]] <- "de_results"

    step_card(
        title = paste0("Step ", step_number, ": Gene Expression Visualization"),
        next_step = "Analysis complete! Use Cell-type Annotation or Advanced Plots for deeper exploration.",

        # Gene input
        layout_columns(
            col_widths = c(3, 6, 3),

            # Source selector
            card(
                card_header("Source"),
                card_body(
                    radioButtons(
                        ns("marker_source"), NULL,
                        choices = marker_source_choices,
                        selected = "upload"
                    )
                )
            ),

            # Gene controls (adapts to source)
            card(
                card_header("Gene controls"),
                card_body(
                    # Upload mode
                    conditionalPanel(
                        condition = sprintf("input['%s'] == 'upload'", ns("marker_source")),
                        input_markers_file(ns("markers_file")),
                        selectInput(ns("header_opt"), "File has header?",
                                    choices = c("Yes", "No"), selected = "Yes"),
                        action_btn(ns("load_markers"), "Load Markers", icon = icon("upload")),
                        tags$small(class = "text-muted d-block mt-1",
                                   "Format: one gene per row; optional 2nd column for group."),
                        conditionalPanel(
                            condition = sprintf("input['%s'] > 0", ns("load_markers")),
                            tags$hr(),
                            with_spinner(uiOutput(ns("group_select"))),
                            with_spinner(uiOutput(ns("upload_gene_count_ui"))),
                            with_spinner(uiOutput(ns("gene_select")))
                        )
                    ),
                    # DE results mode
                    conditionalPanel(
                        condition = sprintf("input['%s'] == 'de_results'", ns("marker_source")),
                        with_spinner(uiOutput(ns("de_cluster_select"))),
                        numericInput(
                            ns("top_n_genes"),
                            "Top N genes (ranked by adj. p-value)",
                            value = 10, min = 1, max = 500, step = 5
                        ),
                        uiOutput(ns("de_genes_preview"))
                    )
                )
            ),

            # Visualize card
            card(
                card_header("Visualize"),
                card_body(
                    radioButtons(ns("slot_selection"),
                                 "Expression values:",
                                 choices = list("Counts" = "counts",
                                                "Normalized" = "data",
                                                 "Scaled" = "scale.data"),
                                 selected = "data"),
                    tags$hr(),
                    action_btn(ns("show_expression"), "Show Expression Plots",
                               icon = icon("chart-area")),
                    uiOutput(ns("plot_generation_status")),
                    uiOutput(ns("plot_refresh_prompt"))
                )
            )
        ),

        # Clustered dot plot
        conditionalPanel(
            condition = sprintf("input['%s'] > 0", ns("show_expression")),
            tags$hr(),
            layout_columns(
                col_widths = c(9, 3),
                card(
                    card_header("Clustered Dot Plot"),
                    card_body(
                        class = "plot-container plot-container-wide",
                        with_spinner(uiOutput(ns("clustered_dot_plot_ui")))
                    )
                ),
                card_plot_download("hm", ns("download_heatmap"),
                                   height_default = 15, width_default = 20,
                                   ns = ns)
            )
        ),

        # Feature + violin plots
        conditionalPanel(
            condition = sprintf("input['%s'] > 0", ns("show_expression")),
            tags$hr(),
            card(
                card_header("Feature Plots"),
                card_body(
                    class = "plot-container plot-container-wide",
                    with_spinner(uiOutput(ns("feature_plots")))
                )
            ),
            card(
                card_header("Violin Plots"),
                card_body(
                    class = "plot-container plot-container-wide",
                    with_spinner(uiOutput(ns("violin_plots")))
                )
            ),
            layout_columns(
                col_widths = c(8, 4),
                card_plot_download("feat", ns("download_features"),
                                   height_default = 10, width_default = 14,
                                   ns = ns),
                card(
                    card_header("Per-gene bundle"),
                    card_body(
                        downloadButton(ns("download_per_gene"),
                                       "Download per-gene bundle (zip)",
                                       class = "btn-download-soft w-100")
                    )
                )
            )
        )
    )
}


#' Visualization Module — Server
#'
#' @param id Module namespace ID.
#' @param seurat_obj Reactive clustered Seurat object.
#' @param markers_de Optional reactive returning the markers data frame from
#'   the DE module (Step 5). When supplied, users can select genes directly
#'   from DE results instead of uploading a file.
#' @param de_step_number Numeric DE step number to use in user messages.
#' @export
#' @keywords internal
mod_visualization_server <- function(id,
                                     seurat_obj,
                                     markers_de = NULL,
                                     de_step_number = 5) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        setBookmarkExclude(c("markers_file"))
        de_step_label <- paste0("Step ", de_step_number)
        plot_request <- reactiveVal(NULL)

        plot_width <- function(output_id, fallback = 1100) {
            width <- session$clientData[[sprintf(
                "output_%s_width",
                session$ns(output_id)
            )]]
            if (is.null(width) || !is.finite(width) || width < 250) {
                fallback
            } else {
                width
            }
        }

        expression_obj <- reactive({
            req(seurat_obj())
            obj <- seurat_obj()
            if ("RNA" %in% names(obj@assays)) {
                SeuratObject::DefaultAssay(obj) <- "RNA"
            }
            obj
        })

        # Load marker gene list (upload path)
        markers_data <- eventReactive(input$load_markers, {
            req(input$markers_file)
            has_header <- input$header_opt == "Yes"
            df <- read_uploaded_table(
                input$markers_file$datapath,
                header = has_header,
                stringsAsFactors = FALSE
            )
            df[[1]] <- trimws(as.character(df[[1]]))
            df <- df[nzchar(df[[1]]), , drop = FALSE]
            rownames(df) <- NULL
            df
        })

        available_uploaded_genes <- reactive({
            req(markers_data())
            genes <- markers_data()[[1]]
            if (!is.null(input$selected_group) && input$selected_group != "All" &&
                ncol(markers_data()) >= 2) {
                group_values <- trimws(as.character(markers_data()[[2]]))
                genes <- markers_data()[[1]][group_values == input$selected_group]
            }
            genes <- unique(genes[nzchar(genes)])
            if (!is.null(expression_obj())) {
                genes <- intersect(genes, rownames(expression_obj()))
            }
            genes
        })

        output$group_select <- renderUI({
            req(markers_data())
            if (ncol(markers_data()) >= 2) {
                groups <- unique(trimws(as.character(markers_data()[[2]])))
                groups <- groups[!is.na(groups) & nzchar(groups)]
                selectInput(ns("selected_group"), "Gene group:",
                            choices = c("All", groups), selected = "All")
            }
        })

        output$upload_gene_count_ui <- renderUI({
            genes <- available_uploaded_genes()
            n_genes <- length(genes)
            req(n_genes > 0)
            numericInput(
                ns("upload_gene_count"),
                "Number of genes to select",
                value = n_genes,
                min = 1,
                max = n_genes,
                step = 1
            )
        })

        output$gene_select <- renderUI({
            genes <- available_uploaded_genes()
            n_select <- input$upload_gene_count %||% length(genes)
            n_select <- suppressWarnings(as.integer(n_select))
            if (!is.finite(n_select)) {
                n_select <- length(genes)
            }
            n_select <- max(1L, min(length(genes), n_select))
            genes <- head(genes, n_select)
            checkboxGroupInput(ns("selected_genes"), "Select genes:",
                               choices = genes, selected = genes, inline = TRUE)
        })

        # DE results path — cluster selector
        output$de_cluster_select <- renderUI({
            if (is.null(markers_de)) {
                return(tags$p(class = "text-muted small",
                              paste0("Run ", de_step_label, " first to use DE results here.")))
            }
            df <- tryCatch(markers_de(), error = function(e) NULL)
            if (is.null(df) || nrow(df) == 0) {
                return(tags$p(class = "text-muted small",
                              paste0("No markers yet - run ", de_step_label, " first.")))
            }
            if ("cluster" %in% colnames(df)) {
                clusters <- sort(unique(as.character(df$cluster)))
                selectInput(ns("de_cluster"), "Cluster:", choices = clusters)
            } else {
                tags$p(class = "text-muted small",
                       "Results are not split by cluster - all markers will be used.")
            }
        })

        # DE results path — preview of selected genes
        output$de_genes_preview <- renderUI({
            if (is.null(markers_de)) return(NULL)
            df <- tryCatch(markers_de(), error = function(e) NULL)
            if (is.null(df) || nrow(df) == 0) return(NULL)

            if ("cluster" %in% colnames(df) && !is.null(input$de_cluster)) {
                df <- df[as.character(df$cluster) == input$de_cluster, , drop = FALSE]
            }
            if ("p_val_adj" %in% colnames(df)) {
                df <- df[order(df$p_val_adj), , drop = FALSE]
            }
            n <- min(input$top_n_genes %||% 10, nrow(df))
            genes <- df$gene[seq_len(n)]
            if (!is.null(expression_obj())) {
                genes <- intersect(genes, rownames(expression_obj()))
            }
            req(length(genes) > 0)
            tags$div(
                class = "mt-2",
                tags$p(class = "text-muted small mb-1",
                       paste0(length(genes), " genes selected:")),
                tags$p(class = "small font-monospace text-break",
                       paste(genes, collapse = ", "))
            )
        })

        # Active genes — unified reactive for both paths
        active_genes <- reactive({
            if (isTRUE(input$marker_source == "de_results")) {
                req(!is.null(markers_de))
                df <- tryCatch(markers_de(), error = function(e) NULL)
                req(!is.null(df), nrow(df) > 0)

                if ("cluster" %in% colnames(df) && !is.null(input$de_cluster)) {
                    df <- df[as.character(df$cluster) == input$de_cluster, , drop = FALSE]
                }
                if ("p_val_adj" %in% colnames(df)) {
                    df <- df[order(df$p_val_adj), , drop = FALSE]
                }
                n <- min(input$top_n_genes %||% 10, nrow(df))
                genes <- df$gene[seq_len(n)]
                if (!is.null(expression_obj())) {
                    genes <- intersect(genes, rownames(expression_obj()))
                }
                genes
            } else {
                genes <- input$selected_genes
                if (is.null(genes) || !length(genes)) {
                    genes <- available_uploaded_genes()
                    n_select <- input$upload_gene_count %||% length(genes)
                    n_select <- suppressWarnings(as.integer(n_select))
                    if (!is.finite(n_select)) {
                        n_select <- length(genes)
                    }
                    n_select <- max(1L, min(length(genes), n_select))
                    genes <- head(genes, n_select)
                }
                req(length(genes) > 0)
                genes
            }
        })

        current_expression_selection <- reactive({
            req(expression_obj(), active_genes())
            genes <- intersect(active_genes(), rownames(expression_obj()))
            req(length(genes) > 0)
            list(
                genes = genes,
                slot = input$slot_selection %||% "data"
            )
        })

        requested_genes <- reactive({
            req(plot_request())
            plot_request()$genes
        })

        requested_slot <- reactive({
            req(plot_request())
            plot_request()$slot
        })

        observeEvent(input$show_expression, {
            selection <- current_expression_selection()
            plot_request(selection)
            showNotification(
                "Generating expression plots. Please wait while the plots render below.",
                type = "message",
                duration = 8
            )
        }, ignoreInit = TRUE)

        observe({
            req((input$show_expression %||% 0) > 0)
            if (!is.null(plot_request())) {
                return(NULL)
            }
            plot_request(current_expression_selection())
        })

        output$plot_generation_status <- renderUI({
            show_count <- max(
                input$show_expression %||% 0,
                input$show_heatmap %||% 0,
                input$show_features %||% 0
            )
            if (show_count <= 0) {
                return(NULL)
            }

            req(plot_request())
            genes <- requested_genes()
            req(length(genes) > 0)

            tags$div(
                class = "alert alert-info small mt-3 mb-0",
                icon("circle-info", class = "me-1"),
                paste0(
                    "Generating expression plots for ",
                    length(genes),
                    " genes. Please wait for the plots below to finish rendering."
                )
            )
        })

        output$plot_refresh_prompt <- renderUI({
            request <- plot_request()
            if (is.null(request)) {
                return(NULL)
            }
            current <- tryCatch(current_expression_selection(), error = function(e) NULL)
            if (is.null(current)) {
                return(NULL)
            }
            changed <- !identical(current$genes, request$genes) ||
                !identical(current$slot, request$slot)
            if (!changed) {
                return(NULL)
            }

            tags$div(
                class = "alert alert-warning small mt-3 mb-0",
                icon("circle-info", class = "me-1"),
                "Selections changed. Click \"Show Expression Plots\" again to update the plots."
            )
        })

        clustered_dot_plot <- reactive({
            req(expression_obj(), requested_genes(), requested_slot())
            obj <- expression_obj()
            genes <- intersect(requested_genes(), rownames(obj))
            req(length(genes) > 0)
            withProgress(message = "Generating clustered dot plot...", value = 0.25, {
                setProgress(0.55, detail = paste(length(genes), "genes"))
                plot <- clustered_expression_dot_plot(
                    obj = obj,
                    genes = genes,
                    layer = requested_slot()
                )
                setProgress(1)
                plot
            })
        })

        output$clustered_dot_plot_ui <- renderUI({
            req(expression_obj(), requested_genes())
            genes <- intersect(requested_genes(), rownames(expression_obj()))
            height <- max(450, ceiling(length(genes) / 2) * 110)
            plotOutput(
                ns("clustered_dot_plot_render"),
                height = paste0(height, "px"),
                width = "100%"
            )
        })

        draw_clustered_dot_plot <- function() {
            plot_obj <- clustered_dot_plot()
            print(plot_obj)
        }

        output$clustered_dot_plot_render <- renderPlot({
            draw_clustered_dot_plot()
        }, width = function() plot_width("clustered_dot_plot_render"))

        # Feature plots
        output$feature_plots <- renderUI({
            req(expression_obj(), requested_genes())
            genes <- intersect(requested_genes(), rownames(expression_obj()))
            n_genes <- length(genes)
            height <- max(350, ceiling(n_genes / 4) * 300)
            plotOutput(
                ns("feature_plot_render"),
                height = paste0(height, "px"),
                width = "100%"
            )
        })

        feature_plot <- reactive({
            req(expression_obj(), requested_genes())
            genes <- intersect(requested_genes(), rownames(expression_obj()))
            req(length(genes) > 0)
            withProgress(message = "Generating feature plots...", value = 0.25, {
                setProgress(0.55, detail = paste(length(genes), "genes"))
                plot_obj <- ensure_expression_layer(
                    expression_obj(),
                    genes,
                    requested_slot()
                )
                plot <- scCustomize::FeaturePlot_scCustom(
                    seurat_object = plot_obj,
                    features = genes,
                    num_columns = min(4, max(1, length(genes))),
                    slot = requested_slot(),
                    na_cutoff = if (identical(requested_slot(), "scale.data")) {
                        NULL
                    } else {
                        1e-09
                    }
                )
                setProgress(1)
                plot
            })
        })

        output$feature_plot_render <- renderPlot({
            feature_plot()
        }, width = function() plot_width("feature_plot_render"))

        # Violin plots
        output$violin_plots <- renderUI({
            req(expression_obj(), requested_genes())
            genes <- intersect(requested_genes(), rownames(expression_obj()))
            n_genes <- length(genes)
            height <- max(350, ceiling(n_genes / 4) * 300)
            plotOutput(
                ns("violin_plot_render"),
                height = paste0(height, "px"),
                width = "100%"
            )
        })

        violin_plot <- reactive({
            req(expression_obj(), requested_genes())
            genes <- intersect(requested_genes(), rownames(expression_obj()))
            req(length(genes) > 0)
            withProgress(message = "Generating violin plots...", value = 0.25, {
                setProgress(0.55, detail = paste(length(genes), "genes"))
                plot_obj <- ensure_expression_layer(
                    expression_obj(),
                    genes,
                    requested_slot()
                )
                plot <- scCustomize::VlnPlot_scCustom(
                    seurat_object = plot_obj,
                    features = genes,
                    pt.size = 0.1,
                    layer = requested_slot()
                )
                setProgress(1)
                plot
            })
        })

        output$violin_plot_render <- renderPlot({
            violin_plot()
        }, width = function() plot_width("violin_plot_render"))

        output$download_heatmap <- plot_download_handler(
            plot_reactive = reactive(draw_clustered_dot_plot),
            filename_base = "clustered_dot_plot",
            input = input,
            prefix = "hm"
        )

        output$download_features <- plot_download_handler(
            plot_reactive = reactive({
                patchwork::wrap_plots(
                    feature_plot(),
                    violin_plot(),
                    ncol = 1
                )
            }),
            filename_base = "feature_and_violin_plots",
            input = input,
            prefix = "feat"
        )

        output$download_per_gene <- downloadHandler(
            filename = function() paste0("per_gene_plots_", Sys.Date(), ".zip"),
            content = function(file) {
                req(expression_obj(), requested_genes())
                genes <- intersect(requested_genes(), rownames(expression_obj()))
                req(length(genes) > 0)

                tmpdir <- tempfile("ascseurat_per_gene_")
                dir.create(tmpdir)
                on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

                withProgress(message = "Generating per-gene plots...",
                             min = 0, max = length(genes), value = 0, {
                    plot_obj <- ensure_expression_layer(
                        expression_obj(),
                        genes,
                        requested_slot()
                    )
                    for (i in seq_along(genes)) {
                        gene <- genes[[i]]
                        feature_single <- scCustomize::FeaturePlot_scCustom(
                            seurat_object = plot_obj,
                            features = gene,
                            slot = requested_slot(),
                            na_cutoff = if (identical(requested_slot(), "scale.data")) {
                                NULL
                            } else {
                                1e-09
                            }
                        )
                        violin_single <- scCustomize::VlnPlot_scCustom(
                            seurat_object = plot_obj,
                            features = gene,
                            pt.size = 0.1,
                            layer = requested_slot()
                        )

                        ggplot2::ggsave(
                            file.path(tmpdir, paste0("feature_", gene, ".png")),
                            feature_single,
                            width = 12,
                            height = 10,
                            units = "cm",
                            dpi = 300
                        )
                        ggplot2::ggsave(
                            file.path(tmpdir, paste0("violin_", gene, ".png")),
                            violin_single,
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
