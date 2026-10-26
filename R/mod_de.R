#' Differential Expression Module — UI
#'
#' @param id Module namespace ID.
#' @param allow_conserved Logical. Whether to expose conserved-marker analysis.
#' @param step_number Numeric step number to show in the workflow title.
#' @export
#' @keywords internal
mod_de_ui <- function(id, allow_conserved = FALSE, step_number = 5) {
    ns <- NS(id)
    if (is.null(step_number) || !is.finite(step_number)) {
        step_number <- 5L
    }
    step_number <- as.integer(step_number)

    analysis_choices <- list(
        "Find markers for all clusters" = "all",
        "Find markers for one cluster" = "one",
        "Compare two clusters" = "compare",
        "DE between conditions" = "conditions"
    )
    if (isTRUE(allow_conserved)) {
        analysis_choices[["Find conserved markers across samples"]] <- "conserved"
    }

    step_card(
        title = paste0(
            "Step ", step_number,
            ": Differential Expression / Marker Identification (Optional)"
        ),
        next_step = paste0(
            "Proceed to Step ", step_number + 1L,
            " - or skip this step and upload your own marker list there directly."
        ),

        layout_columns(
            col_widths = c(3, 3, 3, 3),
            card(
                class = "de-panel-card h-100",
                card_body(
                    class = "d-flex flex-column h-100",
                    radioButtons(
                        ns("de_type"),
                        "Analysis type:",
                        choices = analysis_choices,
                        selected = "all"
                    )
                )
            ),
            div(
                class = "de-selection-panel h-100",
                conditionalPanel(
                    condition = sprintf(
                        "input['%s'] == 'one' || input['%s'] == 'compare' || input['%s'] == 'conditions' || input['%s'] == 'conserved'",
                        ns("de_type"), ns("de_type"), ns("de_type"), ns("de_type")
                    ),
                    card(
                        class = "de-panel-card de-selection-card h-100",
                        card_body(
                            class = "d-flex flex-column h-100",
                            with_spinner(uiOutput(ns("cluster_select_1"))),
                            conditionalPanel(
                                condition = sprintf("input['%s'] == 'compare'", ns("de_type")),
                                with_spinner(uiOutput(ns("cluster_select_2")))
                            ),
                            conditionalPanel(
                                condition = sprintf("input['%s'] == 'conditions'", ns("de_type")),
                                with_spinner(uiOutput(ns("condition_column_ui"))),
                                with_spinner(uiOutput(ns("condition_value_1_ui"))),
                                with_spinner(uiOutput(ns("condition_value_2_ui")))
                            ),
                            conditionalPanel(
                                condition = sprintf("input['%s'] == 'conserved'", ns("de_type")),
                                with_spinner(uiOutput(ns("conserved_grouping_ui")))
                            )
                        )
                    )
                )
            ),
            card(
                class = "de-panel-card h-100",
                card_body(
                    class = "d-flex flex-column h-100",
                    input_de_test(ns("test_use")),
                    numericInput(
                        ns("pval_cutoff"),
                        "Adjusted p-value cutoff",
                        value = 0.05, min = 0, max = 1, step = 0.01
                    ),
                    numericInput(
                        ns("logfc_threshold"),
                        "Log2FC threshold",
                        value = 0.25, min = 0, step = 0.05
                    ),
                    numericInput(
                        ns("min_pct"),
                        "Min fraction of cells expressing gene",
                        value = 0.1, min = 0, max = 1, step = 0.05
                    ),
                    checkboxInput(
                        ns("only_pos"),
                        "Only return positive markers (upregulated in cluster)",
                        value = TRUE
                    )
                )
            ),
            card(
                class = "de-panel-card h-100",
                card_body(
                    class = "d-flex flex-column h-100",
                    action_btn(ns("run_de"), "Find Markers", icon = icon("search")),
                    tags$hr(),
                    conditionalPanel(
                        condition = sprintf("input['%s'] > 0", ns("run_de")),
                        downloadButton(
                            ns("download_markers"),
                            "Download Results (CSV)",
                            class = "btn btn-outline-secondary w-100"
                        )
                    )
                )
            )
        ),

        conditionalPanel(
            condition = sprintf("input['%s'] > 0", ns("run_de")),
            tags$hr(),
            card(
                card_header("Marker Genes"),
                card_body(
                    with_spinner(reactableOutput(ns("markers_table")))
                )
            )
        )
    )
}


#' Differential Expression Module — Server
#'
#' @param id Module namespace ID.
#' @param seurat_obj Reactive clustered Seurat object.
#' @param allow_conserved Logical. Whether to expose conserved-marker analysis.
#' @export
#' @keywords internal
mod_de_server <- function(id, seurat_obj, allow_conserved = FALSE) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        format_marker_results <- function(result, cluster = NULL) {
            if (!"gene" %in% colnames(result)) {
                result$gene <- rownames(result)
            }

            if (!"cluster" %in% colnames(result) && !is.null(cluster)) {
                result$cluster <- cluster
            }

            rownames(result) <- NULL

            leading_cols <- intersect(c("gene", "cluster"), colnames(result))
            result[, c(leading_cols, setdiff(colnames(result), leading_cols)), drop = FALSE]
        }

        sample_metadata_col <- function(obj) {
            candidates <- c("samples", "sample_name", "orig.ident")
            found <- candidates[candidates %in% colnames(obj@meta.data)]
            if (length(found)) found[[1]] else NULL
        }

        # A cleared numericInput arrives as NA.
        pval_cutoff_ok <- function(value) {
            !is.null(value) && !is.na(value) && value > 0 && value < 1
        }

        min_pct_ok <- function(value) {
            !is.null(value) && !is.na(value) && value >= 0 && value <= 1
        }

        logfc_ok <- function(value) {
            !is.null(value) && !is.na(value) && value >= 0
        }

        observe({
            hideFeedback("pval_cutoff")
            if (!is.null(input$pval_cutoff) && !pval_cutoff_ok(input$pval_cutoff)) {
                showFeedbackDanger("pval_cutoff", "P-value must be between 0 and 1")
            }
        })

        observe({
            hideFeedback("min_pct")
            if (!is.null(input$min_pct) && !min_pct_ok(input$min_pct)) {
                showFeedbackDanger("min_pct", "Minimum percent must be between 0 and 1")
            }
        })

        output$cluster_select_1 <- renderUI({
            req(seurat_obj())
            clusters <- levels(Idents(seurat_obj()))
            label <- if (identical(input$de_type, "compare")) "Cluster 1:" else "Cluster:"
            selectInput(ns("cluster_1"), label, choices = clusters)
        })

        output$cluster_select_2 <- renderUI({
            req(seurat_obj())
            clusters <- levels(Idents(seurat_obj()))
            selectInput(ns("cluster_2"), "Cluster 2:", choices = clusters, multiple = TRUE)
        })

        condition_columns <- reactive({
            req(seurat_obj())
            obj <- seurat_obj()
            meta <- obj@meta.data

            candidate_cols <- names(meta)[vapply(meta, function(column) {
                column <- column[!is.na(column)]
                # Continuous values (QC percentages, scores) are not groups;
                # whole-number codes such as time points or batches can be.
                if (is.numeric(column) && any(column != round(column))) {
                    return(FALSE)
                }
                values <- unique(as.character(column))
                length(values) >= 2 && length(values) <= min(50, ncol(obj))
            }, logical(1))]

            candidate_cols <- setdiff(candidate_cols, c(
                "seurat_clusters",
                "celltype.treat",
                "ascseurat_condition"
            ))
            # Cluster assignments and per-cell count QC are not conditions.
            candidate_cols <- grep(
                "(_snn$|_snn_res\\.|^nCount_|^nFeature_)",
                candidate_cols,
                value = TRUE,
                invert = TRUE
            )

            preferred <- c("samples", "sample_name", "treat", "orig.ident")
            c(intersect(preferred, candidate_cols), setdiff(candidate_cols, preferred))
        })

        output$condition_column_ui <- renderUI({
            cols <- condition_columns()
            if (!length(cols)) {
                return(tags$p(
                    class = "text-muted small",
                    "No metadata column with two or more groups is available to compare."
                ))
            }

            sample_col <- sample_metadata_col(seurat_obj())
            default_col <- if (!is.null(sample_col) && sample_col %in% cols) {
                sample_col
            } else {
                cols[[1]]
            }

            labels <- cols
            labels[labels %in% sample_col] <- "samples"
            choices <- stats::setNames(cols, labels)

            selectInput(
                ns("condition_column"),
                "Metadata column:",
                choices = choices,
                selected = default_col
            )
        })

        output$conserved_grouping_ui <- renderUI({
            req(isTRUE(allow_conserved))
            cols <- condition_columns()
            req(length(cols) > 0)

            default_col <- if ("orig.ident" %in% cols) "orig.ident" else cols[[1]]
            selectInput(
                ns("conserved_grouping"),
                "Sample column:",
                choices = cols,
                selected = default_col
            )
        })

        condition_values <- reactive({
            req(seurat_obj(), input$condition_column)
            obj <- seurat_obj()
            values <- obj[[input$condition_column]][, 1]
            if (!is.null(input$cluster_1)) {
                values <- values[as.character(Seurat::Idents(obj)) == input$cluster_1]
            }
            values <- sort(unique(as.character(values[!is.na(values)])))
            req(length(values) >= 2)
            values
        })

        output$condition_value_1_ui <- renderUI({
            values <- condition_values()
            selectInput(
                ns("condition_value_1"),
                "Condition 1:",
                choices = values,
                selected = values[[1]]
            )
        })

        output$condition_value_2_ui <- renderUI({
            values <- condition_values()
            default_value <- if (length(values) >= 2) values[[2]] else values[[1]]
            selectInput(
                ns("condition_value_2"),
                "Condition 2:",
                choices = values,
                selected = default_value
            )
        })

        de_result <- reactiveVal(NULL)

        # Results from an earlier object (before re-clustering, renaming or
        # subsetting) no longer match its clusters.
        observeEvent(tryCatch(seurat_obj(), error = function(e) NULL), {
            de_result(NULL)
        }, ignoreNULL = FALSE, ignoreInit = TRUE)

        observeEvent(input$run_de, {
            de_result(NULL)
            req(seurat_obj())
            if (!pval_cutoff_ok(input$pval_cutoff)) {
                showNotification(
                    "Differential expression was not run: set the adjusted p-value cutoff to a number between 0 and 1.",
                    type = "warning",
                    duration = 12
                )
                return(invisible(NULL))
            }
            if (!min_pct_ok(input$min_pct)) {
                showNotification(
                    "Differential expression was not run: set the minimum fraction of cells expressing the gene to a number between 0 and 1.",
                    type = "warning",
                    duration = 12
                )
                return(invisible(NULL))
            }
            if (!logfc_ok(input$logfc_threshold)) {
                showNotification(
                    "Differential expression was not run: set the log2 fold-change threshold to a number of 0 or more.",
                    type = "warning",
                    duration = 12
                )
                return(invisible(NULL))
            }

            withProgress(message = "Finding markers...", value = 0.3, {
                obj <- seurat_obj()
                if ("RNA" %in% names(obj@assays)) {
                    SeuratObject::DefaultAssay(obj) <- "RNA"
                }
                obj <- tryCatch(
                    SeuratObject::JoinLayers(obj),
                    error = function(e) obj
                )

                result <- tryCatch(
                    switch(input$de_type,
                        "all" = FindAllMarkers(
                            obj,
                            test.use = input$test_use,
                            logfc.threshold = input$logfc_threshold,
                            min.pct = input$min_pct,
                            only.pos = isTRUE(input$only_pos),
                            verbose = FALSE
                        ),
                        "one" = {
                            req(input$cluster_1)
                            result <- FindMarkers(
                                obj,
                                ident.1 = input$cluster_1,
                                test.use = input$test_use,
                                logfc.threshold = input$logfc_threshold,
                                min.pct = input$min_pct,
                                only.pos = isTRUE(input$only_pos),
                                verbose = FALSE
                            )
                            format_marker_results(result, cluster = input$cluster_1)
                        },
                        "compare" = {
                            req(input$cluster_1, input$cluster_2)
                            result <- FindMarkers(
                                obj,
                                ident.1 = input$cluster_1,
                                ident.2 = input$cluster_2,
                                test.use = input$test_use,
                                logfc.threshold = input$logfc_threshold,
                                min.pct = input$min_pct,
                                only.pos = isTRUE(input$only_pos),
                                verbose = FALSE
                            )
                            format_marker_results(result, cluster = input$cluster_1)
                        },
                        "conditions" = {
                            req(
                                input$cluster_1,
                                input$condition_column,
                                input$condition_value_1,
                                input$condition_value_2
                            )

                            if (identical(input$condition_value_1, input$condition_value_2)) {
                                stop("Choose two different condition values for comparison.",
                                     call. = FALSE)
                            }

                            condition_labels <- as.character(obj[[input$condition_column]][, 1])
                            if (all(is.na(condition_labels))) {
                                stop("Selected metadata column does not contain usable values.",
                                     call. = FALSE)
                            }

                            cluster_cells <- names(Seurat::Idents(obj))[
                                as.character(Seurat::Idents(obj)) == input$cluster_1
                            ]
                            if (!length(cluster_cells)) {
                                stop("Selected cluster has no cells.", call. = FALSE)
                            }

                            obj_subset <- subset(obj, cells = cluster_cells)
                            subset_labels <- as.character(obj_subset[[input$condition_column]][, 1])
                            available_values <- unique(subset_labels[!is.na(subset_labels)])
                            missing_values <- setdiff(
                                c(input$condition_value_1, input$condition_value_2),
                                available_values
                            )
                            if (length(missing_values)) {
                                stop(
                                    "The selected cluster does not contain condition value(s): ",
                                    paste(missing_values, collapse = ", "),
                                    ". Choose values present within cluster ", input$cluster_1, ".",
                                    call. = FALSE
                                )
                            }

                            obj_subset$ascseurat_condition <- subset_labels
                            Seurat::Idents(obj_subset) <- "ascseurat_condition"

                            condition_result <- FindMarkers(
                                obj_subset,
                                ident.1 = input$condition_value_1,
                                ident.2 = input$condition_value_2,
                                test.use = input$test_use,
                                logfc.threshold = input$logfc_threshold,
                                min.pct = input$min_pct,
                                only.pos = isTRUE(input$only_pos),
                                verbose = FALSE
                            )

                            condition_result$cluster <- input$cluster_1
                            condition_result$condition_column <- input$condition_column
                            condition_result$comparison <- paste(
                                input$cluster_1,
                                input$condition_value_1,
                                "vs",
                                input$condition_value_2,
                                sep = "_"
                            )
                            format_marker_results(condition_result, cluster = input$cluster_1)
                        },
                        "conserved" = {
                            req(isTRUE(allow_conserved), input$cluster_1, input$conserved_grouping)

                            if (!requireNamespace("metap", quietly = TRUE)) {
                                stop(
                                    "Package 'metap' is required for conserved markers. Install with: install.packages('metap')",
                                    call. = FALSE
                                )
                            }
                            if (!input$conserved_grouping %in% colnames(obj@meta.data)) {
                                stop(
                                    "Grouping column '", input$conserved_grouping,
                                    "' not found in metadata.",
                                    call. = FALSE
                                )
                            }
                            if (identical(input$test_use, "roc")) {
                                stop(
                                    "The 'roc' test reports AUC scores instead of p-values, so markers ",
                                    "cannot be combined across groups. Choose another statistical test ",
                                    "for conserved markers.",
                                    call. = FALSE
                                )
                            }

                            # Seurat skips groups with fewer than min.cells.group cells
                            # of the cluster; conservation needs at least two tested groups.
                            min_cells <- 3
                            group_labels <- as.character(obj[[input$conserved_grouping]][, 1])
                            cluster_group_sizes <- table(
                                group_labels[as.character(Seurat::Idents(obj)) == input$cluster_1]
                            )
                            tested_groups <- names(cluster_group_sizes)[cluster_group_sizes >= min_cells]
                            if (length(tested_groups) < 2) {
                                stop(
                                    "Cluster ", input$cluster_1, " has at least ", min_cells,
                                    " cells in ",
                                    if (length(tested_groups)) {
                                        paste0("only one group ('", tested_groups, "')")
                                    } else {
                                        "no group"
                                    },
                                    " of '", input$conserved_grouping,
                                    "', so conservation across groups cannot be assessed. ",
                                    "Use 'Find markers for one cluster' instead.",
                                    call. = FALSE
                                )
                            }

                            result <- FindConservedMarkers(
                                obj,
                                ident.1 = input$cluster_1,
                                grouping.var = input$conserved_grouping,
                                test.use = input$test_use,
                                logfc.threshold = input$logfc_threshold,
                                min.pct = input$min_pct,
                                min.cells.group = min_cells,
                                only.pos = isTRUE(input$only_pos),
                                verbose = FALSE
                            )
                            # FindConservedMarkers reports only unadjusted max_pval;
                            # filter on the largest per-group adjusted p-value instead.
                            adj_cols <- grep("_p_val_adj$", colnames(result), value = TRUE)
                            result$max_p_val_adj <- apply(
                                result[, adj_cols, drop = FALSE], 1, max
                            )
                            format_marker_results(result, cluster = input$cluster_1)
                        }
                    ),
                    error = function(e) {
                        # req() inside the branches signals an error with no message.
                        if (nzchar(conditionMessage(e))) {
                            showNotification(
                                paste0("Differential expression failed: ", conditionMessage(e)),
                                type = "error",
                                duration = 12
                            )
                        } else {
                            showNotification(
                                "Differential expression was not run: select the cluster(s) and condition values for this comparison first.",
                                type = "warning",
                                duration = 12
                            )
                        }
                        NULL
                    }
                )

                req(!is.null(result))
                setProgress(0.8, message = "Filtering results...")

                if (!"gene" %in% colnames(result)) {
                    result$gene <- rownames(result)
                }
                result <- format_marker_results(result)

                if ("p_val_adj" %in% colnames(result)) {
                    keep <- !is.na(result$p_val_adj) & result$p_val_adj < input$pval_cutoff
                    result <- result[keep, , drop = FALSE]
                } else if ("max_p_val_adj" %in% colnames(result)) {
                    keep <- !is.na(result$max_p_val_adj) & result$max_p_val_adj < input$pval_cutoff
                    result <- result[keep, , drop = FALSE]
                }

                rownames(result) <- NULL
                setProgress(1, message = paste0(nrow(result), " markers found"))
                de_result(result)
            })
        })

        markers <- reactive({
            req(de_result())
            de_result()
        })

        output$markers_table <- renderReactable({
            req(markers())
            df <- markers()
            if (!nrow(df)) {
                return(reactable(data.frame(Message = "No markers found with the current settings")))
            }
            num_cols <- vapply(df, is.numeric, logical(1))
            df[num_cols] <- lapply(df[num_cols], function(x) signif(x, 4))
            reactable(
                df,
                rownames = FALSE,
                searchable = TRUE,
                striped = TRUE,
                highlight = TRUE,
                defaultPageSize = 20
            )
        })

        output$download_markers <- downloadHandler(
            filename = function() paste0("markers_", Sys.Date(), ".csv"),
            content = function(file) {
                write.csv(markers(), file, row.names = FALSE)
            }
        )

        markers
    })
}
