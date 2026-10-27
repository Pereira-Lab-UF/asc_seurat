#' Clustering Module - UI
#'
#' @param id Module namespace ID.
#' @param step_number Numeric step number to show in the workflow title.
#' @param download_label Character label for the clustered object download.
#' @param tsne_option Logical. Show the option to also compute tSNE
#'   (single-sample workflows only).
#' @export
#' @keywords internal
mod_clustering_ui <- function(id,
                              step_number = 4,
                              download_label = "Download Clustered Object",
                              tsne_option = TRUE) {
    ns <- NS(id)
    if (is.null(step_number) || !is.finite(step_number)) {
        step_number <- 4L
    }
    step_number <- as.integer(step_number)

    step_card(
        title = paste0("Step ", step_number, ": Clustering & Visualization"),
        next_step = paste0(
            "After clustering, proceed to Step ", step_number + 1L,
            ": Differential Expression, then Step ", step_number + 2L,
            ": Gene Visualization"
        ),

        layout_columns(
            col_widths = c(6, 6),
            card(
                card_body(
                    numericInput(ns("resolution"),
                                 "Clustering resolution",
                                 value = 0.6, min = 0.05, max = 3, step = 0.1),
                    tags$p(class = "text-muted small",
                           "Uses the selected number of PCs from the previous step. Typical resolution: 0.4-1.2; higher values create more clusters."),
                    if (isTRUE(tsne_option)) {
                        checkboxInput(ns("run_tsne"),
                                      "Also compute tSNE (slow for large datasets)",
                                      value = TRUE)
                    }
                )
            ),
            card(
                card_body(
                    action_btn(ns("run_clustering"),
                               "Run Clustering",
                               icon = icon("project-diagram")),
                    tags$hr(),
                    conditionalPanel(
                        condition = sprintf("output['%s'] == 'true'", ns("clustering_ready")),
                        busy_download_button(
                            ns("download_rds"),
                            download_label
                        )
                    )
                )
            )
        ),

        conditionalPanel(
            condition = sprintf("output['%s'] == 'true'", ns("reanalyze_pending")),
            tags$hr(),
            layout_columns(
                col_widths = c(8, 4),
                card(
                    class = "normalization-elbow-card",
                    card_header("Updated Elbow Plot - Reanalyzed Subset"),
                    card_body(
                        class = "plot-container plot-container-wide",
                        with_spinner(plotOutput(ns("reanalyze_elbow"), height = "440px", width = "100%"))
                    )
                ),
                card(
                    class = "reanalyze-help-card h-100",
                    card_header("Next Step"),
                    card_body(
                        class = "d-flex flex-column justify-content-center",
                        tags$p(
                            "The selected cells were reprocessed as a new subset with a fresh PCA."
                        ),
                        tags$p(
                            class = "mb-0",
                            "Run clustering again using the selected number of PCs, capped to the subset PCA."
                        )
                    )
                )
            )
        ),

        # Dimensionality reduction plots
        conditionalPanel(
            condition = sprintf("output['%s'] == 'true'", ns("clustering_ready")),
            tags$hr(),
            uiOutput(ns("dimred_plots")),
            uiOutput(ns("sample_split_plot_card")),

            # Cluster info
            layout_columns(
                col_widths = c(8, 4),
                card(
                    card_header("Cells per Cluster"),
                    card_body(
                        with_spinner(reactableOutput(ns("cluster_table")))
                    )
                ),
                card_plot_download("umap", ns("download_umap"), ns = ns)
            ),

            uiOutput(ns("phase_plot_ui")),

            # Cluster filtering
            card(
                class = "cluster-filter-card",
                card_header("Cluster Selection / Exclusion (optional)"),
                card_body(
                    layout_columns(
                        col_widths = c(2, 3, 4, 3),
                        radioButtons(ns("filter_clusters"),
                                     "Reanalyze subset?",
                                     choices = list("No" = 0, "Yes" = 1),
                                     selected = 0),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 1", ns("filter_clusters")),
                            radioButtons(
                                ns("reanalyze_strategy"),
                                "After subsetting:",
                                choices = list(
                                    "Recompute PCA and clusters" = "recompute",
                                    "Keep current reduction and clusters" = "keep"
                                ),
                                selected = "recompute"
                            ),
                            radioButtons(ns("filter_mode"),
                                         "Mode:",
                                         choices = list("Select" = "select",
                                                        "Exclude" = "exclude"),
                                         selected = "select")
                        ),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 1", ns("filter_clusters")),
                            with_spinner(uiOutput(ns("cluster_picker")))
                        ),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 1", ns("filter_clusters")),
                            action_btn(ns("reanalyze"),
                                       "Reanalyze Subset",
                                       icon = icon("redo"))
                        )
                    )
                )
            )
        )
    )
}


#' Clustering Module - Server
#'
#' @param id Module namespace ID.
#' @param norm_result Reactive Seurat object or a list of reactives returned by
#'   the normalization module.
#' @return Reactive clustered Seurat object.
#' @export
#' @keywords internal
mod_clustering_server <- function(id, norm_result) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        plot_width <- function(output_id, fallback = 1000) {
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
        get_obj <- if (is.list(norm_result) && !is.function(norm_result)) {
            norm_result$obj
        } else {
            norm_result
        }
        get_pcs <- if (is.list(norm_result) && !is.function(norm_result)) {
            norm_result$n_pcs
        } else {
            reactive(NULL)
        }
        state_obj <- reactiveVal(NULL)
        pending_obj <- reactiveVal(NULL)
        pending_reanalysis <- reactiveVal(FALSE)
        # Subset that later Run Clustering calls start from, until the
        # upstream object changes.
        subset_base <- reactiveVal(NULL)
        min_subset_cells <- 10L

        sample_metadata_col <- function(obj) {
            candidates <- c("samples", "sample_name", "orig.ident")
            found <- candidates[candidates %in% colnames(obj@meta.data)]
            if (length(found)) found[[1]] else NULL
        }

        is_single_sample <- function(obj) {
            sample_col <- sample_metadata_col(obj)
            if (is.null(sample_col)) {
                return(TRUE)
            }

            sample_values <- unique(as.character(obj[[sample_col]][, 1]))
            sample_values <- sample_values[!is.na(sample_values) & nzchar(sample_values)]
            length(sample_values) <= 1
        }

        cluster_sort_levels <- function(values) {
            levels <- unique(as.character(values))
            numeric_levels <- suppressWarnings(as.numeric(levels))
            if (all(!is.na(numeric_levels))) {
                levels[order(numeric_levels)]
            } else {
                sort(levels)
            }
        }

        default_reduction <- function(obj) {
            preferred <- attr(obj, "ascseurat_default_reduction", exact = TRUE)
            candidates <- c(preferred, "integrated.rpca", "harmony", "pca")
            candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
            found <- candidates[candidates %in% names(obj@reductions)]
            if (length(found)) found[[1]] else "pca"
        }

        selected_n_pcs <- function(obj) {
            requested <- tryCatch(get_pcs(), error = function(e) NULL)
            requested <- requested %||%
                attr(obj, "ascseurat_n_pcs", exact = TRUE) %||%
                attr(obj, "ascseurat_clustering_n_pcs", exact = TRUE) %||%
                30L
            requested <- suppressWarnings(as.integer(requested[[1]]))
            if (!is.finite(requested) || requested < 1L) {
                requested <- 30L
            }

            reduction <- default_reduction(obj)
            available <- tryCatch(
                ncol(Seurat::Embeddings(obj, reduction)),
                error = function(e) requested
            )
            if (!is.finite(available) || available < 1L) {
                available <- requested
            }
            available <- as.integer(available)
            if (requested > available) {
                showNotification(
                    paste0(
                        "Clustering used the ", available, " dimensions available in the '",
                        reduction, "' reduction instead of the requested ", requested, "."
                    ),
                    type = "warning",
                    duration = 8
                )
            }

            max(1L, min(requested, available))
        }

        normalization_method_for_subset <- function(obj) {
            if ("SCT" %in% names(obj@assays)) {
                "sctransform"
            } else {
                "lognorm"
            }
        }

        # Repeats the normalization recorded by run_normalization(); objects
        # without a record (e.g. uploaded files) use the defaults.
        normalize_subset <- function(obj) {
            settings <- obj@misc$ascseurat_normalization
            meta_cols <- colnames(obj@meta.data)
            vars <- intersect(settings$vars_to_regress, meta_cols)
            # Regressing a variable with a single value fails for categorical
            # variables (e.g. a subset of one batch) and does nothing otherwise.
            single_value <- vapply(vars, function(v) {
                length(unique(stats::na.omit(obj@meta.data[[v]]))) < 2L
            }, logical(1))
            if (any(single_value)) {
                showNotification(
                    paste0(
                        "Not regressed in the subset because it has a single value: ",
                        paste(vars[single_value], collapse = ", "), "."
                    ),
                    type = "warning",
                    duration = 10
                )
                vars <- vars[!single_value]
            }
            regress_cc <- isTRUE(settings$regress_cc) &&
                (isTRUE(settings$score_cell_cycle) ||
                    all(c("S.Score", "G2M.Score") %in% meta_cols))
            run_normalization(
                obj,
                method = settings$method %||% normalization_method_for_subset(obj),
                nfeatures = settings$nfeatures %||% 3000,
                scale_factor = settings$scale_factor %||% 10000,
                vars_to_regress = if (length(vars)) vars,
                score_cell_cycle = isTRUE(settings$score_cell_cycle),
                cc_organism = settings$cc_organism %||% "human",
                regress_cc = regress_cc
            )
        }

        # Repeats the batch correction on a subset as far as its sample sizes
        # allow: samples too small to be their own batch are grouped with the
        # largest sample. Returns the integrated object, or a character
        # reason when batch correction cannot be repeated.
        reintegrate_subset <- function(obj, integration) {
            sample_col <- integration$sample_col
            # RPCA scores anchors on 30 neighbours within each sample; Harmony
            # and per-sample SCTransform fail on a handful of cells.
            min_batch_cells <- if (identical(integration$method, "rpca")) 31L else min_subset_cells
            batches <- as.character(obj@meta.data[[sample_col]])
            sizes <- table(batches[!is.na(batches) & nzchar(batches)])
            usable <- names(sizes)[sizes >= min_batch_cells]
            if (length(usable) < 2L) {
                return(paste0("fewer than two samples have at least ", min_batch_cells, " cells"))
            }
            small <- setdiff(names(sizes), usable)
            largest <- names(sizes)[which.max(sizes)]
            batches[!batches %in% usable] <- largest
            dims <- integration$dims
            if (identical(integration$method, "rpca")) {
                # Each sample needs more cells than RPCA dimensions.
                dims <- dims[dims < min(table(batches))]
            }
            settings <- obj@misc$ascseurat_normalization

            result <- tryCatch(
                {
                    obj$ascseurat_batch <- batches
                    # run_integration() expects v5 RNA layers split by
                    # sample; v3 assays (e.g. Seurat v4 .rds) have no layers.
                    if (inherits(obj[["RNA"]], "StdAssay")) {
                        obj[["RNA"]] <- split(
                            SeuratObject::JoinLayers(obj[["RNA"]]),
                            f = batches
                        )
                    }
                    integrated <- run_integration(
                        obj,
                        method = integration$method,
                        dims = dims,
                        sample_col = "ascseurat_batch",
                        normalization_method = integration$normalization_method,
                        nfeatures = settings$nfeatures %||% 3000,
                        scale_factor = settings$scale_factor %||% 10000
                    )
                    # RPCA on samples of a few dozen cells can return rows of
                    # NA, which clustering cannot use.
                    reduction <- attr(integrated, "ascseurat_default_reduction", exact = TRUE)
                    if (anyNA(Seurat::Embeddings(integrated, reduction))) {
                        stop("the integrated reduction has missing values", call. = FALSE)
                    }
                    integrated
                },
                error = function(e) conditionMessage(e)
            )
            if (is.character(result)) {
                return(result)
            }
            result$ascseurat_batch <- NULL
            # Later subsets repeat the original settings.
            result@misc$ascseurat_integration <- integration
            if (length(small)) {
                showNotification(
                    paste0(
                        "Samples with fewer than ", min_batch_cells, " cells in the subset (",
                        paste(small, collapse = ", "), ") were not batch-corrected separately; ",
                        "their cells were grouped with the largest sample (", largest, ")."
                    ),
                    type = "warning",
                    duration = 12
                )
            }
            result
        }

        prepare_subset_for_reclustering <- function(obj) {
            # JoinLayers() keeps counts.<s>/data.<s> next to the joined layers
            # when only one sample remains, so drop those copies.
            obj <- tryCatch(
                {
                    rna <- SeuratObject::JoinLayers(obj[["RNA"]])
                    layers <- SeuratObject::Layers(rna)
                    for (joined in intersect(c("counts", "data"), layers)) {
                        for (layer in grep(paste0("^", joined, "\\."), layers, value = TRUE)) {
                            SeuratObject::LayerData(rna, layer = layer) <- NULL
                        }
                    }
                    obj[["RNA"]] <- rna
                    obj
                },
                error = function(e) obj
            )
            obj@reductions <- list()
            obj@graphs <- list()
            if ("neighbors" %in% slotNames(obj)) {
                obj@neighbors <- list()
            }

            cols_to_drop <- grep(
                "(^seurat_clusters$|_snn$|_snn_res\\.)",
                colnames(obj@meta.data),
                value = TRUE
            )
            if (length(cols_to_drop)) {
                obj@meta.data <- obj@meta.data[
                    ,
                    setdiff(colnames(obj@meta.data), cols_to_drop),
                    drop = FALSE
                ]
            }

            sample_col <- sample_metadata_col(obj)
            if (!is.null(sample_col)) {
                Idents(obj) <- factor(as.character(obj[[sample_col]][, 1]))
            } else {
                Idents(obj) <- factor(rep("subset", ncol(obj)))
            }

            attr(obj, "ascseurat_clustering_reduction") <- NULL
            attr(obj, "ascseurat_default_reduction") <- NULL
            obj
        }

        # NULL while the upstream is req()-blocked, so the reset also runs when
        # the previous object stops being available.
        upstream_obj <- reactive(tryCatch(get_obj(), error = function(e) NULL))
        # The upstream re-runs whenever an input it reads changes, even when it
        # returns the same object; only a different object resets the results.
        last_upstream <- NULL

        observeEvent(upstream_obj(), {
            previous <- last_upstream
            last_upstream <<- upstream_obj()
            # Nothing can be stored while no object is available, so there is
            # nothing to clear when one (e.g. the initial one) appears.
            if (is.null(previous) || identical(last_upstream, previous)) {
                return(invisible(NULL))
            }
            state_obj(NULL)
            pending_obj(NULL)
            subset_base(NULL)
            pending_reanalysis(FALSE)
        }, ignoreNULL = FALSE)

        output$clustering_ready <- shiny::renderText({
            if (is.null(state_obj())) "false" else "true"
        })
        shiny::outputOptions(output, "clustering_ready", suspendWhenHidden = FALSE)

        output$reanalyze_pending <- shiny::renderText({
            if (isTRUE(pending_reanalysis())) "true" else "false"
        })
        shiny::outputOptions(output, "reanalyze_pending", suspendWhenHidden = FALSE)

        observeEvent(input$run_clustering, {
            obj_for_clustering <- subset_base()
            if (is.null(obj_for_clustering)) {
                obj_for_clustering <- get_obj()
            }
            req(obj_for_clustering)

            tryCatch(
                withProgress(message = "Clustering...", value = 0.3, {
                    reduction <- default_reduction(obj_for_clustering)
                    n_dims <- selected_n_pcs(obj_for_clustering)
                    # The checkbox only exists in single-sample workflows.
                    tsne_requested <- isTRUE(input$run_tsne)
                    single_sample <- is_single_sample(obj_for_clustering)
                    if (tsne_requested && !single_sample) {
                        showNotification(
                            "tSNE is skipped because the object has several samples; the UMAP split by sample is shown instead.",
                            type = "message",
                            duration = 8
                        )
                    }
                    obj <- run_clustering(
                        obj = obj_for_clustering,
                        dims = seq_len(n_dims),
                        resolution = input$resolution,
                        run_tsne = tsne_requested && single_sample,
                        reduction = reduction
                    )
                    attr(obj, "ascseurat_clustering_n_pcs") <- n_dims
                    setProgress(1, message = "Clustering complete!")
                    state_obj(obj)
                    pending_obj(NULL)
                    pending_reanalysis(FALSE)
                }),
                error = function(e) {
                    showNotification(
                        paste0("Clustering failed: ", conditionMessage(e)),
                        type = "error",
                        duration = 12
                    )
                }
            )
        })

        active_obj <- reactive({
            req(state_obj())
            state_obj()
        })

        observeEvent(input$reanalyze, {
            req(active_obj())
            selected <- input$selected_clusters
            if (is.null(selected) || !length(selected)) {
                showNotification(
                    "Select at least one cluster before reanalyzing the subset.",
                    type = "warning",
                    duration = 6
                )
                return(invisible(NULL))
            }

            obj <- active_obj()
            in_selection <- as.character(Idents(obj)) %in% selected
            keep_cells <- if (identical(input$filter_mode, "exclude")) {
                !in_selection
            } else {
                in_selection
            }
            # Check the size first: subset() errors when no cells match.
            n_subset <- sum(keep_cells)
            if (n_subset < min_subset_cells) {
                showNotification(
                    if (n_subset == 0L) {
                        "No cells remain after applying the requested subset operation."
                    } else {
                        paste0(
                            "The requested subset has only ", n_subset,
                            " cells; at least ", min_subset_cells,
                            " cells are needed to reanalyze it."
                        )
                    },
                    type = "error",
                    duration = 8
                )
                return(invisible(NULL))
            }

            tryCatch(
                withProgress(message = "Reanalyzing subset...", value = 0.3, {
                    obj <- subset(obj, cells = colnames(obj)[keep_cells])

                    if (identical(input$reanalyze_strategy, "keep")) {
                        setProgress(1, message = "Subset kept with current clustering.")
                        pending_obj(NULL)
                        pending_reanalysis(FALSE)
                        subset_base(obj)
                        state_obj(obj)
                        showNotification(
                            "Subset updated with the current reductions and cluster labels.",
                            type = "message",
                            duration = 6
                        )
                        return(invisible(NULL))
                    }

                    setProgress(0.6, message = "Recomputing PCA for the subset...")
                    integration <- obj@misc$ascseurat_integration
                    subset_obj <- prepare_subset_for_reclustering(obj)
                    n_samples <- 0L
                    if (!is.null(integration) &&
                        integration$sample_col %in% colnames(subset_obj@meta.data)) {
                        sample_values <- as.character(subset_obj@meta.data[[integration$sample_col]])
                        n_samples <- length(unique(
                            sample_values[!is.na(sample_values) & nzchar(sample_values)]
                        ))
                    }

                    integrated <- NULL
                    if (n_samples >= 2L) {
                        setProgress(0.6, message = paste0(
                            "Re-integrating the subset with ", integration$method, "..."
                        ))
                        integrated <- reintegrate_subset(subset_obj, integration)
                    }
                    # Without a usable batch correction the subset is
                    # reprocessed like a single sample.
                    subset_obj <- if (inherits(integrated, "Seurat")) {
                        integrated
                    } else {
                        normalize_subset(subset_obj)
                    }

                    pending_obj(subset_obj)
                    subset_base(subset_obj)
                    pending_reanalysis(TRUE)
                    state_obj(NULL)
                    setProgress(1, message = "Subset PCA ready.")
                    if (!is.null(integration) && n_samples < 2L) {
                        showNotification(
                            "Only one sample remains in the subset, so it was reprocessed without batch correction.",
                            type = "warning",
                            duration = 10
                        )
                    }
                    if (is.character(integrated)) {
                        showNotification(
                            paste0(
                                "Batch correction with ", integration$method,
                                " could not be repeated for this subset (", integrated,
                                "), so it was reprocessed without batch correction."
                            ),
                            type = "warning",
                            duration = 12
                        )
                    }
                    showNotification(
                        "Subset PCA updated. Review the new elbow plot, then run clustering again using the selected number of PCs.",
                        type = "message",
                        duration = 8
                    )
                }),
                error = function(e) {
                    showNotification(
                        paste0("Subset reanalysis failed: ", conditionMessage(e)),
                        type = "error",
                        duration = 12
                    )
                }
            )
        })

        observeEvent(input$apply_rename, {
            req(active_obj())
            obj <- active_obj()
            clusters <- levels(Idents(obj))
            new_names <- vapply(seq_along(clusters), function(i) {
                value <- input[[paste0("rename_", i)]]
                if (is.null(value) || !nzchar(trimws(value))) {
                    clusters[[i]]
                } else {
                    trimws(value)
                }
            }, character(1))
            names(new_names) <- clusters

            obj <- do.call(Seurat::RenameIdents, c(list(object = obj), as.list(new_names)))
            obj$seurat_clusters <- as.character(Idents(obj))
            state_obj(obj)
            showNotification("Cluster names updated.", type = "message", duration = 5)
        })

        output$umap_clusters <- renderPlot({
            obj <- active_obj()
            req("umap" %in% names(obj@reductions))
            Seurat::DimPlot(
                object = obj,
                reduction = "umap",
                label = TRUE
            ) +
                ggplot2::labs(color = "Cluster")
        }, height = function() 500, width = function() plot_width("umap_clusters"))

        output$dimred_plots <- renderUI({
            req(active_obj())
            umap_card <- card(
                card_header("UMAP - Clusters"),
                card_body(
                    class = "plot-container plot-container-wide",
                    with_spinner(plotOutput(ns("umap_clusters"), height = "500px", width = "100%"))
                )
            )
            # Single-sample runs without tSNE show the UMAP across the full row.
            if (is_single_sample(active_obj()) &&
                !"tsne" %in% names(active_obj()@reductions)) {
                return(umap_card)
            }
            title <- if (is_single_sample(active_obj())) {
                "tSNE - Clusters"
            } else {
                "UMAP - Samples"
            }

            layout_columns(
                col_widths = c(6, 6),
                umap_card,
                card(
                    card_header(title),
                    card_body(
                        class = "plot-container plot-container-wide",
                        with_spinner(plotOutput(ns("secondary_plot"), height = "500px", width = "100%"))
                    )
                )
            )
        })

        output$secondary_plot <- renderPlot({
            req(active_obj())
            if (is_single_sample(active_obj())) {
                req("tsne" %in% names(active_obj()@reductions))
                Seurat::DimPlot(
                    object = active_obj(),
                    reduction = "tsne",
                    label = TRUE
                ) +
                    ggplot2::labs(color = "Cluster")
            } else {
                sample_col <- sample_metadata_col(active_obj())
                req(!is.null(sample_col))
                Seurat::DimPlot(
                    object = active_obj(),
                    reduction = "umap",
                    group.by = sample_col,
                    shuffle = TRUE,
                    seed = 42,
                    pt.size = 0.35,
                    alpha = 0.6
                ) + ggplot2::labs(color = "samples")
            }
        }, height = function() 500, width = function() plot_width("secondary_plot"))

        output$sample_split_plot_card <- renderUI({
            req(active_obj())
            if (is_single_sample(active_obj())) {
                return(NULL)
            }
            tags$div(
                class = "mt-3",
                card(
                    card_header("UMAP - Clusters Split by Samples"),
                    card_body(
                        class = "plot-container plot-container-wide",
                        with_spinner(plotOutput(ns("sample_split_umap"), height = "520px", width = "100%"))
                    )
                )
            )
        })

        output$sample_split_umap <- renderPlot({
            req(active_obj())
            obj <- active_obj()
            sample_col <- sample_metadata_col(obj)
            req(!is.null(sample_col), "umap" %in% names(obj@reductions))
            group_col <- if ("seurat_clusters" %in% colnames(obj@meta.data)) {
                "seurat_clusters"
            } else {
                NULL
            }
            plot_args <- list(
                object = obj,
                reduction = "umap",
                split.by = sample_col,
                label = TRUE
            )
            if (!is.null(group_col)) {
                plot_args$group.by <- group_col
            }
            do.call(Seurat::DimPlot, plot_args) +
                ggplot2::labs(color = "Cluster")
        }, height = function() 520, width = function() plot_width("sample_split_umap"))

        # Cluster table
        output$cluster_table <- renderReactable({
            req(active_obj())
            obj <- active_obj()
            clusters <- as.character(Idents(obj))
            cluster_levels <- cluster_sort_levels(clusters)
            clusters <- factor(clusters, levels = cluster_levels)
            cluster_counts <- table(clusters)
            df <- data.frame(
                Cluster = cluster_levels,
                check.names = FALSE
            )
            sample_col <- sample_metadata_col(obj)
            if (!is.null(sample_col)) {
                sample_values <- as.character(obj[[sample_col]][, 1])
                sample_counts <- as.data.frame.matrix(table(clusters, sample_values))
                sample_counts <- sample_counts[cluster_levels, , drop = FALSE]
                sample_counts[is.na(sample_counts)] <- 0L
                sample_counts[] <- lapply(sample_counts, as.integer)
                colnames(sample_counts) <- paste0(colnames(sample_counts), " (n of cells)")
                df <- cbind(df, sample_counts)
            }
            df[["Total (n of cells)"]] <- as.integer(cluster_counts)
            rownames(df) <- NULL

            reactable(
                df,
                defaultPageSize = 20,
                striped = TRUE,
                highlight = TRUE,
                rownames = FALSE
            )
        })

        # Dynamic cluster picker
        output$cluster_picker <- renderUI({
            req(active_obj())
            clusters <- levels(Idents(active_obj()))
            checkboxGroupInput(ns("selected_clusters"),
                               "Select clusters:",
                               choices = clusters,
                               inline = TRUE)
        })

        output$rename_inputs <- renderUI({
            req(active_obj())
            clusters <- levels(Idents(active_obj()))
            tagList(lapply(seq_along(clusters), function(i) {
                textInput(
                    ns(paste0("rename_", i)),
                    label = paste("Cluster", clusters[[i]]),
                    value = clusters[[i]]
                )
            }))
        })

        # Download RDS
        output$download_rds <- downloadHandler(
            filename = function() {
                paste0("seurat_clustered_", Sys.Date(), ".rds")
            },
            content = function(file) {
                on.exit(
                    session$sendCustomMessage(
                        "ascDownloadReady",
                        list(id = ns("download_rds"))
                    ),
                    add = TRUE
                )
                save_download_rds(active_obj(), file)
            }
        )

        # Download UMAP
        output$download_umap <- plot_download_handler(
            plot_reactive = reactive({
                Seurat::DimPlot(
                    object = active_obj(),
                    reduction = "umap",
                    label = TRUE
                ) +
                    ggplot2::labs(color = "Cluster")
            }),
            filename_base = "umap_clusters",
            input = input,
            prefix = "umap"
        )

        output$phase_plot_ui <- renderUI({
            req(active_obj())
            if (!"Phase" %in% colnames(active_obj()@meta.data)) {
                return(NULL)
            }

            card(
                card_header("UMAP - Cell-cycle phase"),
                card_body(
                    class = "plot-container plot-container-wide",
                    with_spinner(plotOutput(ns("phase_umap"), height = "500px", width = "100%"))
                )
            )
        })

        output$phase_umap <- renderPlot({
            req(active_obj(), "Phase" %in% colnames(active_obj()@meta.data))
            Seurat::DimPlot(
                object = active_obj(),
                reduction = "umap",
                group.by = "Phase"
            )
        }, height = function() 500, width = function() plot_width("phase_umap"))

        output$reanalyze_elbow <- renderPlot({
            req(pending_obj())
            max_dims <- min(50, ncol(Seurat::Embeddings(pending_obj(), "pca")))
            old_par <- graphics::par(no.readonly = TRUE)
            on.exit(graphics::par(old_par), add = TRUE)
            graphics::par(mar = c(4.5, 4.5, 1.5, 1))
            ElbowPlot(pending_obj(), ndims = max_dims)
        }, height = function() 440, width = function() plot_width("reanalyze_elbow"))

        output$rename_section <- renderUI({
            req(active_obj())
            card(
                class = "rename-section-card",
                card_body(
                    tags$details(
                        class = "asc-collapse-panel",
                        tags$summary(
                            class = "asc-collapse-summary",
                            icon("tags", class = "me-2"),
                            "Rename clusters (optional)"
                        ),
                        div(
                            class = "mt-3",
                            with_spinner(uiOutput(ns("rename_inputs"))),
                            action_btn(ns("apply_rename"),
                                       "Apply New Names",
                                       icon = icon("tag"))
                        )
                    )
                )
            )
        })

        return(active_obj)
    })
}
