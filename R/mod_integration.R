#' Integration Module — UI
#'
#' @param id Module namespace ID.
#' @export
#' @keywords internal
mod_integration_ui <- function(id) {
    ns <- NS(id)

    step_card(
        title = "Step 1: Multi-Sample Integration",

        layout_columns(
            col_widths = c(4, 4, 4),
            card(
                card_header("Load Samples"),
                card_body(
                    tags$div(
                        class = "alert alert-info asc-qc-readiness-note",
                        icon("circle-info", class = "me-1"),
                        "Before integration, inspect each sample individually to define QC parameters, ",
                        "or use samples that have already undergone quality control."
                    ),
                    data_path_help(),
                    numericInput(ns("n_samples"),
                                 "Number of samples",
                                 value = 2, min = 1, max = 24, step = 1),
                    uiOutput(ns("sample_table_ui"))
                )
            ),
            card(
                card_header("QC Parameters per sample"),
                card_body(
                    uiOutput(ns("qc_params_ui"))
                )
            ),
            card(
                card_header("Integration Method"),
                card_body(
                    radioButtons(ns("int_method"),
                                 "Method:",
                                 choices = list(
                                     "RPCA (Seurat v5)" = "rpca",
                                     "Harmony" = "harmony"
                                 ),
                                 selected = "rpca"),
                    input_norm_method(ns("norm_method")),
                    numericInput(ns("n_dims"), "Number of PCs for integration",
                                 value = 30, min = 5, max = 100),
                    tags$hr(),
                    action_btn(ns("run_integration"),
                               "Load & Integrate",
                               icon = icon("layer-group"))
                )
            )
        ),

        # Post-integration summary
        conditionalPanel(
            condition = sprintf("input['%s'] > 0", ns("run_integration")),
            tags$hr(),
            card(
                card_header("Integration Summary"),
                card_body(
                    with_spinner(verbatimTextOutput(ns("int_summary"))),
                    busy_download_button(
                        ns("download_int_rds"),
                        "Download Integrated Object (not clustered).",
                        class = "btn btn-outline-secondary",
                        busy_label = "Download Integrated Object (not clustered)."
                    )
                )
            ),
            layout_columns(
                col_widths = c(7, 5),
                card(
                    card_header("Integration QC Metrics"),
                    card_body(
                        class = "plot-container",
                        with_spinner(plotOutput(ns("int_qc_violin"), height = "420px"))
                    )
                ),
                card(
                    card_header("Additional Cell Filtering"),
                    card_body(
                        radioButtons(
                            ns("post_filter_mode"),
                            "Filtering mode",
                            choices = list(
                                "No additional filtering" = "none",
                                "Use filtering thresholds" = "filter"
                            ),
                            selected = "none"
                        ),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 'none'", ns("post_filter_mode")),
                            tags$p(
                                class = "text-muted small mb-0",
                                "Clustering will use the integrated object without additional cell filtering."
                            )
                        ),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 'filter'", ns("post_filter_mode")),
                            uiOutput(ns("post_filter_ui")),
                            layout_columns(
                                col_widths = c(6, 6),
                                action_btn(ns("apply_int_filter"),
                                           "Apply Filter",
                                           icon = icon("filter")),
                                actionButton(ns("reset_int_filter"),
                                             "Reset Filter",
                                             icon = icon("rotate-left"),
                                             class = "btn btn-outline-secondary w-100")
                            )
                        ),
                        tags$hr(),
                        verbatimTextOutput(ns("post_filter_summary"))
                    )
                )
            ),
            card(
                card_header("Elbow Plot - Select PCs for Clustering"),
                card_body(
                    class = "plot-container plot-container-wide",
                    layout_columns(
                        col_widths = c(8, 4),
                        with_spinner(plotOutput(ns("int_elbow"), height = "420px", width = "100%")),
                        tags$div(
                            class = "pt-2",
                            numericInput(
                                ns("cluster_n_dims"),
                                "Number of PCs to use for clustering",
                                value = 30,
                                min = 1,
                                max = 100
                            ),
                            tags$p(
                                class = "text-muted small mb-0",
                                "Choose this after reviewing the elbow plot."
                            )
                        )
                    )
                )
            )
        )
    )
}


#' Integration Module — Server
#'
#' @param id Module namespace ID.
#' @return Reactive integrated Seurat object.
#' @export
#' @keywords internal
mod_integration_server <- function(id) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        filtered_obj <- reactiveVal(NULL)
        qc_preview_cache <- reactiveVal(list())

        qc_defaults <- list(
            mito_pattern = "",
            min_genes = 200,
            max_genes = 5000,
            max_mito = 5
        )

        config_value <- function(config, col, i, default) {
            if (!col %in% colnames(config)) {
                return(default)
            }
            value <- config[[col]][i]
            if (length(value) == 0 || is.na(value) || !nzchar(trimws(as.character(value)))) {
                return(default)
            }
            value
        }

        qc_input_value <- function(id, default, numeric = FALSE) {
            value <- input[[id]]
            if (is.null(value) || length(value) == 0 ||
                (is.character(value) && !nzchar(trimws(value)))) {
                value <- default
            }
            if (isTRUE(numeric)) {
                value <- suppressWarnings(as.numeric(value))
            } else {
                value <- trimws(as.character(value))
            }
            value
        }

        sanitize_qc_config <- function(config) {
            for (col in names(qc_defaults)) {
                if (!col %in% colnames(config)) {
                    config[[col]] <- qc_defaults[[col]]
                }
            }

            config$mito_pattern <- trimws(as.character(config$mito_pattern))
            config$mito_pattern[!nzchar(config$mito_pattern)] <- qc_defaults$mito_pattern

            for (col in c("min_genes", "max_genes", "max_mito")) {
                config[[col]] <- suppressWarnings(as.numeric(config[[col]]))
                invalid <- is.na(config[[col]]) | !is.finite(config[[col]])
                config[[col]][invalid] <- qc_defaults[[col]]
            }

            config
        }

        sample_metadata_col <- function(obj) {
            candidates <- c("samples", "sample_name", "orig.ident")
            found <- candidates[candidates %in% colnames(obj@meta.data)]
            if (length(found)) found[[1]] else NULL
        }

        ensure_samples_metadata <- function(obj) {
            sample_col <- sample_metadata_col(obj)
            if (is.null(sample_col)) {
                obj$samples <- "sample_1"
            } else {
                obj$samples <- as.character(obj[[sample_col]][, 1])
            }
            obj
        }

        mito_target_line <- function(target, max_genes = 60) {
            genes <- target$genes %||% character()
            pattern <- trimws(target$pattern %||% "")
            pattern <- if (nzchar(pattern)) pattern else "no mito pattern"
            gene_text <- if (!length(genes)) {
                "none"
            } else {
                listed <- head(genes, max_genes)
                suffix <- if (length(genes) > max_genes) {
                    paste0(", ... +", length(genes) - max_genes, " more")
                } else {
                    ""
                }
                paste0(paste(listed, collapse = ", "), suffix)
            }
            paste0(
                target$sample, " (", pattern, "): ",
                length(genes), " genes - ", gene_text
            )
        }

        filter_value <- function(id, default) {
            value <- suppressWarnings(as.numeric(input[[id]]))
            if (length(value) == 0 || is.na(value) || !is.finite(value)) {
                default
            } else {
                value
            }
        }

        integration_error_type <- function(message) {
            warning_patterns <- c(
                "does not exist",
                "needs both a sample name and a data path",
                "missing required columns",
                "must be unique",
                "must satisfy",
                "must be between"
            )
            if (any(grepl(paste(warning_patterns, collapse = "|"), message, ignore.case = TRUE))) {
                "warning"
            } else {
                "error"
            }
        }

        notify_integration_error <- function(e) {
            message <- conditionMessage(e)
            showNotification(
                paste0("Integration was not run: ", message),
                type = integration_error_type(message),
                duration = 15
            )
            NULL
        }

        path_suggestion <- function(raw_path) {
            raw_path <- trimws(as.character(raw_path))
            if (!nzchar(raw_path) || grepl("^(/|[A-Za-z]:[/\\\\])", raw_path)) {
                return("")
            }

            candidates <- unique(c(
                file.path("data", raw_path),
                file.path("RDS_files", raw_path)
            ))
            for (candidate in candidates) {
                candidate_path <- resolve_input_path(candidate)
                if (file.exists(candidate_path) || dir.exists(candidate_path)) {
                    return(paste0(" Did you mean '", candidate, "'?"))
                }
            }

            ""
        }

        validate_sample_paths <- function(config) {
            resolved_paths <- vapply(
                config$data_path,
                resolve_input_path,
                character(1)
            )
            exists <- file.exists(resolved_paths) | dir.exists(resolved_paths)
            if (all(exists)) {
                return(resolved_paths)
            }

            missing <- which(!exists)
            details <- vapply(missing, function(i) {
                paste0(
                    config$sample_name[[i]], " -> ",
                    resolved_paths[[i]],
                    path_suggestion(config$data_path[[i]])
                )
            }, character(1))
            stop(
                "Data path(s) do not exist: ",
                paste(details, collapse = "; "),
                ". In Docker, confirm that the host data folder was mounted and use a container-visible path such as data/sample_name.",
                call. = FALSE
            )
        }

        apply_qc_inputs <- function(config) {
            config <- sanitize_qc_config(config)
            for (i in seq_len(nrow(config))) {
                config$mito_pattern[[i]] <- qc_input_value(
                    paste0("qc_mito_pattern_", i),
                    config$mito_pattern[[i]]
                )
                config$min_genes[[i]] <- qc_input_value(
                    paste0("qc_min_genes_", i),
                    config$min_genes[[i]],
                    numeric = TRUE
                )
                config$max_genes[[i]] <- qc_input_value(
                    paste0("qc_max_genes_", i),
                    config$max_genes[[i]],
                    numeric = TRUE
                )
                config$max_mito[[i]] <- qc_input_value(
                    paste0("qc_max_mito_", i),
                    config$max_mito[[i]],
                    numeric = TRUE
                )
            }
            sanitize_qc_config(config)
        }

        sample_qc_defaults <- function(sample_name, raw_path, mito_pattern) {
            raw_path <- trimws(as.character(raw_path %||% ""))
            if (!nzchar(raw_path)) {
                return(NULL)
            }

            resolved_path <- resolve_input_path(raw_path)
            if (!file.exists(resolved_path) && !dir.exists(resolved_path)) {
                return(NULL)
            }

            info <- file.info(resolved_path)
            modified <- if (nrow(info) && !is.na(info$mtime[[1]])) {
                as.character(info$mtime[[1]])
            } else {
                ""
            }
            cache_key <- paste(resolved_path, modified, mito_pattern, sep = "\r")
            cache <- qc_preview_cache()
            cached <- cache[[cache_key]]
            if (!is.null(cached)) {
                return(cached)
            }

            defaults <- tryCatch({
                obj <- load_input_data(
                    resolved_path,
                    project = sample_name,
                    min.cells = 0,
                    min.features = 0
                )
                if (nzchar(trimws(mito_pattern))) {
                    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
                        obj,
                        pattern = mito_pattern
                    )
                }
                values <- qc_filter_defaults(
                    obj@meta.data,
                    defaults = c(qc_defaults, list(min_counts = 0, max_counts = NA_real_))
                )
                values[c("min_genes", "max_genes", "max_mito")]
            }, error = function(e) NULL)

            if (!is.null(defaults)) {
                cache[[cache_key]] <- defaults
                qc_preview_cache(cache)
            }
            defaults
        }

        add_qc_suggestions <- function(config) {
            if (!"data_path" %in% colnames(config)) {
                return(config)
            }
            config <- sanitize_qc_config(config)
            for (i in seq_len(nrow(config))) {
                defaults <- sample_qc_defaults(
                    sample_name = config$sample_name[[i]],
                    raw_path = config$data_path[[i]],
                    mito_pattern = config$mito_pattern[[i]]
                )
                if (is.null(defaults)) {
                    next
                }
                config$min_genes[[i]] <- defaults$min_genes
                config$max_genes[[i]] <- defaults$max_genes
                config$max_mito[[i]] <- defaults$max_mito
            }
            config
        }

        qc_controls <- function(config) {
            config <- sanitize_qc_config(config)
            tagList(lapply(seq_len(nrow(config)), function(i) {
                sample_label <- config$sample_name[[i]]
                if (is.na(sample_label) || !nzchar(trimws(as.character(sample_label)))) {
                    sample_label <- paste("Sample", i)
                }
                tags$div(
                    class = "asc-qc-sample",
                    tags$div(class = "asc-qc-sample-title", sample_label),
                    tags$div(
                        class = "asc-qc-grid",
                        textInput(
                            ns(paste0("qc_mito_pattern_", i)),
                            "Mito pattern",
                            value = config_value(config, "mito_pattern", i, qc_defaults$mito_pattern)
                        ),
                        numericInput(
                            ns(paste0("qc_min_genes_", i)),
                            "Min genes",
                            value = config_value(config, "min_genes", i, qc_defaults$min_genes),
                            min = 0
                        ),
                        numericInput(
                            ns(paste0("qc_max_genes_", i)),
                            "Max genes",
                            value = config_value(config, "max_genes", i, qc_defaults$max_genes),
                            min = 0
                        ),
                        numericInput(
                            ns(paste0("qc_max_mito_", i)),
                            "Max mito %",
                            value = config_value(config, "max_mito", i, qc_defaults$max_mito),
                            min = 0,
                            max = 100
                        )
                    )
                )
            }))
        }

        sample_config_for_qc <- reactive({
            n_samples <- input$n_samples %||% 2L
            n_samples <- suppressWarnings(as.integer(n_samples))
            if (!is.finite(n_samples)) {
                n_samples <- 2L
            }
            n_samples <- max(1L, min(24L, n_samples))

            data.frame(
                sample_name = vapply(seq_len(n_samples), function(i) {
                    sample_name <- input[[paste0("sample_name_", i)]] %||% paste0("sample_", i)
                    sample_name <- trimws(as.character(sample_name))
                    if (nzchar(sample_name)) sample_name else paste0("sample_", i)
                }, character(1)),
                data_path = vapply(seq_len(n_samples), function(i) {
                    trimws(as.character(input[[paste0("sample_path_", i)]] %||% ""))
                }, character(1)),
                stringsAsFactors = FALSE
            )
        })

        qc_preview_config <- shiny::debounce(sample_config_for_qc, millis = 1000)

        output$sample_table_ui <- renderUI({
            n_samples <- input$n_samples %||% 2L
            n_samples <- suppressWarnings(as.integer(n_samples))
            if (!is.finite(n_samples)) {
                n_samples <- 2L
            }
            n_samples <- max(1L, min(24L, n_samples))

            tagList(
                tags$div(
                    class = "asc-sample-grid asc-sample-grid-header",
                    tags$strong("Sample name"),
                    tags$strong("10X data path")
                ),
                lapply(seq_len(n_samples), function(i) {
                    tags$div(
                        class = "asc-sample-grid",
                        textInput(
                            ns(paste0("sample_name_", i)),
                            label = NULL,
                            value = paste0("sample_", i),
                            placeholder = paste0("sample_", i)
                        ),
                        textInput(
                            ns(paste0("sample_path_", i)),
                            label = NULL,
                            value = "",
                            placeholder = file.path("data", paste0("sample_", i))
                        )
                    )
                })
            )
        })

        output$qc_params_ui <- renderUI({
            config <- sample_config_for_qc()
            qc_controls(config)
        })

        observeEvent(qc_preview_config(), {
            config <- add_qc_suggestions(qc_preview_config())
            config <- sanitize_qc_config(config)
            for (i in seq_len(nrow(config))) {
                updateNumericInput(
                    session,
                    paste0("qc_min_genes_", i),
                    value = config$min_genes[[i]]
                )
                updateNumericInput(
                    session,
                    paste0("qc_max_genes_", i),
                    value = config$max_genes[[i]]
                )
                updateNumericInput(
                    session,
                    paste0("qc_max_mito_", i),
                    value = config$max_mito[[i]]
                )
            }
        }, ignoreInit = TRUE)

        manual_sample_config <- function() {
            n_samples <- input$n_samples %||% 2L
            n_samples <- suppressWarnings(as.integer(n_samples))
            if (!is.finite(n_samples)) {
                n_samples <- 2L
            }
            n_samples <- max(1L, min(24L, n_samples))

            config <- data.frame(
                sample_name = character(n_samples),
                data_path = character(n_samples),
                stringsAsFactors = FALSE
            )

            for (i in seq_len(n_samples)) {
                sample_name <- trimws(input[[paste0("sample_name_", i)]] %||% "")
                data_path <- trimws(input[[paste0("sample_path_", i)]] %||% "")
                if (!nzchar(sample_name) || !nzchar(data_path)) {
                    stop(
                        "Sample ", i,
                        " needs both a sample name and a data path.",
                        call. = FALSE
                    )
                }
                config$sample_name[[i]] <- sample_name
                config$data_path[[i]] <- data_path
            }

            apply_qc_inputs(config)
        }

        load_sample_config <- function(config) {
            required_cols <- c("sample_name", "data_path")
            missing_cols <- setdiff(required_cols, colnames(config))
            if (length(missing_cols)) {
                stop(
                    "Sample configuration is missing required columns: ",
                    paste(missing_cols, collapse = ", "),
                    call. = FALSE
                )
            }

            config$sample_name <- trimws(as.character(config$sample_name))
            config$data_path <- trimws(as.character(config$data_path))
            if (any(!nzchar(config$sample_name)) || any(!nzchar(config$data_path))) {
                stop("Every configured sample needs a sample_name and data_path.",
                     call. = FALSE)
            }

            if (anyDuplicated(config$sample_name)) {
                stop("Each sample_name must be unique.", call. = FALSE)
            }

            config <- sanitize_qc_config(config)
            resolved_paths <- validate_sample_paths(config)

            setProgress(0.2, message = "Loading individual samples...")
            mito_targets <- vector("list", nrow(config))
            obj_list <- lapply(seq_len(nrow(config)), function(i) {
                sample_name <- config$sample_name[i]
                data_path <- resolved_paths[[i]]
                min_genes <- config$min_genes[[i]]
                max_genes <- config$max_genes[[i]]
                max_mito <- config$max_mito[[i]]
                mito <- config$mito_pattern[[i]]

                if (min_genes < 0 || max_genes < 0 || max_genes <= min_genes) {
                    stop(
                        "QC gene thresholds for sample '", sample_name,
                        "' must satisfy 0 <= min_genes < max_genes.",
                        call. = FALSE
                    )
                }
                if (max_mito < 0 || max_mito > 100) {
                    stop(
                        "Max mitochondrial percentage for sample '", sample_name,
                        "' must be between 0 and 100.",
                        call. = FALSE
                    )
                }

                message("Loading: ", sample_name)

                data <- load_input_data(
                    path = data_path,
                    project = sample_name,
                    min.cells = 3,
                    min.features = min_genes
                )

                mito_genes <- if (nzchar(trimws(mito))) {
                    grep(mito, rownames(data), value = TRUE)
                } else {
                    character()
                }
                mito_targets[[i]] <<- list(
                    sample = sample_name,
                    pattern = mito,
                    genes = mito_genes
                )
                data[["percent.mt"]] <- if (length(mito_genes)) {
                    PercentageFeatureSet(data, pattern = mito)
                } else {
                    0
                }

                data <- subset(
                    data,
                    subset = nFeature_RNA > min_genes &
                        nFeature_RNA < max_genes &
                        percent.mt < max_mito
                )

                extra_cols <- setdiff(
                    colnames(config),
                    c("sample_name", "data_path", names(qc_defaults))
                )
                if (length(extra_cols)) {
                    for (col in extra_cols) {
                        data[[col]] <- config[[col]][i]
                    }
                }
                data$sample_name <- sample_name
                data$samples <- sample_name

                data
            })

            setProgress(0.4, message = "Merging samples...")
            obj <- if (length(obj_list) == 1L) {
                obj_list[[1]]
            } else {
                merge(
                    obj_list[[1]],
                    y = obj_list[-1],
                    add.cell.ids = config$sample_name
                )
            }
            attr(obj, "ascseurat_mito_targets") <- mito_targets
            ensure_samples_metadata(obj)
        }

        integrated_obj <- eventReactive(input$run_integration, {
            tryCatch(
                withProgress(message = "Loading and integrating samples...", value = 0.1, {
                    config <- manual_sample_config()
                    config <- apply_qc_inputs(config)
                    obj <- load_sample_config(config)
                    obj <- ensure_samples_metadata(obj)
                    mito_targets <- attr(obj, "ascseurat_mito_targets", exact = TRUE)

                    # Integrate
                    setProgress(0.7, message = paste0("Integrating with ", input$int_method, "..."))
                    obj <- run_integration(
                        obj,
                        method = input$int_method,
                        dims = 1:input$n_dims,
                        sample_col = "samples",
                        normalization_method = input$norm_method
                    )
                    obj <- ensure_samples_metadata(obj)
                    attr(obj, "ascseurat_mito_targets") <- mito_targets
                    attr(obj, "ascseurat_n_pcs") <- input$n_dims

                    setProgress(1, message = "Integration complete!")
                    obj
                }),
                error = notify_integration_error
            )
        })

        observeEvent(integrated_obj(), {
            filtered_obj(NULL)
            updateRadioButtons(session, "post_filter_mode", selected = "none")
            obj <- integrated_obj()
            available <- tryCatch(
                ncol(Seurat::Embeddings(obj, "pca")),
                error = function(e) input$n_dims %||% 30L
            )
            available <- suppressWarnings(as.integer(available))
            if (!is.finite(available) || available < 1L) {
                available <- 30L
            }
            selected <- suppressWarnings(as.integer(input$n_dims %||% 30L))
            if (!is.finite(selected) || selected < 1L) {
                selected <- min(30L, available)
            }
            updateNumericInput(
                session,
                "cluster_n_dims",
                value = max(1L, min(selected, available)),
                max = available
            )
        }, ignoreInit = TRUE)

        active_integrated_obj <- reactive({
            req(integrated_obj())
            if (identical(input$post_filter_mode, "filter")) {
                filtered_obj() %||% integrated_obj()
            } else {
                integrated_obj()
            }
        })

        output$integration_ready <- shiny::renderText({
            obj <- tryCatch(integrated_obj(), error = function(e) NULL)
            if (is.null(obj)) "false" else "true"
        })
        shiny::outputOptions(output, "integration_ready", suspendWhenHidden = FALSE)

        selected_cluster_n_dims <- reactive({
            value <- input$cluster_n_dims %||% input$n_dims %||% 30L
            value <- suppressWarnings(as.integer(value[[1]]))
            if (!is.finite(value) || value < 1L) {
                value <- 30L
            }
            value
        })

        output$int_summary <- renderPrint({
            req(active_integrated_obj())
            obj <- active_integrated_obj()
            sample_col <- sample_metadata_col(obj)
            cat("Total cells:", ncol(obj), "\n")
            cat("Total features:", nrow(obj), "\n")
            if (!is.null(sample_col)) {
                cat("Samples:", length(unique(obj[[sample_col]][, 1])), "\n")
            }
            cat("Method:", input$int_method, "\n")
            filter_summary <- attr(obj, "ascseurat_filter_summary", exact = TRUE)
            if (!is.null(filter_summary)) {
                cat(
                    "Post-integration filter:",
                    filter_summary$n_after, "of", filter_summary$n_before,
                    "cells retained\n"
                )
            }
            mito_targets <- attr(obj, "ascseurat_mito_targets", exact = TRUE)
            if (length(mito_targets)) {
                cat("\nMitochondrial tag targets:\n")
                for (target in mito_targets) {
                    cat("- ", mito_target_line(target), "\n", sep = "")
                }
            }
        })

        output$post_filter_ui <- renderUI({
            req(integrated_obj())
            obj <- integrated_obj()
            meta <- obj@meta.data
            defaults <- qc_filter_defaults(meta)
            controls <- list()

            if ("nFeature_RNA" %in% colnames(meta)) {
                controls <- c(controls, list(
                    numericInput(
                        ns("filter_min_features"),
                        "Min genes",
                        value = defaults$min_genes,
                        min = 0
                    ),
                    numericInput(
                        ns("filter_max_features"),
                        "Max genes",
                        value = defaults$max_genes,
                        min = 0
                    )
                ))
            }
            if ("nCount_RNA" %in% colnames(meta)) {
                controls <- c(controls, list(
                    numericInput(
                        ns("filter_min_counts"),
                        "Min UMIs",
                        value = defaults$min_counts,
                        min = 0
                    ),
                    numericInput(
                        ns("filter_max_counts"),
                        "Max UMIs",
                        value = defaults$max_counts,
                        min = 0
                    )
                ))
            }
            if ("percent.mt" %in% colnames(meta)) {
                controls <- c(controls, list(
                    numericInput(
                        ns("filter_max_mito"),
                        "Max mito %",
                        value = defaults$max_mito,
                        min = 0,
                        max = 100
                    )
                ))
            }

            if (!length(controls)) {
                return(tags$p(class = "text-muted small mb-0",
                              "No QC metadata columns are available for filtering."))
            }
            do.call(tagList, controls)
        })

        observeEvent(input$apply_int_filter, {
            req(identical(input$post_filter_mode, "filter"))
            req(integrated_obj())
            obj <- integrated_obj()
            meta <- obj@meta.data
            keep <- rep(TRUE, ncol(obj))

            if ("nFeature_RNA" %in% colnames(meta)) {
                keep <- keep &
                    meta$nFeature_RNA >= filter_value("filter_min_features", -Inf) &
                    meta$nFeature_RNA <= filter_value("filter_max_features", Inf)
            }
            if ("nCount_RNA" %in% colnames(meta)) {
                keep <- keep &
                    meta$nCount_RNA >= filter_value("filter_min_counts", -Inf) &
                    meta$nCount_RNA <= filter_value("filter_max_counts", Inf)
            }
            if ("percent.mt" %in% colnames(meta)) {
                keep <- keep & meta$percent.mt <= filter_value("filter_max_mito", Inf)
            }

            cells <- colnames(obj)[keep & !is.na(keep)]
            if (!length(cells)) {
                showNotification(
                    "No cells remain after applying the integration filter.",
                    type = "error",
                    duration = 8
                )
                return(invisible(NULL))
            }

            filtered <- subset(obj, cells = cells)
            attr(filtered, "ascseurat_mito_targets") <-
                attr(obj, "ascseurat_mito_targets", exact = TRUE)
            attr(filtered, "ascseurat_default_reduction") <-
                attr(obj, "ascseurat_default_reduction", exact = TRUE)
            attr(filtered, "ascseurat_n_pcs") <-
                attr(obj, "ascseurat_n_pcs", exact = TRUE)
            attr(filtered, "ascseurat_filter_summary") <- list(
                n_before = ncol(obj),
                n_after = ncol(filtered)
            )
            filtered_obj(filtered)
            showNotification(
                paste0("Integration filter retained ", ncol(filtered), " of ", ncol(obj), " cells."),
                type = "message",
                duration = 6
            )
        })

        observeEvent(input$reset_int_filter, {
            req(integrated_obj())
            filtered_obj(NULL)
            showNotification("Integration filter reset.", type = "message", duration = 5)
        })

        observeEvent(input$post_filter_mode, {
            if (identical(input$post_filter_mode, "none")) {
                filtered_obj(NULL)
            }
        }, ignoreInit = TRUE)

        output$post_filter_summary <- renderPrint({
            req(active_integrated_obj())
            obj <- active_integrated_obj()
            filter_summary <- attr(obj, "ascseurat_filter_summary", exact = TRUE)
            if (is.null(filter_summary)) {
                cat("Current object:", ncol(obj), "cells\n")
            } else {
                cat("Current object:", filter_summary$n_after, "cells\n")
                cat("Original integrated object:", filter_summary$n_before, "cells\n")
            }
        })

        output$int_qc_violin <- renderPlot({
            req(active_integrated_obj())
            obj <- active_integrated_obj()
            meta <- obj@meta.data
            features <- intersect(
                c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                colnames(meta)
            )
            req(length(features) > 0)
            sample_col <- sample_metadata_col(obj)
            sample_values <- if (!is.null(sample_col) && sample_col %in% colnames(meta)) {
                as.character(meta[[sample_col]])
            } else {
                rep("all cells", nrow(meta))
            }

            plot_data <- do.call(rbind, lapply(features, function(feature) {
                data.frame(
                    sample = sample_values,
                    metric = feature,
                    value = suppressWarnings(as.numeric(meta[[feature]])),
                    stringsAsFactors = FALSE
                )
            }))
            plot_data <- plot_data[is.finite(plot_data$value), , drop = FALSE]
            req(nrow(plot_data) > 0)

            ggplot2::ggplot(plot_data, ggplot2::aes(x = sample, y = value, fill = sample)) +
                ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
                ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.6) +
                ggplot2::facet_wrap(~metric, scales = "free_y") +
                ggplot2::labs(x = NULL, y = NULL) +
                ggplot2::theme_minimal(base_size = 12) +
                ggplot2::theme(
                    legend.position = "none",
                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
                )
        })

        output$int_elbow <- renderPlot({
            req(active_integrated_obj())
            obj <- active_integrated_obj()
            req("pca" %in% names(obj@reductions))
            max_dims <- min(50, ncol(Seurat::Embeddings(obj, "pca")))
            old_par <- graphics::par(no.readonly = TRUE)
            on.exit(graphics::par(old_par), add = TRUE)
            graphics::par(mar = c(4.5, 4.5, 1.5, 1))
            ElbowPlot(obj, ndims = max_dims)
        })

        output$download_int_rds <- downloadHandler(
            filename = function() paste0("integrated_", Sys.Date(), ".rds"),
            content = function(file) {
                on.exit(
                    session$sendCustomMessage(
                        "ascDownloadReady",
                        list(id = ns("download_int_rds"))
                    ),
                    add = TRUE
                )
                saveRDS(active_integrated_obj(), file)
            }
        )

        return(list(
            obj = active_integrated_obj,
            n_pcs = selected_cluster_n_dims
        ))
    })
}
