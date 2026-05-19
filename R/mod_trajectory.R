#' Trajectory Inference Module - UI
#'
#' @param id Module namespace ID.
#' @export
#' @keywords internal
mod_trajectory_ui <- function(id) {
    ns <- NS(id)

    tagList(
        step_card(
            title = "Trajectory Inference",

            layout_columns(
                col_widths = c(4, 4, 4),
                card(
                    card_header("Data"),
                    card_body(
                        radioButtons(ns("rds_source"),
                                     "RDS source:",
                                     choices = list(
                                         "Upload from browser" = "upload",
                                         "Path on this computer/server" = "path"
                                     ),
                                     selected = "upload"),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 'path'", ns("rds_source")),
                            textInput(ns("rds_path"),
                                      "Path to processed Seurat RDS:",
                                      value = file.path("RDS_files", "seurat_clustered.rds"),
                                      placeholder = "RDS_files/seurat_clustered.rds")
                        ),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 'upload'", ns("rds_source")),
                            fileInput(ns("rds_file"),
                                      "Upload processed Seurat RDS:",
                                      accept = ".rds")
                        ),
                        tags$p(class = "text-muted small",
                               "Use a path for large files. Relative paths start from ",
                               tags$code(ascseurat_workdir()),
                               ". In Docker, launch with a bind mount such as ",
                               tags$code('-v "$HOME:$HOME:ro"'),
                               " before entering host paths.")
                    )
                ),
                card(
                    card_header("Method"),
                    card_body(
                        radioButtons(
                            ns("trajectory_method"),
                            "Trajectory method:",
                            choices = list(
                                "Slingshot" = "slingshot",
                                "PAGA" = "paga",
                                "Monocle 3" = "monocle3"
                            ),
                            selected = "slingshot"
                        ),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 'slingshot'", ns("trajectory_method")),
                            tags$p(
                                class = "text-muted small mb-0",
                                "Fits lineages on the clustered UMAP embedding. ",
                                "Start and end clusters can guide lineage orientation."
                            )
                        ),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 'paga'", ns("trajectory_method")),
                            tags$p(
                                class = "text-muted small mb-0",
                                "Estimates cluster connectivity and pseudotime from the ",
                                "selected root/start cluster."
                            )
                        ),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 'monocle3'", ns("trajectory_method")),
                            tags$p(
                                class = "text-muted small mb-0",
                                "Learns a Monocle 3 principal graph on the existing UMAP ",
                                "and orders cells from the selected root cluster."
                            )
                        )
                    )
                ),
                card(
                    card_header("Options"),
                    card_body(
                        radioButtons(ns("set_start"),
                                     "Set root/start cluster?",
                                     choices = list("Yes (strongly recommended)" = 1, "No" = 0),
                                     selected = 1),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 1", ns("set_start")),
                            with_spinner(uiOutput(ns("start_cluster_ui"))),
                            conditionalPanel(
                                condition = sprintf("input['%s'] == 'slingshot'", ns("trajectory_method")),
                                with_spinner(uiOutput(ns("end_cluster_ui")))
                            ),
                            conditionalPanel(
                                condition = sprintf("input['%s'] != 'slingshot'", ns("trajectory_method")),
                                tags$p(
                                    class = "text-muted small",
                                    "For PAGA, the selected root/start cluster chooses the ",
                                    "pseudotime root; for Monocle 3, it orders cells from the root. ",
                                    "End clusters are used only by Slingshot."
                                )
                            )
                        ),
                        tags$hr(),
                        action_btn(ns("run_ti"),
                                   "Run Trajectory Inference",
                                   icon = icon("timeline"))
                    )
                )
            )
        ),

        conditionalPanel(
            condition = sprintf("input['%s'] > 0", ns("run_ti")),
            step_card(
                title = "Trajectory Results",

                layout_columns(
                    col_widths = c(4, 4, 4),
                    card(
                        card_header("Pseudotime"),
                        card_body(
                            class = "plot-container",
                            with_spinner(plotOutput(ns("pseudotime_plot"), height = "450px"))
                        )
                    ),
                    card(
                        card_header("Trajectory Graph"),
                        card_body(
                            class = "plot-container",
                            with_spinner(plotOutput(ns("trajectory_plot"), height = "450px"))
                        )
                    ),
                    card(
                        card_header("Trajectory Info"),
                        card_body(
                            tags$div(
                                class = "asc-pre-wrap",
                                with_spinner(verbatimTextOutput(ns("lineage_info")))
                            ),
                            card_plot_download(
                                "ti",
                                ns("download_ti_plot"),
                                ns = ns,
                                plot_choices = list(
                                    "Pseudotime plot" = "pseudotime",
                                    "Trajectory graph" = "trajectory"
                                ),
                                plot_selected = "pseudotime"
                            )
                        )
                    )
                ),

                layout_columns(
                    col_widths = c(6, 6),
                    downloadButton(ns("download_ti_rds"),
                                   "Download Trajectory Object",
                                   class = "btn btn-outline-secondary"),
                    tags$p(class = "text-muted",
                           "Saves the Seurat object with trajectory metadata.")
                )
            )
        ),

        step_card(
            title = "Custom Trajectory Heatmap",

            layout_columns(
                col_widths = c(4, 4, 4),
                card(
                    card_header("Trajectory"),
                    card_body(
                        radioButtons(
                            ns("custom_hm_source"),
                            "Trajectory source:",
                            choices = list(
                                "Use current trajectory" = "current",
                                "Load trajectory RDS" = "rds"
                            ),
                            selected = "current"
                        ),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 'current'", ns("custom_hm_source")),
                            tags$p(
                                class = "text-muted small mb-0",
                                "Uses the trajectory generated above in this session."
                            )
                        ),
                        conditionalPanel(
                            condition = sprintf("input['%s'] == 'rds'", ns("custom_hm_source")),
                            radioButtons(
                                ns("custom_hm_rds_source"),
                                "RDS source:",
                                choices = list(
                                    "Upload from browser" = "upload",
                                    "Path on this computer/server" = "path"
                                ),
                                selected = "upload"
                            ),
                            conditionalPanel(
                                condition = sprintf("input['%s'] == 'path'", ns("custom_hm_rds_source")),
                                textInput(
                                    ns("custom_hm_rds_path"),
                                    "Path to trajectory RDS:",
                                    value = file.path("RDS_files", "trajectory.rds"),
                                    placeholder = "RDS_files/trajectory.rds"
                                )
                            ),
                            conditionalPanel(
                                condition = sprintf("input['%s'] == 'upload'", ns("custom_hm_rds_source")),
                                fileInput(
                                    ns("custom_hm_rds_file"),
                                    "Upload trajectory RDS:",
                                    accept = ".rds"
                                ),
                                preserve_file_input_scroll_js(ns("custom_hm_rds_file"))
                            ),
                            tags$p(
                                class = "text-muted small mb-0",
                                "The RDS must include trajectory_path and pseudotime metadata."
                            )
                        )
                    )
                ),
                card(
                    card_header("Genes"),
                    card_body(
                        fileInput(
                            ns("custom_hm_genes_file"),
                            "Upload gene list CSV/TSV:",
                            accept = c(
                                "text/csv", ".csv", "text/tsv", ".tsv", ".txt",
                                "text/comma-separated-values",
                                "text/tab-separated-values"
                            )
                        ),
                        preserve_file_input_scroll_js(ns("custom_hm_genes_file")),
                        radioButtons(
                            ns("custom_hm_gene_header"),
                            "File has header?",
                            choices = list("Yes" = 1, "No" = 0),
                            selected = 1,
                            inline = TRUE
                        ),
                        numericInput(
                            ns("custom_hm_max_genes"),
                            "Maximum genes to plot",
                            value = 20,
                            min = 1,
                            max = 500,
                            step = 1
                        ),
                        tags$p(
                            class = "text-muted small mb-0",
                            "The first column should contain gene names."
                        )
                    )
                ),
                card(
                    card_header("Run"),
                    card_body(
                        action_btn(
                            ns("custom_hm_run"),
                            "Generate Heatmap",
                            icon = icon("table-cells")
                        )
                    )
                )
            ),

            conditionalPanel(
                condition = sprintf("input['%s'] > 0", ns("custom_hm_run")),
                tags$hr(),
                layout_columns(
                    col_widths = c(8, 4),
                    card(
                        card_header("Trajectory Heatmap"),
                        card_body(
                            class = "plot-container",
                            tags$p(
                                class = "text-muted small mb-2",
                                "Cells are grouped by trajectory path and ordered by pseudotime within each path."
                            ),
                            with_spinner(plotOutput(ns("custom_traj_heatmap"), height = "600px"))
                        )
                    ),
                    card(
                        card_header("Heatmap Details"),
                        card_body(
                            with_spinner(verbatimTextOutput(ns("custom_traj_heatmap_info"))),
                            card_plot_download(
                                "custom_traj_hm",
                                ns("download_custom_traj_heatmap"),
                                ns = ns,
                                height_default = 10,
                                width_default = 15
                            )
                        )
                    )
                )
            )
        ),

        conditionalPanel(
            condition = sprintf(
                "input['%s'] > 0 && input['%s'] == 'slingshot'",
                ns("run_ti"), ns("trajectory_method")
            ),
            step_card(
                title = "Trajectory Differential Expression (Slingshot only)",

                layout_columns(
                    col_widths = c(6, 6),
                    card(
                        card_body(
                            selectInput(
                                ns("traj_de_method"),
                                "DE engine:",
                                choices = list(
                                    "PseudotimeDE-fast (recommended)" = "pseudotimede",
                                    "scMaSigPro" = "scmasigpro",
                                    "tradeSeq (long runtime)" = "tradeseq"
                                ),
                                selected = "pseudotimede"
                            ),
                            conditionalPanel(
                                condition = sprintf("input['%s'] == 'scmasigpro'", ns("traj_de_method")),
                                tags$p(
                                    class = "text-muted small",
                                    "scMaSigPro identifies genes that change ",
                                    "significantly along Slingshot pseudotime."
                                )
                            ),
                            conditionalPanel(
                                condition = sprintf("input['%s'] == 'tradeseq'", ns("traj_de_method")),
                                tags$p(
                                    class = "text-muted small",
                                    "tradeSeq fits spline models along Slingshot lineages. ",
                                    "It can take a long time on large datasets."
                                )
                            ),
                            numericInput(ns("n_top_genes"),
                                         "Top genes to display/plot",
                                         value = 10, min = 10, max = 500),
                            tags$p(
                                class = "text-muted small",
                                "All variable genes are tested. This setting only ",
                                "limits how many top genes are shown in the table ",
                                "and trajectory plots."
                            ),
                            conditionalPanel(
                                condition = sprintf("input['%s'] == 'tradeseq'", ns("traj_de_method")),
                                numericInput(ns("tradeseq_nknots"),
                                             "tradeSeq knots",
                                             value = 6, min = 3, max = 10)
                            )
                        )
                    ),
                    card(
                        card_body(
                            action_btn(ns("run_traj_de"),
                                       "Find Trajectory DE Genes",
                                       icon = icon("dna")),
                            conditionalPanel(
                                condition = sprintf("input['%s'] > 0", ns("run_traj_de")),
                                downloadButton(ns("download_traj_genes"),
                                               "Download Gene List (CSV)",
                                               class = "btn btn-outline-secondary w-100 mt-2")
                            )
                        )
                    )
                ),

                conditionalPanel(
                    condition = sprintf("input['%s'] > 0", ns("run_traj_de")),
                    tags$hr(),
                    layout_columns(
                        col_widths = c(6, 6),
                        card(
                            card_header("Trajectory Heatmap"),
                            card_body(
                                class = "plot-container",
                                tags$p(
                                    class = "text-muted small mb-2",
                                    "Path1, Path2, and so on correspond to Slingshot Lineage 1, Lineage 2, and so on. ",
                                    "Cells are grouped by their assigned lineage path and ordered by pseudotime within each path. ",
                                    "Genes are ranked globally by the selected DE method, so one path can dominate the heatmap ",
                                    "when it has the strongest trajectory-associated signal."
                                ),
                                with_spinner(plotOutput(ns("traj_heatmap"), height = "600px"))
                            )
                        ),
                        card(
                            card_header("Top Trajectory Genes"),
                            card_body(
                                with_spinner(reactableOutput(ns("traj_genes_table")))
                            )
                        )
                    ),
                    card(
                        card_header("Top Gene Feature Plots"),
                        card_body(
                            class = "plot-container",
                            with_spinner(plotOutput(ns("traj_feature_plot"), height = "620px"))
                        )
                    )
                )
            )
        ),

        conditionalPanel(
            condition = sprintf(
                "input['%s'] > 0 && input['%s'] == 'paga'",
                ns("run_ti"), ns("trajectory_method")
            ),
            step_card(
                title = "Connected Cluster Gene Discovery",

                layout_columns(
                    col_widths = c(4, 4, 4),
                    card(
                        card_body(
                            radioButtons(
                                ns("paga_gene_mode"),
                                "Gene discovery mode:",
                                choices = list(
                                    "Connected cluster markers" = "edge_markers",
                                    "Pseudotime-associated genes" = "pseudotime"
                                ),
                                selected = "edge_markers"
                            ),
                            tags$div(
                                class = "text-muted small mb-3",
                                tags$p(
                                    class = "mb-1",
                                    tags$strong("Connected cluster markers: "),
                                    "compares the two clusters joined by a PAGA edge to find genes enriched on either side of that connection."
                                ),
                                tags$p(
                                    class = "mb-0",
                                    tags$strong("Pseudotime-associated genes: "),
                                    "uses cells from the connected clusters and ranks genes whose expression changes with PAGA pseudotime."
                                )
                            ),
                            radioButtons(
                                ns("paga_edge_scope"),
                                "PAGA edges:",
                                choices = list(
                                    "Selected edge" = "selected",
                                    "Top connected edges" = "top"
                                ),
                                selected = "selected"
                            ),
                            tags$p(
                                class = "text-muted small",
                                "Selected edge focuses on one connection. Top connected edges runs the same test across the strongest PAGA connections."
                            ),
                            with_spinner(uiOutput(ns("paga_edge_ui"))),
                            conditionalPanel(
                                condition = sprintf("input['%s'] == 'top'", ns("paga_edge_scope")),
                                numericInput(ns("paga_top_edges"),
                                             "Number of top edges",
                                             value = 5, min = 1, max = 25)
                            )
                        )
                    ),
                    card(
                        card_body(
                            radioButtons(
                                ns("paga_genes_to_test"),
                                "Genes to test:",
                                choices = list(
                                    "Variable genes (faster)" = "variable",
                                    "All genes" = "all"
                                ),
                                selected = "variable"
                            ),
                            numericInput(ns("paga_top_genes"),
                                         "Top genes to show/plot",
                                         value = 10, min = 5, max = 500),
                            tags$p(
                                class = "text-muted small",
                                "Variable genes are faster and usually enough for exploration. ",
                                "All genes is more complete but slower. ",
                                "Top genes controls the table and plots; the CSV includes the full result."
                            )
                        )
                    ),
                    card(
                        card_body(
                            action_btn(ns("run_paga_genes"),
                                       "Find PAGA Genes",
                                       icon = icon("dna")),
                            conditionalPanel(
                                condition = sprintf("input['%s'] > 0", ns("run_paga_genes")),
                                downloadButton(ns("download_paga_genes"),
                                               "Download Gene List (CSV)",
                                               class = "btn btn-outline-secondary w-100 mt-2")
                            )
                        )
                    )
                ),

                conditionalPanel(
                    condition = sprintf("input['%s'] > 0", ns("run_paga_genes")),
                    tags$hr(),
                    layout_columns(
                        col_widths = c(7, 5),
                        card(
                            card_header("PAGA Gene Results"),
                            card_body(
                                with_spinner(reactableOutput(ns("paga_genes_table")))
                            )
                        ),
                        card(
                            card_header("Connected Cluster Dot Plot"),
                            card_body(
                                class = "plot-container",
                                with_spinner(plotOutput(ns("paga_dot_plot"), height = "520px"))
                            )
                        )
                    ),
                    card(
                        card_header("PAGA Feature Plots"),
                        card_body(
                            class = "plot-container",
                            with_spinner(plotOutput(ns("paga_feature_plot"), height = "620px"))
                        )
                    )
                )
            )
        ),

        conditionalPanel(
            condition = sprintf(
                "input['%s'] > 0 && input['%s'] == 'monocle3'",
                ns("run_ti"), ns("trajectory_method")
            ),
            step_card(
                title = "Trajectory-Variable Genes",

                layout_columns(
                    col_widths = c(4, 4, 4),
                    card(
                        card_body(
                            radioButtons(
                                ns("monocle_genes_to_test"),
                                "Genes to test:",
                                choices = list(
                                    "Variable genes (faster)" = "variable",
                                    "All genes" = "all"
                                ),
                                selected = "variable"
                            ),
                            numericInput(ns("monocle_top_genes"),
                                         "Top genes to show/plot",
                                         value = 50, min = 5, max = 500),
                            tags$p(
                                class = "text-muted small",
                                "Ranks genes that vary over the learned principal graph."
                            )
                        )
                    ),
                    card(
                        card_body(
                            action_btn(ns("run_monocle_genes"),
                                       "Find Trajectory Genes",
                                       icon = icon("dna")),
                            conditionalPanel(
                                condition = sprintf("input['%s'] > 0", ns("run_monocle_genes")),
                                downloadButton(ns("download_monocle_genes"),
                                               "Download Gene List (CSV)",
                                               class = "btn btn-outline-secondary w-100 mt-2")
                            )
                        )
                    ),
                    card(
                        card_body(
                            action_btn(ns("run_monocle_modules"),
                                       "Find Gene Modules",
                                       icon = icon("diagram-project")),
                            tags$p(
                                class = "text-muted small",
                                "Uses significant trajectory genes from the current result."
                            ),
                            conditionalPanel(
                                condition = sprintf("input['%s'] > 0", ns("run_monocle_modules")),
                                downloadButton(ns("download_monocle_modules"),
                                               "Download Module Genes (CSV)",
                                               class = "btn btn-outline-secondary w-100 mt-2")
                            )
                        )
                    )
                ),

                conditionalPanel(
                    condition = sprintf("input['%s'] > 0", ns("run_monocle_genes")),
                    tags$hr(),
                    card(
                        card_header("Trajectory-Variable Genes"),
                        card_body(
                            with_spinner(reactableOutput(ns("monocle_genes_table")))
                        )
                    ),
                    layout_columns(
                        col_widths = c(6, 6),
                        card(
                            card_header("Trajectory Gene Heatmap"),
                            card_body(
                                class = "plot-container",
                                tags$p(
                                    class = "text-muted small mb-2",
                                    "Cells are ordered by Monocle 3 pseudotime and grouped by trajectory path."
                                ),
                                with_spinner(plotOutput(ns("monocle_gene_heatmap"), height = "620px"))
                            )
                        ),
                        card(
                            card_header("Top Gene Feature Plots"),
                            card_body(
                                class = "plot-container",
                                with_spinner(plotOutput(ns("monocle_feature_plot"), height = "620px"))
                            )
                        )
                    ),
                    conditionalPanel(
                        condition = sprintf("input['%s'] > 0", ns("run_monocle_modules")),
                        layout_columns(
                            col_widths = c(5, 7),
                            card(
                                card_header("Gene Module Summary"),
                                card_body(
                                    with_spinner(reactableOutput(ns("monocle_modules_table")))
                                )
                            ),
                            card(
                                card_header("Gene Module Heatmap"),
                                card_body(
                                    class = "plot-container",
                                    with_spinner(plotOutput(ns("monocle_modules_heatmap"),
                                                            height = "520px"))
                                )
                            )
                        )
                    )
                )
            )
        )
    )
}


sort_cluster_labels <- function(labels) {
    labels <- as.character(labels)
    numeric_labels <- suppressWarnings(as.numeric(labels))
    if (length(labels) && all(!is.na(numeric_labels))) {
        labels[order(numeric_labels)]
    } else {
        sort(labels)
    }
}


cluster_levels <- function(obj) {
    levels <- levels(Idents(obj))
    if (!length(levels)) {
        levels <- unique(as.character(Idents(obj)))
    }
    sort_cluster_labels(levels)
}


seurat_clusters <- function(obj) {
    clusters <- as.character(Idents(obj))
    cluster_names <- names(Idents(obj))
    if (is.null(cluster_names) || !length(cluster_names)) {
        cluster_names <- colnames(obj)
    }
    names(clusters) <- cluster_names
    clusters
}


selected_cluster_or_null <- function(value) {
    if (!is.null(value) && length(value) && nzchar(value[[1]])) {
        as.character(value[[1]])
    } else {
        NULL
    }
}


trajectory_rds_available <- function(rds_source = "upload", rds_path = "",
                                     rds_file = NULL) {
    rds_source <- rds_source %||% "upload"
    if (identical(rds_source, "upload")) {
        datapath <- rds_file$datapath %||% ""
        return(length(datapath) > 0 && nzchar(datapath[[1]]) &&
                   file.exists(datapath[[1]]))
    }

    rds_path <- rds_path %||% ""
    if (!nzchar(rds_path)) {
        return(FALSE)
    }

    file.exists(resolve_input_path(rds_path))
}


cluster_control_placeholder <- function(message = "Load a processed Seurat RDS to choose a root/start cluster.") {
    tags$p(class = "text-muted small mb-0", message)
}


trajectory_start_label <- function(method = "slingshot") {
    if (identical(method, "slingshot")) {
        "Start cluster:"
    } else {
        "Root/start cluster:"
    }
}


read_trajectory_gene_list <- function(path, header = TRUE) {
    df <- read_uploaded_table(
        path,
        header = header,
        stringsAsFactors = FALSE
    )
    if (!ncol(df)) {
        return(character())
    }

    genes <- trimws(as.character(df[[1]]))
    unique(genes[!is.na(genes) & nzchar(genes)])
}


validate_trajectory_heatmap_object <- function(obj) {
    metadata <- obj[[]]
    required <- c("trajectory_path", "pseudotime")
    missing <- setdiff(required, colnames(metadata))
    if (length(missing)) {
        stop(
            "The trajectory object must include trajectory_path and pseudotime metadata. ",
            "Use a trajectory object saved from this tab or run trajectory inference first.",
            call. = FALSE
        )
    }

    invisible(TRUE)
}


trajectory_heatmap_plot <- function(obj, genes) {
    validate_trajectory_heatmap_object(obj)

    genes <- unique(trimws(as.character(genes)))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    genes <- intersect(genes, rownames(obj))
    if (!length(genes)) {
        stop("None of the selected genes are present in the trajectory object.",
             call. = FALSE)
    }

    if ("RNA" %in% names(obj@assays)) {
        SeuratObject::DefaultAssay(obj) <- "RNA"
    }

    assay <- SeuratObject::DefaultAssay(obj)
    layers <- tryCatch(
        SeuratObject::Layers(obj[[assay]]),
        error = function(e) character()
    )
    data_layers <- grep("^data($|[.])", layers, value = TRUE)
    if (length(data_layers) > 1L) {
        obj <- tryCatch(
            SeuratObject::JoinLayers(obj, assay = assay),
            error = function(e) obj
        )
    }

    metadata <- obj[[]]
    cell_metadata <- data.frame(
        cell = colnames(obj),
        trajectory_path = as.character(metadata$trajectory_path),
        pseudotime = suppressWarnings(as.numeric(metadata$pseudotime)),
        stringsAsFactors = FALSE
    )
    cell_metadata <- cell_metadata[
        !is.na(cell_metadata$trajectory_path) &
            nzchar(cell_metadata$trajectory_path) &
            is.finite(cell_metadata$pseudotime),
        ,
        drop = FALSE
    ]
    if (nrow(cell_metadata) < 2L) {
        stop(
            "The trajectory object does not contain enough cells with finite pseudotime ",
            "and trajectory path assignments.",
            call. = FALSE
        )
    }

    cell_order <- cell_metadata$cell[order(
        cell_metadata$trajectory_path,
        cell_metadata$pseudotime,
        na.last = TRUE
    )]

    Seurat::DoHeatmap(
        obj,
        features = genes,
        cells = cell_order,
        slot = "data",
        group.by = "trajectory_path"
    ) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(title = NULL) +
        ggplot2::theme(plot.title = ggplot2::element_blank())
}


trajectory_feature_plot <- function(obj, genes, max_genes = 12) {
    obj <- prepare_trajectory_de_object(obj)
    genes <- unique(trimws(as.character(genes)))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    genes <- intersect(genes, rownames(obj))
    max_genes <- suppressWarnings(as.integer(max_genes %||% 12))
    if (!is.finite(max_genes) || max_genes < 1L) {
        max_genes <- 12L
    }
    genes <- utils::head(genes, max_genes)
    if (!length(genes)) {
        stop("None of the selected genes are present in the Seurat object.",
             call. = FALSE)
    }

    scCustomize::FeaturePlot_scCustom(
        seurat_object = obj,
        features = genes,
        order = TRUE
    )
}


resolve_start_cluster <- function(obj, start_cluster = NULL) {
    start_cluster <- selected_cluster_or_null(start_cluster)
    if (!is.null(start_cluster)) {
        return(start_cluster)
    }
    levels <- cluster_levels(obj)
    if (!length(levels)) {
        stop("No cluster identities were found in the Seurat object.", call. = FALSE)
    }
    levels[[1]]
}


resolve_reduction_name <- function(obj, reduction) {
    reductions <- names(obj@reductions)
    if (reduction %in% reductions) {
        return(reduction)
    }

    matched <- reductions[tolower(reductions) == tolower(reduction)]
    if (length(matched)) {
        return(matched[[1]])
    }

    stop("Required reduction not found in Seurat object: ", reduction, call. = FALSE)
}


get_reduction_embeddings <- function(obj, reduction) {
    SeuratObject::Embeddings(obj, resolve_reduction_name(obj, reduction))
}


get_counts_matrix <- function(obj) {
    assay <- if ("RNA" %in% names(obj@assays)) {
        "RNA"
    } else {
        SeuratObject::DefaultAssay(obj)
    }

    assay_obj <- obj[[assay]]
    layers <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character())
    count_layers <- grep("^counts($|[.])", layers, value = TRUE)
    if (length(count_layers) > 1L) {
        obj <- tryCatch(
            SeuratObject::JoinLayers(obj, assay = assay),
            error = function(e) {
                stop(
                    "This Seurat v5 object has multiple count layers in assay '",
                    assay,
                    "', but they could not be joined for trajectory analysis: ",
                    conditionMessage(e),
                    call. = FALSE
                )
            }
        )
    }

    counts <- tryCatch(
        Seurat::GetAssayData(obj, assay = assay, layer = "counts"),
        error = function(e) {
            Seurat::GetAssayData(obj, assay = assay, slot = "counts")
        }
    )

    if (!nrow(counts) || !ncol(counts)) {
        stop("Could not retrieve a non-empty counts matrix from assay '", assay, "'.",
             call. = FALSE)
    }

    counts
}


trajectory_sce_from_metadata <- function(obj) {
    umap <- get_reduction_embeddings(obj, "umap")
    cells <- colnames(obj)
    if (!all(cells %in% rownames(umap))) {
        stop(
            "The UMAP embedding does not contain all cells in this Seurat object. ",
            "Run clustering/UMAP again before trajectory inference.",
            call. = FALSE
        )
    }
    umap <- umap[cells, , drop = FALSE]

    metadata <- obj@meta.data
    if (!all(cells %in% rownames(metadata))) {
        if (nrow(metadata) == length(cells)) {
            rownames(metadata) <- cells
        } else {
            stop(
                "The Seurat metadata does not match the object cell barcodes.",
                call. = FALSE
            )
        }
    }
    metadata <- metadata[cells, , drop = FALSE]

    counts <- get_counts_matrix(obj)
    if (!all(cells %in% colnames(counts))) {
        stop(
            "The count matrix does not contain all trajectory cells. ",
            "Run normalization and clustering again before trajectory inference.",
            call. = FALSE
        )
    }
    counts <- counts[, cells, drop = FALSE]

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = metadata
    )
    SingleCellExperiment::reducedDim(sce, "UMAP") <- umap
    sce$ident <- unname(seurat_clusters(obj)[cells])
    sce
}


representative_cell_for_cluster <- function(embedding, clusters, cluster) {
    clusters <- as.character(clusters)
    cells <- rownames(embedding)[clusters == cluster]
    if (!length(cells)) {
        stop("Start cluster was not found in the trajectory embedding: ", cluster,
             call. = FALSE)
    }

    cluster_embedding <- embedding[cells, , drop = FALSE]
    centroid <- colMeans(cluster_embedding)
    distances <- rowSums(sweep(cluster_embedding, 2, centroid, "-")^2)
    cells[[which.min(distances)]]
}


assign_common_trajectory_metadata <- function(obj, method, pseudotime, path) {
    obj$trajectory_method <- method
    obj$pseudotime <- as.numeric(pseudotime)
    obj$trajectory_path <- as.character(path)
    obj
}


#' Extract Hard Lineage Assignments from Slingshot
#'
#' @param sce A SingleCellExperiment with Slingshot results.
#' @param path_prefix Prefix for lineage labels.
#' @return A list with pseudotime and lineage assignment metadata.
#' @keywords internal
extract_slingshot_assignments <- function(sce, path_prefix = "Path") {
    pseudotime <- as.matrix(slingshot::slingPseudotime(sce, na = FALSE))
    weights <- as.matrix(slingshot::slingCurveWeights(sce))

    if (!ncol(pseudotime)) {
        stop("Slingshot did not return any lineage pseudotime values.", call. = FALSE)
    }

    path_names <- paste0(path_prefix, seq_len(ncol(pseudotime)))
    colnames(pseudotime) <- path_names
    colnames(weights) <- path_names

    weights[!is.finite(weights)] <- -Inf
    assigned_index <- max.col(weights, ties.method = "first")
    assigned_path <- path_names[assigned_index]
    assigned_pseudotime <- pseudotime[cbind(seq_len(nrow(pseudotime)), assigned_index)]

    missing_pseudotime <- !is.finite(assigned_pseudotime)
    if (any(missing_pseudotime)) {
        fallback_index <- apply(pseudotime[missing_pseudotime, , drop = FALSE], 1, function(values) {
            valid <- which(is.finite(values))
            if (!length(valid)) {
                return(NA_integer_)
            }
            valid[[1]]
        })

        keep_fallback <- !is.na(fallback_index)
        missing_rows <- which(missing_pseudotime)
        assigned_index[missing_rows[keep_fallback]] <- fallback_index[keep_fallback]
        assigned_path[missing_rows[keep_fallback]] <- path_names[fallback_index[keep_fallback]]
        assigned_pseudotime[missing_rows[keep_fallback]] <- pseudotime[
            cbind(missing_rows[keep_fallback], fallback_index[keep_fallback])
        ]
    }

    keep_cells <- !is.na(assigned_path) & is.finite(assigned_pseudotime)

    list(
        path_names = path_names,
        pseudotime = pseudotime,
        weights = weights,
        assigned_index = assigned_index,
        assigned_path = assigned_path,
        assigned_pseudotime = assigned_pseudotime,
        keep_cells = keep_cells
    )
}


#' Run scMaSigPro on Slingshot Output
#'
#' @param sce A SingleCellExperiment with trajectory_path and trajectory_pseudotime columns.
#' @param gene_subset Optional character vector of genes to test.
#' @param p_value Numeric significance cutoff.
#' @param rsq Numeric R-squared cutoff.
#' @param poly_degree Integer polynomial degree for the model.
#' @return A list containing the fitted scMaSigPro object and a result table.
#' @keywords internal
run_scmasigpro <- function(sce, gene_subset = NULL,
                           p_value = 0.05, rsq = 0.6, poly_degree = 2) {
    if (!requireNamespace("scMaSigPro", quietly = TRUE)) {
        stop(
            "Package 'scMaSigPro' is not installed. Install with: remotes::install_github('BioBam/scMaSigPro')",
            call. = FALSE
        )
    }

    keep_cells <- !is.na(sce$trajectory_path) &
        is.finite(sce$trajectory_pseudotime)
    sce <- sce[, keep_cells]

    if (length(unique(as.character(sce$trajectory_path))) < 2) {
        stop("scMaSigPro requires at least two Slingshot lineages.", call. = FALSE)
    }

    scmasigpro_sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = SingleCellExperiment::counts(sce))
    )
    scmasigpro_sce$trajectory_path <- as.character(sce$trajectory_path)
    scmasigpro_sce$trajectory_pseudotime <- as.numeric(sce$trajectory_pseudotime)

    if (!is.null(gene_subset)) {
        gene_subset <- unique(intersect(as.character(gene_subset), rownames(scmasigpro_sce)))
        if (length(gene_subset)) {
            scmasigpro_sce <- scmasigpro_sce[gene_subset, ]
        }
    }

    scmp_obj <- scMaSigPro::as_scmp(
        object = scmasigpro_sce,
        from = "sce",
        ptime_col = "trajectory_pseudotime",
        path_col = "trajectory_path",
        verbose = FALSE,
        additional_params = list(
            labels_exist = TRUE,
            exist_ptime_col = "trajectory_pseudotime",
            exist_path_col = "trajectory_path"
        )
    )
    empty_scmasigpro_result <- function(message) {
        list(
            engine = "scMaSigPro",
            scmp = scmp_obj,
            table = data.frame(
                gene = character(),
                p_value = numeric(),
                rsquared = numeric(),
                lineages = character(),
                n_lineages = integer(),
                stringsAsFactors = FALSE
            ),
            message = message
        )
    }

    scmp_obj <- scMaSigPro::sc.squeeze(
        scmp_obj,
        ptime_col = "trajectory_pseudotime",
        path_col = "trajectory_path",
        verbose = FALSE,
        additional_params = list(use_unique_time_points = TRUE)
    )
    scmp_obj <- scMaSigPro::sc.set.poly(
        scmp_obj,
        poly_degree = poly_degree,
        path_col = "trajectory_path"
    )
    scmasigpro_error <- NULL
    scmp_obj <- tryCatch({
        obj_fit <- scMaSigPro::sc.p.vector(
            scmp_obj,
            p_value = p_value,
            verbose = FALSE,
            parallel = FALSE
        )
        obj_fit <- scMaSigPro::sc.t.fit(
            obj_fit,
            p_value = p_value,
            verbose = FALSE,
            parallel = FALSE
        )
        scMaSigPro::sc.filter(
            obj_fit,
            rsq = rsq,
            p_value = p_value,
            vars = "groups",
            includeInflu = TRUE
        )
    }, error = function(e) {
        scmasigpro_error <<- conditionMessage(e)
        scmp_obj
    })
    if (!is.null(scmasigpro_error)) {
        return(empty_scmasigpro_result(
            paste(
                "scMaSigPro completed, but no significant branching genes were detected",
                "with the current settings."
            )
        ))
    }

    solution_table <- tryCatch(
        scMaSigPro::showSol(
            scmp_obj,
            view = FALSE,
            return = TRUE,
            includeInflu = TRUE
        ),
        error = function(e) data.frame()
    )
    solution_table <- as.data.frame(solution_table)
    significant_genes <- tryCatch(
        scmp_obj@Significant@genes,
        error = function(e) list()
    )

    path_membership <- do.call(rbind, lapply(names(significant_genes), function(path_name) {
        genes <- significant_genes[[path_name]]
        if (!length(genes)) {
            return(NULL)
        }

        data.frame(
            gene = genes,
            lineage = path_name,
            stringsAsFactors = FALSE
        )
    }))

    if (is.null(path_membership) || !nrow(path_membership)) {
        return(empty_scmasigpro_result(
            paste(
                "scMaSigPro completed, but no significant branching genes were detected",
                "with the current settings."
            )
        ))
    }

    lineage_summary <- stats::aggregate(
        lineage ~ gene,
        data = path_membership,
        FUN = function(values) paste(unique(values), collapse = ", ")
    )
    lineage_counts <- stats::aggregate(
        lineage ~ gene,
        data = path_membership,
        FUN = function(values) length(unique(values))
    )
    colnames(lineage_summary)[2] <- "lineages"
    colnames(lineage_counts)[2] <- "n_lineages"

    result <- data.frame(
        gene = rownames(solution_table),
        stringsAsFactors = FALSE
    )
    result$p_value <- if ("p-value" %in% colnames(solution_table)) {
        solution_table[["p-value"]]
    } else {
        NA_real_
    }
    result$rsquared <- if ("R-squared" %in% colnames(solution_table)) {
        solution_table[["R-squared"]]
    } else {
        NA_real_
    }

    result <- merge(result, lineage_summary, by = "gene", all.x = TRUE)
    result <- merge(result, lineage_counts, by = "gene", all.x = TRUE)
    result <- result[result$gene %in% unique(path_membership$gene), , drop = FALSE]
    result <- result[order(result$p_value, -result$rsquared, result$gene), , drop = FALSE]
    rownames(result) <- NULL

    list(
        engine = "scMaSigPro",
        scmp = scmp_obj,
        table = result
    )
}


#' Run PseudotimeDE-fast on Slingshot Output
#'
#' @param sce A SingleCellExperiment with trajectory_path and trajectory_pseudotime columns.
#' @param gene_subset Optional character vector of genes to test.
#' @param cores Number of workers to use.
#' @return A list containing the engine name and result table.
#' @keywords internal
run_pseudotimede_fast <- function(sce, gene_subset = NULL, cores = 2) {
    assert_pseudotimede_available()

    keep_cells <- !is.na(sce$trajectory_path) &
        is.finite(sce$trajectory_pseudotime)
    sce <- sce[, keep_cells]

    if (ncol(sce) < 5) {
        stop("PseudotimeDE-fast requires at least five cells with finite pseudotime.",
             call. = FALSE)
    }

    if (!is.null(gene_subset)) {
        gene_subset <- unique(intersect(as.character(gene_subset), rownames(sce)))
        if (length(gene_subset)) {
            sce <- sce[gene_subset, ]
        }
    }

    cores <- as.integer(cores %||% 2L)
    if (is.na(cores) || cores < 1L) {
        cores <- 1L
    }

    ori_tbl <- data.frame(
        cell = colnames(sce),
        pseudotime = as.numeric(sce$trajectory_pseudotime),
        stringsAsFactors = FALSE
    )

    result <- as.data.frame(PseudotimeDE::runTauStarDE(
        gene.vec = rownames(sce),
        ori.tbl = ori_tbl,
        mat = SingleCellExperiment::counts(sce),
        assay.use = "counts",
        mode = "asymptotic",
        mc.cores = cores,
        SIMPLIFY = TRUE
    ))

    if ("taustar.pv" %in% colnames(result)) {
        result$p_value <- as.numeric(result$taustar.pv)
        result$adjusted_p_value <- stats::p.adjust(result$p_value, method = "BH")
        result <- result[order(result$p_value, result$gene), , drop = FALSE]
    }
    rownames(result) <- NULL

    list(
        engine = "PseudotimeDE-fast",
        table = result
    )
}


pseudotimede_available <- function() {
    requireNamespace("PseudotimeDE", quietly = TRUE)
}


pseudotimede_install_message <- function() {
    paste(
        "PseudotimeDE-fast is not installed; it is required for trajectory differential expression.",
        "The Docker image includes it.",
        "On macOS, install the official R GNU Fortran toolchain first: https://mac.r-project.org/tools/.",
        "For a local R package install, run:",
        "remotes::install_git('https://github.com/dsong-lab/PseudotimeDE.git')"
    )
}


assert_pseudotimede_available <- function() {
    if (!pseudotimede_available()) {
        stop(pseudotimede_install_message(), call. = FALSE)
    }
    invisible(TRUE)
}


#' Run tradeSeq on Slingshot Output
#'
#' @param sce A SingleCellExperiment with Slingshot results.
#' @param assignment Slingshot assignment metadata from extract_slingshot_assignments().
#' @param gene_subset Optional character vector of genes to test.
#' @param nknots Number of knots for tradeSeq GAM fitting.
#' @return A list containing the fitted tradeSeq object and a result table.
#' @keywords internal
run_tradeseq <- function(sce, assignment, gene_subset = NULL, nknots = 6) {
    if (!requireNamespace("tradeSeq", quietly = TRUE)) {
        stop(
            "Package 'tradeSeq' is not installed. Install with: BiocManager::install('tradeSeq')",
            call. = FALSE
        )
    }

    keep_cells <- assignment$keep_cells
    sce <- sce[, keep_cells]
    pseudotime <- assignment$pseudotime[keep_cells, , drop = FALSE]
    cell_weights <- assignment$weights[keep_cells, , drop = FALSE]
    cell_weights[!is.finite(cell_weights)] <- 0
    cell_weights[cell_weights < 0] <- 0

    if (!is.null(gene_subset)) {
        gene_subset <- unique(intersect(as.character(gene_subset), rownames(sce)))
        if (length(gene_subset)) {
            sce <- sce[gene_subset, ]
        }
    }

    fit <- tradeSeq::fitGAM(
        counts = SingleCellExperiment::counts(sce),
        pseudotime = pseudotime,
        cellWeights = cell_weights,
        nknots = as.integer(nknots),
        verbose = FALSE
    )
    association <- as.data.frame(tradeSeq::associationTest(fit))

    result <- data.frame(
        gene = rownames(association),
        association,
        row.names = NULL,
        check.names = FALSE
    )
    p_col <- intersect(c("pvalue", "p.value", "p_value"), colnames(result))
    if (length(p_col)) {
        result$p_value <- result[[p_col[[1]]]]
        result <- result[order(result$p_value, result$gene), , drop = FALSE]
    }
    rownames(result) <- NULL

    list(
        engine = "tradeSeq",
        fit = fit,
        table = result
    )
}


#' Check Required Trajectory Dependencies
#'
#' @param method Trajectory method identifier.
#' @return Invisibly TRUE when trajectory dependencies are available.
#' @keywords internal
check_trajectory_dependencies <- function(method = "slingshot") {
    required <- switch(
        method,
        slingshot = c("slingshot", "SingleCellExperiment", "DelayedMatrixStats"),
        paga = "reticulate",
        monocle3 = "monocle3",
        c("slingshot", "SingleCellExperiment", "DelayedMatrixStats")
    )
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
        stop(
            "Trajectory inference requires missing package(s): ",
            paste(missing, collapse = ", "),
            ". Install the full Asc-Seurat dependency stack or use the Docker image.",
            call. = FALSE
        )
    }

    if (identical(method, "paga")) {
        declare_paga_python_requirements()
        if (!paga_available()) {
            stop(
                "PAGA is not available. Run ascseurat::setup_paga() ",
                "or use the Asc-Seurat Docker image.",
                call. = FALSE
            )
        }
    }

    invisible(TRUE)
}


#' Run Slingshot Trajectory Inference
#'
#' @keywords internal
run_slingshot_trajectory <- function(obj, start_cluster = NULL, end_cluster = NULL) {
    check_trajectory_dependencies("slingshot")

    start_clust <- selected_cluster_or_null(start_cluster)
    end_clust <- selected_cluster_or_null(end_cluster)

    sce <- trajectory_sce_from_metadata(obj)

    sce <- slingshot::slingshot(
        sce,
        clusterLabels = "ident",
        reducedDim = "UMAP",
        start.clus = start_clust,
        end.clus = end_clust
    )

    assignment <- extract_slingshot_assignments(sce)
    obj <- assign_common_trajectory_metadata(
        obj,
        method = "Slingshot",
        pseudotime = assignment$assigned_pseudotime,
        path = assignment$assigned_path
    )
    for (i in seq_along(assignment$path_names)) {
        obj[[paste0("pseudotime_", assignment$path_names[[i]])]] <- assignment$pseudotime[, i]
    }

    sce$trajectory_path <- assignment$assigned_path
    sce$trajectory_pseudotime <- assignment$assigned_pseudotime
    sds <- slingshot::SlingshotDataSet(sce)

    list(
        obj = obj,
        sce = sce,
        sds = sds,
        assignment = assignment,
        method = "slingshot",
        method_label = "Slingshot",
        root_cluster = start_clust,
        end_cluster = end_clust
    )
}


#' Run PAGA Trajectory Inference
#'
#' @keywords internal
run_paga_trajectory <- function(obj, start_cluster = NULL, end_cluster = NULL) {
    check_trajectory_dependencies("paga")

    pca <- get_reduction_embeddings(obj, "pca")
    umap <- get_reduction_embeddings(obj, "umap")
    cells <- intersect(rownames(pca), rownames(umap))
    if (length(cells) < 3) {
        stop("PAGA requires at least three cells with PCA and UMAP embeddings.",
             call. = FALSE)
    }
    pca <- pca[cells, , drop = FALSE]
    umap <- umap[cells, , drop = FALSE]
    clusters <- seurat_clusters(obj)[cells]
    cluster_order <- sort_cluster_labels(unique(clusters))
    root_cluster <- resolve_start_cluster(obj, start_cluster)
    root_cell <- representative_cell_for_cluster(pca, clusters, root_cluster)
    root_index <- match(root_cell, cells) - 1L

    scanpy <- reticulate::import("scanpy", convert = FALSE)
    anndata <- reticulate::import("anndata", convert = FALSE)
    pandas <- reticulate::import("pandas", convert = FALSE)

    obs <- pandas$DataFrame(
        reticulate::dict(seurat_cluster = clusters),
        index = cells
    )
    adata <- anndata$AnnData(
        X = reticulate::r_to_py(unname(pca)),
        obs = obs
    )
    adata$obsm[["X_umap"]] <- reticulate::r_to_py(unname(umap))
    cluster_dtype <- pandas$api$types$CategoricalDtype(
        categories = reticulate::r_to_py(cluster_order),
        ordered = FALSE
    )
    adata$obs[["seurat_cluster"]] <- adata$obs[["seurat_cluster"]]$astype(cluster_dtype)

    n_neighbors <- as.integer(min(15, max(2, nrow(pca) - 1L)))
    n_dcs <- as.integer(min(10, ncol(pca), nrow(pca) - 1L))
    if (n_dcs < 2) {
        stop("PAGA requires at least two usable PCA dimensions.", call. = FALSE)
    }

    scanpy$pp$neighbors(adata, n_neighbors = n_neighbors, use_rep = "X")
    scanpy$tl$paga(adata, groups = "seurat_cluster")
    scanpy$tl$diffmap(adata, n_comps = n_dcs)
    adata$uns[["iroot"]] <- as.integer(root_index)
    scanpy$tl$dpt(adata, n_dcs = n_dcs, n_branchings = 0L)

    pseudotime <- as.numeric(reticulate::py_to_r(
        adata$obs[["dpt_pseudotime"]]$to_numpy()
    ))
    names(pseudotime) <- cells

    paga_uns <- adata$uns[["paga"]]
    connectivities_py <- paga_uns[["connectivities"]]
    connectivities <- tryCatch(
        reticulate::py_to_r(connectivities_py$toarray()),
        error = function(e) reticulate::py_to_r(connectivities_py)
    )
    connectivities <- as.matrix(connectivities)
    dimnames(connectivities) <- list(cluster_order, cluster_order)

    full_pseudotime <- rep(NA_real_, ncol(obj))
    names(full_pseudotime) <- colnames(obj)
    full_pseudotime[cells] <- pseudotime
    full_path <- rep(NA_character_, ncol(obj))
    names(full_path) <- colnames(obj)
    full_path[cells] <- clusters

    obj <- assign_common_trajectory_metadata(
        obj,
        method = "PAGA",
        pseudotime = full_pseudotime,
        path = full_path
    )
    obj$paga_cluster <- full_path

    list(
        obj = obj,
        adata = adata,
        method = "paga",
        method_label = "PAGA",
        root_cluster = root_cluster,
        root_cell = root_cell,
        end_cluster = NULL,
        paga = list(
            clusters = cluster_order,
            connectivities = connectivities,
            umap = umap,
            cell_clusters = clusters
        )
    )
}


#' Run Monocle 3 Trajectory Inference
#'
#' @keywords internal
run_monocle3_trajectory <- function(obj, start_cluster = NULL, end_cluster = NULL) {
    check_trajectory_dependencies("monocle3")

    counts <- get_counts_matrix(obj)
    counts <- counts[, colnames(obj), drop = FALSE]
    cell_metadata <- obj@meta.data
    cell_metadata$seurat_cluster <- seurat_clusters(obj)[rownames(cell_metadata)]
    gene_metadata <- data.frame(
        gene_short_name = rownames(counts),
        row.names = rownames(counts),
        stringsAsFactors = FALSE
    )

    cds <- monocle3::new_cell_data_set(
        counts,
        cell_metadata = cell_metadata,
        gene_metadata = gene_metadata
    )
    num_dim <- min(50L, ncol(cds) - 1L, nrow(cds) - 1L)
    if (num_dim < 2L) {
        stop("Monocle 3 requires at least two dimensions for preprocessing.",
             call. = FALSE)
    }

    cds <- monocle3::preprocess_cds(cds, num_dim = as.integer(num_dim))
    umap <- get_reduction_embeddings(obj, "umap")
    umap <- umap[colnames(cds), , drop = FALSE]
    SingleCellExperiment::reducedDims(cds)$UMAP <- umap
    cds <- monocle3::cluster_cells(
        cds,
        reduction_method = "UMAP",
        cluster_method = "louvain",
        k = as.integer(min(20L, ncol(cds) - 1L)),
        verbose = FALSE
    )
    cds <- monocle3::learn_graph(cds, use_partition = TRUE, verbose = FALSE)

    root_cluster <- resolve_start_cluster(obj, start_cluster)
    root_cells <- names(which(seurat_clusters(obj) == root_cluster))
    if (!length(root_cells)) {
        stop("Start cluster was not found for Monocle 3: ", root_cluster,
             call. = FALSE)
    }
    cds <- monocle3::order_cells(
        cds,
        reduction_method = "UMAP",
        root_cells = root_cells
    )

    pseudotime <- as.numeric(monocle3::pseudotime(cds))
    names(pseudotime) <- colnames(cds)
    partitions <- tryCatch(
        as.character(monocle3::partitions(cds, reduction_method = "UMAP")),
        error = function(e) rep("1", ncol(cds))
    )
    names(partitions) <- colnames(cds)

    obj <- assign_common_trajectory_metadata(
        obj,
        method = "Monocle 3",
        pseudotime = pseudotime[colnames(obj)],
        path = paste0("Partition ", partitions[colnames(obj)])
    )
    obj$monocle3_partition <- partitions[colnames(obj)]

    list(
        obj = obj,
        cds = cds,
        method = "monocle3",
        method_label = "Monocle 3",
        root_cluster = root_cluster,
        root_cells = root_cells,
        end_cluster = NULL
    )
}


strongest_paga_connections <- function(connectivities, n = 10, min_weight = 0.05) {
    if (is.null(connectivities) || !length(connectivities)) {
        return(data.frame())
    }

    connectivities <- as.matrix(connectivities)
    edges <- which(upper.tri(connectivities) & connectivities >= min_weight,
                   arr.ind = TRUE)
    if (!nrow(edges)) {
        return(data.frame())
    }

    result <- data.frame(
        from = rownames(connectivities)[edges[, 1]],
        to = colnames(connectivities)[edges[, 2]],
        connectivity = connectivities[edges],
        stringsAsFactors = FALSE
    )
    result <- result[order(-result$connectivity, result$from, result$to), , drop = FALSE]
    utils::head(result, n)
}


paga_connections_table <- function(connectivities, min_weight = 0.05) {
    if (is.null(connectivities) || !length(connectivities)) {
        return(data.frame(
            edge = character(),
            source_cluster = character(),
            target_cluster = character(),
            connectivity = numeric(),
            stringsAsFactors = FALSE
        ))
    }

    connectivities <- as.matrix(connectivities)
    edges <- which(upper.tri(connectivities) & connectivities >= min_weight,
                   arr.ind = TRUE)
    if (!nrow(edges)) {
        return(data.frame(
            edge = character(),
            source_cluster = character(),
            target_cluster = character(),
            connectivity = numeric(),
            stringsAsFactors = FALSE
        ))
    }

    result <- data.frame(
        source_cluster = rownames(connectivities)[edges[, 1]],
        target_cluster = colnames(connectivities)[edges[, 2]],
        connectivity = as.numeric(connectivities[edges]),
        stringsAsFactors = FALSE
    )
    result <- result[order(-result$connectivity, result$source_cluster,
                           result$target_cluster), , drop = FALSE]
    result$edge <- paste(result$source_cluster, result$target_cluster, sep = "||")
    result <- result[, c("edge", "source_cluster", "target_cluster", "connectivity"),
                     drop = FALSE]
    rownames(result) <- NULL
    result
}


trajectory_gene_subset <- function(obj, mode = "variable") {
    mode <- mode %||% "variable"
    if (identical(mode, "all")) {
        return(rownames(obj))
    }

    genes <- Seurat::VariableFeatures(obj)
    genes <- intersect(genes, rownames(obj))
    if (!length(genes)) {
        genes <- rownames(obj)
    }
    genes
}


prepare_trajectory_de_object <- function(obj) {
    if ("RNA" %in% names(obj@assays)) {
        SeuratObject::DefaultAssay(obj) <- "RNA"
    }
    tryCatch(
        SeuratObject::JoinLayers(obj),
        error = function(e) obj
    )
}


paga_edge_choices <- function(res) {
    edges <- paga_connections_table(res$paga$connectivities)
    if (!nrow(edges)) {
        return(c("No connected PAGA edges available" = ""))
    }

    labels <- sprintf(
        "%s -> %s (%.3f)",
        edges$source_cluster,
        edges$target_cluster,
        edges$connectivity
    )
    stats::setNames(edges$edge, labels)
}


selected_paga_edges <- function(res, scope = "selected", selected_edge = NULL,
                                top_n = 5) {
    edges <- paga_connections_table(res$paga$connectivities)
    if (!nrow(edges)) {
        stop("No connected PAGA edges were found for gene discovery.",
             call. = FALSE)
    }

    if (identical(scope, "top")) {
        top_n <- as.integer(top_n %||% 5L)
        if (is.na(top_n) || top_n < 1L) {
            top_n <- 5L
        }
        return(utils::head(edges, top_n))
    }

    selected_edge <- selected_edge %||% edges$edge[[1]]
    selected <- edges[edges$edge == selected_edge, , drop = FALSE]
    if (!nrow(selected)) {
        stop("Selected PAGA edge is no longer available.", call. = FALSE)
    }
    selected
}


format_paga_marker_table <- function(markers, edge) {
    if (!nrow(markers)) {
        return(data.frame())
    }

    logfc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))
    if (!length(logfc_col)) {
        markers$avg_log2FC <- NA_real_
        logfc_col <- "avg_log2FC"
    }
    p_col <- intersect(c("p_val", "p.value", "p_value"), colnames(markers))
    adj_col <- intersect(c("p_val_adj", "p.adjust", "adjusted_p_value"),
                         colnames(markers))

    result <- data.frame(
        edge = edge$edge,
        source_cluster = edge$source_cluster,
        target_cluster = edge$target_cluster,
        gene = rownames(markers),
        direction = ifelse(
            markers[[logfc_col[[1]]]] >= 0,
            paste0("enriched in ", edge$source_cluster),
            paste0("enriched in ", edge$target_cluster)
        ),
        avg_log2FC = as.numeric(markers[[logfc_col[[1]]]]),
        p_value = if (length(p_col)) as.numeric(markers[[p_col[[1]]]]) else NA_real_,
        p_val_adj = if (length(adj_col)) as.numeric(markers[[adj_col[[1]]]]) else NA_real_,
        pct_source = if ("pct.1" %in% colnames(markers)) as.numeric(markers$pct.1) else NA_real_,
        pct_target = if ("pct.2" %in% colnames(markers)) as.numeric(markers$pct.2) else NA_real_,
        connectivity = edge$connectivity,
        stringsAsFactors = FALSE
    )
    result <- result[order(result$p_val_adj, -abs(result$avg_log2FC),
                           result$gene, na.last = TRUE), , drop = FALSE]
    rownames(result) <- NULL
    result
}


empty_paga_gene_table <- function() {
    data.frame(
        edge = character(),
        source_cluster = character(),
        target_cluster = character(),
        gene = character(),
        direction = character(),
        avg_log2FC = numeric(),
        p_value = numeric(),
        p_val_adj = numeric(),
        pct_source = numeric(),
        pct_target = numeric(),
        connectivity = numeric(),
        rho = numeric(),
        stringsAsFactors = FALSE
    )
}


run_paga_edge_markers <- function(res, edges, gene_subset = NULL,
                                  min_pct = 0.05) {
    obj <- prepare_trajectory_de_object(res$obj)
    paga_clusters <- if ("paga_cluster" %in% colnames(res$obj@meta.data)) {
        res$obj$paga_cluster
    } else {
        seurat_clusters(res$obj)
    }
    obj$paga_cluster <- as.character(paga_clusters)
    SeuratObject::Idents(obj) <- factor(obj$paga_cluster)

    if (!is.null(gene_subset)) {
        gene_subset <- unique(intersect(as.character(gene_subset), rownames(obj)))
    }
    if (!length(gene_subset)) {
        gene_subset <- rownames(obj)
    }

    tables <- lapply(seq_len(nrow(edges)), function(i) {
        edge <- edges[i, , drop = FALSE]
        markers <- Seurat::FindMarkers(
            obj,
            ident.1 = edge$source_cluster,
            ident.2 = edge$target_cluster,
            features = gene_subset,
            test.use = "wilcox",
            logfc.threshold = 0,
            min.pct = min_pct,
            only.pos = FALSE,
            verbose = FALSE
        )
        format_paga_marker_table(markers, edge)
    })
    result <- do.call(rbind, tables)
    if (is.null(result)) {
        result <- empty_paga_gene_table()
    } else if (!"rho" %in% colnames(result)) {
        result$rho <- NA_real_
    }
    rownames(result) <- NULL

    list(
        engine = "PAGA connected cluster markers",
        mode = "edge_markers",
        edges = edges,
        table = result
    )
}


trajectory_expression_matrix <- function(obj, genes, cells, layer = "data") {
    genes <- unique(intersect(as.character(genes), rownames(obj)))
    cells <- unique(intersect(as.character(cells), colnames(obj)))
    if (!length(genes) || !length(cells)) {
        stop("No matching genes or cells were available for expression testing.",
             call. = FALSE)
    }

    expr <- tryCatch(
        Seurat::FetchData(
            object = obj,
            vars = genes,
            cells = cells,
            layer = layer,
            clean = FALSE
        ),
        error = function(e) {
            Seurat::FetchData(
                object = obj,
                vars = genes,
                cells = cells,
                layer = "counts",
                clean = FALSE
            )
        }
    )
    genes <- intersect(genes, colnames(expr))
    if (!length(genes)) {
        stop("The selected genes are not available in this Seurat object.",
             call. = FALSE)
    }

    expr <- as.matrix(expr[, genes, drop = FALSE])
    storage.mode(expr) <- "numeric"
    expr
}


run_paga_pseudotime_genes <- function(res, edges, gene_subset = NULL) {
    obj <- prepare_trajectory_de_object(res$obj)
    paga_clusters <- if ("paga_cluster" %in% colnames(res$obj@meta.data)) {
        res$obj$paga_cluster
    } else {
        seurat_clusters(res$obj)
    }
    obj$paga_cluster <- as.character(paga_clusters)
    pseudotime <- as.numeric(obj$pseudotime)
    names(pseudotime) <- colnames(obj)
    clusters <- stats::setNames(as.character(obj$paga_cluster), colnames(obj))

    if (!is.null(gene_subset)) {
        gene_subset <- unique(intersect(as.character(gene_subset), rownames(obj)))
    }
    if (!length(gene_subset)) {
        gene_subset <- rownames(obj)
    }

    tables <- lapply(seq_len(nrow(edges)), function(i) {
        edge <- edges[i, , drop = FALSE]
        edge_cells <- colnames(obj)[clusters %in%
            c(edge$source_cluster, edge$target_cluster)]
        edge_cells <- edge_cells[is.finite(pseudotime[edge_cells])]
        if (length(edge_cells) < 5L) {
            return(NULL)
        }
        if (length(unique(pseudotime[edge_cells])) < 3L) {
            return(NULL)
        }

        expr <- trajectory_expression_matrix(obj, gene_subset, edge_cells, layer = "data")
        genes <- colnames(expr)
        pt <- pseudotime[rownames(expr)]
        rho <- vapply(seq_along(genes), function(j) {
            suppressWarnings(stats::cor(expr[, j], pt, method = "spearman",
                                        use = "complete.obs"))
        }, numeric(1))
        rho[!is.finite(rho)] <- NA_real_
        n <- sum(is.finite(pt))
        denom <- pmax(1 - rho^2, .Machine$double.eps)
        statistic <- rho * sqrt((n - 2) / denom)
        p_value <- 2 * stats::pt(-abs(statistic), df = n - 2)
        p_value[!is.finite(p_value)] <- NA_real_

        source_cells <- edge_cells[clusters[edge_cells] == edge$source_cluster]
        target_cells <- edge_cells[clusters[edge_cells] == edge$target_cluster]
        pct_source <- if (length(source_cells)) {
            colMeans(expr[source_cells, , drop = FALSE] > 0, na.rm = TRUE)
        } else {
            rep(NA_real_, length(genes))
        }
        pct_target <- if (length(target_cells)) {
            colMeans(expr[target_cells, , drop = FALSE] > 0, na.rm = TRUE)
        } else {
            rep(NA_real_, length(genes))
        }

        result <- data.frame(
            edge = edge$edge,
            source_cluster = edge$source_cluster,
            target_cluster = edge$target_cluster,
            gene = genes,
            direction = ifelse(rho >= 0, "increases with pseudotime",
                               "decreases with pseudotime"),
            avg_log2FC = NA_real_,
            p_value = p_value,
            p_val_adj = stats::p.adjust(p_value, method = "BH"),
            pct_source = as.numeric(pct_source),
            pct_target = as.numeric(pct_target),
            connectivity = edge$connectivity,
            rho = rho,
            stringsAsFactors = FALSE
        )
        result[order(result$p_val_adj, -abs(result$rho), result$gene,
                     na.last = TRUE), , drop = FALSE]
    })

    result <- do.call(rbind, tables)
    if (is.null(result)) {
        result <- empty_paga_gene_table()
    }
    rownames(result) <- NULL

    list(
        engine = "PAGA pseudotime-associated genes",
        mode = "pseudotime",
        edges = edges,
        table = result
    )
}


run_monocle3_graph_test <- function(res, gene_subset = NULL, cores = 2) {
    if (!requireNamespace("monocle3", quietly = TRUE)) {
        stop(
            "Monocle 3 gene discovery requires the monocle3 package.",
            call. = FALSE
        )
    }

    cds <- res$cds
    if (!is.null(gene_subset)) {
        gene_subset <- unique(intersect(as.character(gene_subset), rownames(cds)))
        if (length(gene_subset)) {
            cds <- cds[gene_subset, ]
        }
    }

    graph_result <- as.data.frame(monocle3::graph_test(
        cds,
        neighbor_graph = "principal_graph",
        cores = as.integer(cores %||% 2L)
    ))
    graph_result$gene <- rownames(graph_result)
    if (!"morans_I" %in% colnames(graph_result)) {
        graph_result$morans_I <- NA_real_
    }
    if (!"p_value" %in% colnames(graph_result)) {
        graph_result$p_value <- NA_real_
    }
    if (!"q_value" %in% colnames(graph_result)) {
        graph_result$q_value <- stats::p.adjust(graph_result$p_value, method = "BH")
    }
    graph_result$status <- ifelse(
        is.finite(graph_result$q_value) & graph_result$q_value <= 0.05,
        "significant",
        "not significant"
    )

    result <- graph_result[, c("gene", "morans_I", "p_value", "q_value", "status"),
                           drop = FALSE]
    result <- result[order(result$q_value, -result$morans_I, result$gene,
                           na.last = TRUE), , drop = FALSE]
    rownames(result) <- NULL

    list(
        engine = "Monocle 3 graph_test",
        table = result
    )
}


run_monocle3_gene_modules <- function(res, genes, cores = 2) {
    if (!requireNamespace("monocle3", quietly = TRUE)) {
        stop(
            "Monocle 3 gene modules require the monocle3 package.",
            call. = FALSE
        )
    }

    genes <- unique(intersect(as.character(genes), rownames(res$cds)))
    if (length(genes) < 2L) {
        stop("At least two trajectory genes are required to find modules.",
             call. = FALSE)
    }

    modules <- as.data.frame(monocle3::find_gene_modules(
        res$cds[genes, ],
        resolution = 1e-2,
        cores = as.integer(cores %||% 2L)
    ))
    gene_col <- intersect(c("gene_id", "id", "gene_short_name"), colnames(modules))
    if (!length(gene_col)) {
        modules$gene_id <- rownames(modules)
        gene_col <- "gene_id"
    }
    if (!"module" %in% colnames(modules)) {
        modules$module <- "Module 1"
    }

    summary <- stats::aggregate(
        modules[[gene_col[[1]]]],
        by = list(module = modules$module),
        FUN = function(values) paste(utils::head(unique(values), 10), collapse = ", ")
    )
    colnames(summary)[2] <- "top_genes"
    counts <- stats::aggregate(
        modules[[gene_col[[1]]]],
        by = list(module = modules$module),
        FUN = function(values) length(unique(values))
    )
    colnames(counts)[2] <- "n_genes"
    summary <- merge(counts, summary, by = "module", all.x = TRUE)
    summary <- summary[order(summary$module), , drop = FALSE]
    rownames(summary) <- NULL

    list(
        modules = modules,
        summary = summary
    )
}


plot_monocle_module_heatmap <- function(res, module_result,
                                        genes_per_module = 8,
                                        bins = 20) {
    obj <- prepare_trajectory_de_object(res$obj)
    modules <- module_result$modules
    if (is.null(modules) || !nrow(modules)) {
        stop("No Monocle 3 gene modules were available to plot.", call. = FALSE)
    }

    gene_col <- intersect(
        c("gene_id", "id", "gene_short_name", "gene"),
        colnames(modules)
    )
    if (!length(gene_col)) {
        stop("Monocle 3 module results did not include gene identifiers.",
             call. = FALSE)
    }
    if (!"module" %in% colnames(modules)) {
        modules$module <- "Module 1"
    }

    modules$gene <- as.character(modules[[gene_col[[1]]]])
    modules$module <- as.character(modules$module)
    modules <- modules[modules$gene %in% rownames(obj), , drop = FALSE]
    if (!nrow(modules)) {
        stop("None of the Monocle 3 module genes were found in the Seurat object.",
             call. = FALSE)
    }
    modules <- modules[order(modules$module, modules$gene), , drop = FALSE]

    genes <- unlist(
        tapply(
            modules$gene,
            modules$module,
            function(values) utils::head(unique(values), genes_per_module)
        ),
        use.names = FALSE
    )
    genes <- unique(intersect(genes, rownames(obj)))

    pseudotime <- as.numeric(obj$pseudotime)
    names(pseudotime) <- colnames(obj)
    cells <- names(pseudotime)[is.finite(pseudotime)]
    if (length(cells) < 5L || length(unique(pseudotime[cells])) < 3L) {
        stop("Pseudotime values are not informative enough for a module heatmap.",
             call. = FALSE)
    }

    probs <- seq(0, 1, length.out = as.integer(bins %||% 20L) + 1L)
    breaks <- unique(stats::quantile(
        pseudotime[cells],
        probs = probs,
        na.rm = TRUE,
        names = FALSE
    ))
    if (length(breaks) < 3L) {
        breaks <- unique(pretty(pseudotime[cells], n = min(5L, length(cells))))
    }
    labels <- paste0("PT", seq_len(length(breaks) - 1L))
    pt_bin <- cut(pseudotime[cells], breaks = breaks, include.lowest = TRUE,
                  labels = labels)
    cells <- cells[!is.na(pt_bin)]
    pt_bin <- pt_bin[!is.na(pt_bin)]
    if (!length(cells)) {
        stop("No cells could be assigned to pseudotime bins.", call. = FALSE)
    }

    expr <- trajectory_expression_matrix(obj, genes, cells, layer = "data")
    scaled <- scale(expr)
    scaled[!is.finite(scaled)] <- 0
    module_lookup <- modules$module[match(colnames(scaled), modules$gene)]
    gene_levels <- rev(colnames(scaled))

    long <- data.frame(
        gene = rep(colnames(scaled), each = nrow(scaled)),
        module = rep(module_lookup, each = nrow(scaled)),
        bin = rep(as.character(pt_bin), times = ncol(scaled)),
        value = as.vector(scaled),
        stringsAsFactors = FALSE
    )
    avg <- stats::aggregate(value ~ gene + module + bin, long, mean)
    avg$gene <- factor(avg$gene, levels = gene_levels)
    avg$bin <- factor(avg$bin, levels = labels)

    ggplot2::ggplot(avg, ggplot2::aes(x = bin, y = gene, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::facet_grid(module ~ ., scales = "free_y", space = "free_y") +
        ggplot2::scale_fill_gradient2(
            low = "#2166AC",
            mid = "white",
            high = "#B2182B",
            midpoint = 0,
            name = "Scaled expression"
        ) +
        ggplot2::labs(
            title = "Monocle 3 gene modules across pseudotime",
            x = "Pseudotime bin",
            y = NULL
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid = ggplot2::element_blank(),
            strip.text.y = ggplot2::element_text(angle = 0)
        )
}


plot_paga_graph <- function(res) {
    umap <- res$paga$umap
    clusters <- res$paga$cell_clusters
    cluster_order <- res$paga$clusters
    connectivities <- res$paga$connectivities
    colors <- grDevices::hcl.colors(length(cluster_order), palette = "Dark 3")
    names(colors) <- cluster_order

    plot(
        umap,
        col = colors[clusters],
        pch = 16,
        cex = 0.5,
        xlab = "",
        ylab = "",
        xaxt = "n",
        yaxt = "n",
        main = ""
    )

    centers <- do.call(rbind, lapply(cluster_order, function(cluster) {
        colMeans(umap[clusters == cluster, , drop = FALSE])
    }))
    rownames(centers) <- cluster_order

    edges <- which(connectivities > 0.05 & upper.tri(connectivities), arr.ind = TRUE)
    if (nrow(edges)) {
        edge_weights <- connectivities[edges]
        for (i in seq_len(nrow(edges))) {
            from <- rownames(connectivities)[edges[i, 1]]
            to <- colnames(connectivities)[edges[i, 2]]
            graphics::segments(
                centers[from, 1], centers[from, 2],
                centers[to, 1], centers[to, 2],
                lwd = 1 + 4 * edge_weights[[i]],
                col = grDevices::adjustcolor("black", alpha.f = 0.45)
            )
        }
    }

    graphics::points(centers, pch = 21, bg = "white", col = "black", cex = 1.7)
    graphics::text(centers, labels = rownames(centers), cex = 0.8)
}


slingshot_lineage_edges <- function(res) {
    lineages <- tryCatch(
        slingshot::slingLineages(res$sds),
        error = function(e) list()
    )
    empty_edges <- data.frame(
        from = character(),
        to = character(),
        lineage = character(),
        step = integer(),
        stringsAsFactors = FALSE
    )

    if (!length(lineages)) {
        return(list(lineages = lineages, edges = empty_edges, lineage_edges = empty_edges))
    }

    edges <- do.call(rbind, lapply(seq_along(lineages), function(i) {
        lineage <- as.character(lineages[[i]])
        lineage <- lineage[nzchar(lineage)]
        if (length(lineage) < 2) {
            return(NULL)
        }

        data.frame(
            from = head(lineage, -1),
            to = utils::tail(lineage, -1),
            lineage = paste("Lineage", i),
            step = seq_len(length(lineage) - 1L),
            stringsAsFactors = FALSE
        )
    }))

    if (is.null(edges) || !nrow(edges)) {
        return(list(lineages = lineages, edges = empty_edges, lineage_edges = empty_edges))
    }

    edge_key <- mapply(
        function(from, to) paste(sort(c(from, to)), collapse = "||"),
        edges$from,
        edges$to,
        USE.NAMES = FALSE
    )
    unique_edges <- edges[!duplicated(edge_key), c("from", "to"), drop = FALSE]
    unique_edges$lineage <- vapply(unique(edge_key), function(key) {
        paste(unique(edges$lineage[edge_key == key]), collapse = ", ")
    }, character(1))
    unique_edges$step <- seq_len(nrow(unique_edges))
    rownames(edges) <- NULL
    rownames(unique_edges) <- NULL

    list(lineages = lineages, edges = unique_edges, lineage_edges = edges)
}


slingshot_mst_layout <- function(res) {
    graph <- slingshot_lineage_edges(res)
    lineages <- graph$lineages
    edges <- graph$edges
    lineage_edges <- graph$lineage_edges

    cluster_counts <- table(as.character(SeuratObject::Idents(res$obj)))
    clusters <- sort_cluster_labels(unique(c(
        names(cluster_counts),
        unlist(lapply(lineages, as.character), use.names = FALSE),
        edges$from,
        edges$to
    )))
    clusters <- clusters[nzchar(clusters)]
    if (!length(clusters)) {
        return(list(
            nodes = data.frame(),
            segments = data.frame(),
            edges = edges
        ))
    }

    root <- res$root_cluster
    if (is.null(root) || !nzchar(as.character(root))) {
        lineage_roots <- unlist(lapply(lineages, function(lineage) {
            lineage <- as.character(lineage)
            if (length(lineage)) lineage[[1]] else character()
        }), use.names = FALSE)
        root <- if (length(lineage_roots)) lineage_roots[[1]] else clusters[[1]]
    }
    root <- as.character(root)
    clusters <- unique(c(root, clusters))

    directed_edges <- lineage_edges[, c("from", "to", "lineage"), drop = FALSE]
    directed_edges <- directed_edges[
        !duplicated(paste(directed_edges$from, directed_edges$to, sep = "||")),
        ,
        drop = FALSE
    ]
    if (!nrow(directed_edges) && nrow(edges)) {
        directed_edges <- edges[, c("from", "to", "lineage"), drop = FALSE]
    }

    parent <- stats::setNames(rep(NA_character_, length(clusters)), clusters)
    if (nrow(directed_edges)) {
        for (i in seq_len(nrow(directed_edges))) {
            child <- directed_edges$to[[i]]
            if (child %in% names(parent) && is.na(parent[[child]]) &&
                !identical(child, root)) {
                parent[[child]] <- directed_edges$from[[i]]
            }
        }
    }

    children <- split(names(parent)[!is.na(parent)], parent[!is.na(parent)])
    depth <- stats::setNames(rep(NA_integer_, length(clusters)), clusters)
    depth[[root]] <- 0L
    queue <- root
    while (length(queue)) {
        current <- queue[[1]]
        queue <- queue[-1]
        child_nodes <- children[[current]]
        if (!length(child_nodes)) {
            next
        }
        for (child in child_nodes) {
            depth[[child]] <- depth[[current]] + 1L
        }
        queue <- c(queue, child_nodes)
    }
    missing_depth <- names(depth)[is.na(depth)]
    if (length(missing_depth)) {
        depth[missing_depth] <- max(depth[!is.na(depth)], 0L) + seq_along(missing_depth)
    }

    y_pos <- stats::setNames(rep(NA_real_, length(clusters)), clusters)
    assign_y <- function(node) {
        child_nodes <- children[[node]]
        if (!length(child_nodes)) {
            return(y_pos[[node]])
        }
        child_y <- vapply(child_nodes, assign_y, numeric(1))
        y_pos[[node]] <<- mean(child_y)
        y_pos[[node]]
    }
    leaves <- names(parent)[!names(parent) %in% names(children)]
    lineage_leaves <- unique(vapply(lineages, function(lineage) {
        lineage <- as.character(lineage)
        lineage <- lineage[nzchar(lineage)]
        if (length(lineage)) {
            utils::tail(lineage, 1L)
        } else {
            NA_character_
        }
    }, character(1)))
    lineage_leaves <- lineage_leaves[!is.na(lineage_leaves)]
    leaves <- unique(c(lineage_leaves, sort_cluster_labels(leaves)))
    leaves <- leaves[leaves %in% names(y_pos)]
    y_pos[leaves] <- seq_along(leaves)
    invisible(assign_y(root))
    missing_y <- names(y_pos)[!is.finite(y_pos)]
    if (length(missing_y)) {
        y_pos[missing_y] <- max(y_pos[is.finite(y_pos)], 0) + seq_along(missing_y)
    }

    nodes <- data.frame(
        cluster = clusters,
        x = as.numeric(depth[clusters]),
        y = as.numeric(y_pos[clusters]),
        parent = unname(parent[clusters]),
        stringsAsFactors = FALSE
    )
    nodes$cells <- as.integer(cluster_counts[nodes$cluster])
    nodes$cells[is.na(nodes$cells)] <- 0L
    nodes$role <- ifelse(nodes$cluster == root, "Root",
                         ifelse(nodes$cluster %in% as.character(res$end_cluster),
                                "End", "Cluster"))
    rownames(nodes) <- NULL

    segments <- data.frame()
    lineage_paths <- data.frame()
    node_lookup <- stats::setNames(seq_len(nrow(nodes)), nodes$cluster)

    if (nrow(directed_edges)) {
        segments <- do.call(rbind, lapply(seq_len(nrow(directed_edges)), function(i) {
            from <- node_lookup[[directed_edges$from[[i]]]]
            to <- node_lookup[[directed_edges$to[[i]]]]
            if (is.null(from) || is.null(to) || is.na(from) || is.na(to)) {
                return(NULL)
            }

            data.frame(
                x = nodes$x[[from]],
                xend = nodes$x[[to]],
                y = nodes$y[[to]],
                yend = nodes$y[[to]],
                lineages = directed_edges$lineage[[i]],
                stringsAsFactors = FALSE
            )
        }))
    }
    if (is.null(segments)) {
        segments <- data.frame()
    }

    vertical_segments <- data.frame()
    if (length(children)) {
        vertical_segments <- do.call(rbind, lapply(names(children), function(parent_node) {
            child_nodes <- children[[parent_node]]
            child_nodes <- child_nodes[child_nodes %in% nodes$cluster]
            if (length(child_nodes) < 2L) {
                return(NULL)
            }
            parent_idx <- node_lookup[[parent_node]]
            child_idx <- node_lookup[child_nodes]
            data.frame(
                x = nodes$x[[parent_idx]],
                xend = nodes$x[[parent_idx]],
                y = min(nodes$y[child_idx]),
                yend = max(nodes$y[child_idx]),
                stringsAsFactors = FALSE
            )
        }))
    }
    if (is.null(vertical_segments)) {
        vertical_segments <- data.frame()
    }

    if (length(lineages)) {
        lineage_paths <- do.call(rbind, lapply(seq_along(lineages), function(i) {
            lineage <- as.character(lineages[[i]])
            lineage <- lineage[lineage %in% names(node_lookup)]
            if (length(lineage) < 2) {
                return(NULL)
            }
            idx <- node_lookup[lineage]
            data.frame(
                cluster = lineage,
                x = nodes$x[idx],
                y = nodes$y[idx],
                lineage = paste("Lineage", i),
                order = seq_along(lineage),
                stringsAsFactors = FALSE
            )
        }))
    }
    if (is.null(lineage_paths)) {
        lineage_paths <- data.frame()
    }

    list(
        nodes = nodes,
        segments = segments,
        vertical_segments = vertical_segments,
        lineage_paths = lineage_paths,
        edges = edges,
        lineage_edges = lineage_edges
    )
}


plot_slingshot_mst <- function(res) {
    layout <- slingshot_mst_layout(res)
    nodes <- layout$nodes
    segments <- layout$segments
    vertical_segments <- layout$vertical_segments

    if (!nrow(nodes)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No Slingshot lineage graph available")
        return(invisible(NULL))
    }

    p <- ggplot2::ggplot()
    if (nrow(segments)) {
        p <- p + ggplot2::geom_segment(
            data = segments,
            ggplot2::aes(
                x = .data$x,
                xend = .data$xend,
                y = .data$y,
                yend = .data$yend
            ),
            linewidth = 2.8,
            color = "#c7c1b6",
            lineend = "round"
        )
    }
    if (nrow(vertical_segments)) {
        p <- p + ggplot2::geom_segment(
            data = vertical_segments,
            ggplot2::aes(
                x = .data$x,
                xend = .data$xend,
                y = .data$y,
                yend = .data$yend
            ),
            linewidth = 2.8,
            color = "#c7c1b6",
            lineend = "round"
        )
    }

    p +
        ggplot2::geom_point(
            data = nodes,
            ggplot2::aes(
                x = .data$x,
                y = .data$y,
                fill = .data$role
            ),
            shape = 21,
            color = "#2d2a23",
            size = 8,
            stroke = 0.75
        ) +
        ggplot2::geom_text(
            data = nodes,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$cluster),
            size = 3.2,
            fontface = "bold",
            color = "#202020"
        ) +
        ggplot2::scale_fill_manual(
            values = c(Cluster = "#f8f4e8", Root = "#b7e4c7", End = "#f7c59f"),
            breaks = c("Root", "End", "Cluster")
        ) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.14))) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
        ggplot2::labs(
            x = NULL,
            y = NULL,
            fill = NULL
        ) +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::theme(
            plot.background = ggplot2::element_rect(fill = "white", color = NA),
            panel.background = ggplot2::element_rect(fill = "white", color = NA),
            legend.background = ggplot2::element_rect(fill = "white", color = NA),
            legend.position = "none"
        )
}


#' Trajectory Inference Module - Server
#'
#' @param id Module namespace ID.
#' @export
#' @keywords internal
mod_trajectory_server <- function(id) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        setBookmarkExclude(c("rds_file", "custom_hm_rds_file", "custom_hm_genes_file"))

        ti_obj_available <- reactive({
            trajectory_rds_available(
                input$rds_source %||% "upload",
                input$rds_path %||% "",
                input$rds_file
            )
        })

        ti_obj <- reactive({
            if ((input$rds_source %||% "upload") == "upload") {
                req(input$rds_file)
                load_seurat_rds(input$rds_file$datapath)
            } else {
                rds_path <- resolve_input_path(input$rds_path %||% "")
                if (!file.exists(rds_path)) {
                    stop("Processed Seurat RDS does not exist: ", rds_path,
                         call. = FALSE)
                }
                load_seurat_rds(rds_path)
            }
        })

        output$start_cluster_ui <- renderUI({
            method <- input$trajectory_method %||% "slingshot"
            if (!isTRUE(ti_obj_available())) {
                return(cluster_control_placeholder())
            }

            obj <- ti_obj()
            clusters <- cluster_levels(obj)
            selectInput(ns("start_cluster"), trajectory_start_label(method),
                        choices = c("Auto" = "", clusters))
        })

        output$end_cluster_ui <- renderUI({
            if (!isTRUE(ti_obj_available())) {
                return(cluster_control_placeholder(
                    "Load a processed Seurat RDS to choose an optional end cluster."
                ))
            }

            obj <- ti_obj()
            clusters <- cluster_levels(obj)
            selectInput(ns("end_cluster"), "End cluster (optional):",
                        choices = c("Auto" = "", clusters))
        })

        output$paga_edge_ui <- renderUI({
            req(ti_result())
            res <- ti_result()
            if (!identical(res$method, "paga")) {
                return(NULL)
            }

            choices <- paga_edge_choices(res)
            selectInput(
                ns("paga_edge"),
                "Connected PAGA edge:",
                choices = choices,
                selected = unname(choices[[1]])
            )
        })

        ti_result <- eventReactive(input$run_ti, {
            if (!isTRUE(ti_obj_available())) {
                message <- if (identical(input$rds_source %||% "upload", "path")) {
                    paste0(
                        "Processed Seurat RDS does not exist: ",
                        resolve_input_path(input$rds_path %||% "")
                    )
                } else {
                    "Upload a processed Seurat RDS before running trajectory inference."
                }
                showNotification(message, type = "warning", duration = 10)
                return(NULL)
            }

            withProgress(message = "Running trajectory inference...", value = 0.2, {
                tryCatch({
                    obj <- ti_obj()
                    method <- input$trajectory_method %||% "slingshot"
                    start_cluster <- if (isTRUE(as.numeric(input$set_start %||% 1) == 1)) {
                        input$start_cluster
                    } else {
                        NULL
                    }
                    end_cluster <- if (
                        identical(method, "slingshot") &&
                            isTRUE(as.numeric(input$set_start %||% 1) == 1)
                    ) {
                        input$end_cluster
                    } else {
                        NULL
                    }

                    setProgress(0.4, message = paste("Running", method, "..."))
                    res <- switch(
                        method,
                        slingshot = run_slingshot_trajectory(obj, start_cluster, end_cluster),
                        paga = run_paga_trajectory(obj, start_cluster, end_cluster),
                        monocle3 = run_monocle3_trajectory(obj, start_cluster, end_cluster),
                        run_slingshot_trajectory(obj, start_cluster, end_cluster)
                    )

                    setProgress(1, message = paste(res$method_label, "complete!"))
                    res
                }, error = function(e) {
                    showNotification(
                        paste0("Trajectory inference failed: ", conditionMessage(e)),
                        type = "error",
                        duration = 15
                    )
                    NULL
                })
            })
        })

        pseudotime_plot <- reactive({
            req(ti_result())
            res <- ti_result()

            scCustomize::FeaturePlot_scCustom(
                seurat_object = res$obj,
                features = "pseudotime",
                order = TRUE
            ) +
                ggplot2::labs(title = NULL) +
                ggplot2::theme_void() +
                ggplot2::theme(plot.title = ggplot2::element_blank())
        })

        trajectory_plot <- reactive({
            req(ti_result())
            res <- ti_result()

            function() {
                if (identical(res$method, "slingshot")) {
                    print(plot_slingshot_mst(res))
                } else if (identical(res$method, "paga")) {
                    plot_paga_graph(res)
                } else if (identical(res$method, "monocle3")) {
                    print(monocle3::plot_cells(
                        res$cds,
                        color_cells_by = "pseudotime",
                        label_groups_by_cluster = FALSE,
                        label_leaves = FALSE,
                        label_branch_points = FALSE
                    ) +
                        ggplot2::labs(title = NULL) +
                        ggplot2::theme(plot.title = ggplot2::element_blank()))
                }
            }
        })

        output$pseudotime_plot <- renderPlot({
            req(ti_result())
            pseudotime_plot()
        })

        output$trajectory_plot <- renderPlot({
            req(ti_result())
            trajectory_plot()()
        })

        output$lineage_info <- renderPrint({
            req(ti_result())
            res <- ti_result()

            cat("Method:", res$method_label, "\n")
            if (!is.null(res$root_cluster)) {
                if (identical(res$method, "slingshot")) {
                    cat("Start cluster:", res$root_cluster, "\n")
                } else {
                    cat("Root/start cluster:", res$root_cluster, "\n")
                }
            }
            if (!is.null(res$end_cluster)) {
                cat("End cluster:", res$end_cluster, "\n")
            }

            if (identical(res$method, "slingshot")) {
                lineages <- slingshot::slingLineages(res$sds)
                lineage_counts <- table(res$obj$trajectory_path, useNA = "ifany")

                cat("\nLineages:\n")
                for (i in seq_along(lineages)) {
                    cat("Lineage", i, ":", paste(lineages[[i]], collapse = " -> "), "\n")
                }
                cat("\nCells assigned per lineage:\n")
                print(lineage_counts)
            } else if (identical(res$method, "paga")) {
                cat("Root cell:", res$root_cell, "\n\n")
                cat("Interpretation:\n")
                cat("PAGA estimates cluster connectivity/topology, not explicit lineages.\n")
                cat("Pseudotime is calculated from the selected root/start cluster.\n")
                cat("Use connected clusters or pseudotime-ordered cells for downstream gene discovery.\n\n")
                cat("Cells per PAGA cluster:\n")
                print(table(res$obj$paga_cluster, useNA = "ifany"))
                cat("\nStrongest PAGA cluster connections:\n")
                connections <- strongest_paga_connections(res$paga$connectivities)
                if (nrow(connections)) {
                    connections$connectivity <- round(connections$connectivity, 3)
                    print(connections, row.names = FALSE)
                } else {
                    cat("No cluster connections above the display threshold.\n")
                }
            } else if (identical(res$method, "monocle3")) {
                cat("Root cells:", length(res$root_cells), "\n\n")
                cat("Cells per Monocle 3 partition:\n")
                print(table(res$obj$monocle3_partition, useNA = "ifany"))
            }
        })

        output$download_ti_plot <- downloadHandler(
            filename = function() {
                fmt <- input$ti_format %||% "png"
                plot_name <- if (identical(input$ti_plot_choice %||% "pseudotime", "trajectory")) {
                    "trajectory_graph"
                } else {
                    "trajectory_pseudotime"
                }
                paste0(plot_name, "_", Sys.Date(), ".", fmt)
            },
            content = function(file) {
                req(ti_result())
                fmt <- input$ti_format %||% "png"
                w <- input$ti_width %||% 15
                h <- input$ti_height %||% 10
                res <- as.numeric(input$ti_res %||% 300)
                selected_plot <- if (identical(input$ti_plot_choice %||% "pseudotime", "trajectory")) {
                    trajectory_plot()
                } else {
                    pseudotime_plot()
                }

                write_plot_file(
                    plot = selected_plot,
                    file = file,
                    format = fmt,
                    width = w,
                    height = h,
                    dpi = res
                )
            }
        )

        output$download_ti_rds <- downloadHandler(
            filename = function() paste0("trajectory_", Sys.Date(), ".rds"),
            content = function(file) {
                req(ti_result())
                saveRDS(ti_result()$obj, file)
            }
        )

        custom_traj_heatmap_request <- eventReactive(input$custom_hm_run, {
            tryCatch({
                if (is.null(input$custom_hm_genes_file)) {
                    stop("Upload a gene list before generating the trajectory heatmap.",
                         call. = FALSE)
                }

                obj <- if (identical(input$custom_hm_source %||% "current", "current")) {
                    current <- tryCatch(ti_result(), error = function(e) NULL)
                    if (is.null(current)) {
                        stop(
                            "Run trajectory inference first, or choose Load trajectory RDS.",
                            call. = FALSE
                        )
                    }
                    current$obj
                } else if (identical(input$custom_hm_rds_source %||% "upload", "upload")) {
                    if (is.null(input$custom_hm_rds_file) ||
                        !file.exists(input$custom_hm_rds_file$datapath)) {
                        stop("Upload a trajectory RDS before generating the heatmap.",
                             call. = FALSE)
                    }
                    load_seurat_rds(input$custom_hm_rds_file$datapath)
                } else {
                    rds_path <- resolve_input_path(input$custom_hm_rds_path %||% "")
                    if (!file.exists(rds_path)) {
                        stop("Trajectory RDS does not exist: ", rds_path, call. = FALSE)
                    }
                    load_seurat_rds(rds_path)
                }

                validate_trajectory_heatmap_object(obj)

                has_header <- identical(as.character(input$custom_hm_gene_header %||% "1"), "1")
                genes <- read_trajectory_gene_list(
                    input$custom_hm_genes_file$datapath,
                    header = has_header
                )
                if (!length(genes)) {
                    stop("The uploaded gene list did not contain any gene names.",
                         call. = FALSE)
                }

                max_genes <- suppressWarnings(as.integer(input$custom_hm_max_genes %||% 20))
                if (!is.finite(max_genes) || max_genes < 1L) {
                    max_genes <- length(genes)
                }
                genes <- utils::head(genes, max_genes)
                present_genes <- intersect(genes, rownames(obj))
                missing_genes <- setdiff(genes, present_genes)
                if (!length(present_genes)) {
                    stop("None of the uploaded genes are present in the trajectory object.",
                         call. = FALSE)
                }
                if (length(missing_genes)) {
                    showNotification(
                        paste0(
                            length(missing_genes),
                            " gene(s) from the list were not found in the trajectory object."
                        ),
                        type = "warning",
                        duration = 10
                    )
                }

                list(
                    obj = obj,
                    genes = present_genes,
                    requested_genes = length(genes),
                    missing_genes = length(missing_genes)
                )
            }, error = function(e) {
                showNotification(
                    paste0("Custom trajectory heatmap failed: ", conditionMessage(e)),
                    type = "error",
                    duration = 15
                )
                NULL
            })
        }, ignoreInit = TRUE)

        custom_traj_heatmap_plot <- reactive({
            request <- custom_traj_heatmap_request()
            req(request)
            trajectory_heatmap_plot(request$obj, request$genes)
        })

        output$custom_traj_heatmap <- renderPlot({
            custom_traj_heatmap_plot()
        })

        output$custom_traj_heatmap_info <- renderPrint({
            request <- custom_traj_heatmap_request()
            req(request)
            metadata <- request$obj[[]]
            paths <- sort_cluster_labels(unique(as.character(metadata$trajectory_path)))
            paths <- paths[!is.na(paths) & nzchar(paths)]

            cat("Genes requested:", request$requested_genes, "\n")
            cat("Genes plotted:", length(request$genes), "\n")
            if (request$missing_genes > 0) {
                cat("Genes not found:", request$missing_genes, "\n")
            }
            cat("Cells:", ncol(request$obj), "\n")
            cat("Trajectory paths:", paste(paths, collapse = ", "), "\n")
        })

        output$download_custom_traj_heatmap <- plot_download_handler(
            plot_reactive = reactive(custom_traj_heatmap_plot()),
            filename_base = "custom_trajectory_heatmap",
            input = input,
            prefix = "custom_traj_hm"
        )

        paga_gene_result <- eventReactive(input$run_paga_genes, {
            req(ti_result())
            res <- ti_result()
            if (!identical(res$method, "paga")) {
                showNotification(
                    "Run PAGA trajectory inference before finding PAGA genes.",
                    type = "warning",
                    duration = 10
                )
                return(NULL)
            }

            if (identical(input$paga_genes_to_test %||% "variable", "all")) {
                showNotification(
                    "Testing all genes can take longer. Variable genes are recommended for teaching-sized runs.",
                    type = "message",
                    duration = 8
                )
            }

            withProgress(message = "Finding PAGA genes...", value = 0.2, {
                setProgress(0.35, message = "Preparing connected cluster edges...")
                edges <- selected_paga_edges(
                    res,
                    scope = input$paga_edge_scope %||% "selected",
                    selected_edge = input$paga_edge,
                    top_n = input$paga_top_edges %||% 5
                )
                gene_subset <- trajectory_gene_subset(
                    res$obj,
                    mode = input$paga_genes_to_test %||% "variable"
                )

                result <- tryCatch({
                    if (identical(input$paga_gene_mode %||% "edge_markers", "pseudotime")) {
                        setProgress(0.65, message = "Ranking genes by pseudotime association...")
                        run_paga_pseudotime_genes(res, edges, gene_subset = gene_subset)
                    } else {
                        setProgress(0.65, message = "Comparing connected clusters...")
                        run_paga_edge_markers(res, edges, gene_subset = gene_subset)
                    }
                }, error = function(e) {
                    showNotification(
                        paste0("PAGA gene discovery failed: ", conditionMessage(e)),
                        type = "error",
                        duration = 15
                    )
                    NULL
                })

                if (!is.null(result)) {
                    setProgress(1, message = paste0(nrow(result$table), " PAGA gene rows found"))
                }
                result
            })
        })

        output$paga_genes_table <- renderReactable({
            req(paga_gene_result())
            df <- utils::head(paga_gene_result()$table, input$paga_top_genes)
            if (nrow(df) == 0) {
                return(reactable(data.frame(Message = "No PAGA genes found with the current settings.")))
            }

            num_cols <- sapply(df, is.numeric)
            df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))

            columns <- list()
            if ("edge" %in% colnames(df)) {
                columns$edge <- reactable::colDef(name = "Edge", minWidth = 95)
            }
            if ("source_cluster" %in% colnames(df)) {
                columns$source_cluster <- reactable::colDef(name = "Source", width = 85)
            }
            if ("target_cluster" %in% colnames(df)) {
                columns$target_cluster <- reactable::colDef(name = "Target", width = 85)
            }
            if ("gene" %in% colnames(df)) {
                columns$gene <- reactable::colDef(name = "Gene", minWidth = 110)
            }
            if ("direction" %in% colnames(df)) {
                columns$direction <- reactable::colDef(
                    name = "Direction",
                    minWidth = 150,
                    style = list(whiteSpace = "normal", wordBreak = "break-word")
                )
            }
            if ("avg_log2FC" %in% colnames(df)) {
                columns$avg_log2FC <- reactable::colDef(name = "Avg log2FC", width = 100)
            }
            if ("rho" %in% colnames(df)) {
                columns$rho <- reactable::colDef(name = "Spearman rho", width = 110)
            }
            if ("p_value" %in% colnames(df)) {
                columns$p_value <- reactable::colDef(name = "P value", width = 90)
            }
            if ("p_val_adj" %in% colnames(df)) {
                columns$p_val_adj <- reactable::colDef(name = "Adj. P value", width = 105)
            }
            if ("pct_source" %in% colnames(df)) {
                columns$pct_source <- reactable::colDef(name = "Pct source", width = 100)
            }
            if ("pct_target" %in% colnames(df)) {
                columns$pct_target <- reactable::colDef(name = "Pct target", width = 100)
            }
            if ("connectivity" %in% colnames(df)) {
                columns$connectivity <- reactable::colDef(name = "Connectivity", width = 110)
            }

            reactable(
                df,
                striped = TRUE,
                highlight = TRUE,
                defaultPageSize = 25,
                defaultColDef = reactable::colDef(minWidth = 90),
                columns = columns
            )
        })

        output$paga_dot_plot <- renderPlot({
            req(paga_gene_result(), ti_result())
            result <- paga_gene_result()
            obj <- prepare_trajectory_de_object(ti_result()$obj)
            obj$paga_cluster <- as.character(ti_result()$obj$paga_cluster)
            genes <- utils::head(unique(result$table$gene), input$paga_top_genes)
            genes <- intersect(genes, rownames(obj))
            req(length(genes) > 0)

            clusters <- stats::setNames(as.character(obj$paga_cluster), colnames(obj))
            edge_clusters <- unique(c(result$edges$source_cluster, result$edges$target_cluster))
            cells <- names(clusters)[clusters %in% edge_clusters]
            req(length(cells) > 1)
            plot_obj <- subset(obj, cells = cells)
            plot_obj$paga_cluster <- factor(
                as.character(plot_obj$paga_cluster),
                levels = sort_cluster_labels(unique(as.character(plot_obj$paga_cluster)))
            )

            Seurat::DotPlot(
                plot_obj,
                features = genes,
                group.by = "paga_cluster"
            ) +
                ggplot2::coord_flip() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                    axis.text.y = ggplot2::element_text(angle = 0, hjust = 1),
                    plot.title = ggplot2::element_blank()
                ) +
                ggplot2::labs(title = NULL)
        })

        output$paga_feature_plot <- renderPlot({
            req(paga_gene_result(), ti_result())
            trajectory_feature_plot(
                ti_result()$obj,
                paga_gene_result()$table$gene,
                max_genes = min(input$paga_top_genes %||% 10, 12)
            )
        })

        output$download_paga_genes <- downloadHandler(
            filename = function() paste0("paga_genes_", Sys.Date(), ".csv"),
            content = function(file) {
                write.csv(paga_gene_result()$table, file, row.names = FALSE)
            }
        )

        monocle_gene_result <- eventReactive(input$run_monocle_genes, {
            req(ti_result())
            res <- ti_result()
            if (!identical(res$method, "monocle3")) {
                showNotification(
                    "Run Monocle 3 trajectory inference before finding trajectory genes.",
                    type = "warning",
                    duration = 10
                )
                return(NULL)
            }

            if (identical(input$monocle_genes_to_test %||% "variable", "all")) {
                showNotification(
                    "Testing all genes can take longer. Variable genes are recommended for teaching-sized runs.",
                    type = "message",
                    duration = 8
                )
            }

            withProgress(message = "Finding Monocle 3 trajectory genes...", value = 0.3, {
                gene_subset <- trajectory_gene_subset(
                    res$obj,
                    mode = input$monocle_genes_to_test %||% "variable"
                )
                result <- tryCatch({
                    setProgress(0.65, message = "Running graph-based gene test...")
                    run_monocle3_graph_test(res, gene_subset = gene_subset, cores = 2)
                }, error = function(e) {
                    showNotification(
                        paste0("Monocle 3 gene discovery failed: ", conditionMessage(e)),
                        type = "error",
                        duration = 15
                    )
                    NULL
                })
                if (!is.null(result)) {
                    setProgress(1, message = paste0(nrow(result$table), " trajectory genes tested"))
                }
                result
            })
        })

        monocle_module_result <- eventReactive(input$run_monocle_modules, {
            req(ti_result(), monocle_gene_result())
            res <- ti_result()
            genes <- monocle_gene_result()$table
            genes <- genes[is.finite(genes$q_value) & genes$q_value <= 0.05, , drop = FALSE]
            if (nrow(genes) < 2L) {
                genes <- utils::head(monocle_gene_result()$table, 100)
            }

            withProgress(message = "Finding Monocle 3 gene modules...", value = 0.4, {
                tryCatch(
                    run_monocle3_gene_modules(res, genes$gene, cores = 2),
                    error = function(e) {
                        showNotification(
                            paste0("Monocle 3 gene modules failed: ", conditionMessage(e)),
                            type = "error",
                            duration = 15
                        )
                        NULL
                    }
                )
            })
        })

        output$monocle_genes_table <- renderReactable({
            req(monocle_gene_result())
            df <- utils::head(monocle_gene_result()$table, input$monocle_top_genes)
            if (nrow(df) == 0) {
                return(reactable(data.frame(Message = "No Monocle 3 trajectory genes found.")))
            }

            num_cols <- sapply(df, is.numeric)
            df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))
            reactable(df, striped = TRUE, highlight = TRUE, defaultPageSize = 25)
        })

        output$monocle_feature_plot <- renderPlot({
            req(monocle_gene_result(), ti_result())
            trajectory_feature_plot(
                ti_result()$obj,
                monocle_gene_result()$table$gene,
                max_genes = min(input$monocle_top_genes %||% 50, 12)
            )
        })

        output$monocle_gene_heatmap <- renderPlot({
            req(monocle_gene_result(), ti_result())
            genes <- utils::head(
                unique(monocle_gene_result()$table$gene),
                input$monocle_top_genes %||% 50
            )
            trajectory_heatmap_plot(ti_result()$obj, genes)
        })

        output$monocle_modules_table <- renderReactable({
            req(monocle_module_result())
            result <- monocle_module_result()
            if (is.null(result) || !nrow(result$summary)) {
                return(reactable(data.frame(Message = "No Monocle 3 gene modules found.")))
            }

            reactable(result$summary, striped = TRUE, highlight = TRUE,
                      defaultPageSize = 25)
        })

        output$monocle_modules_heatmap <- renderPlot({
            req(ti_result(), monocle_module_result())
            plot_monocle_module_heatmap(ti_result(), monocle_module_result())
        })

        output$download_monocle_genes <- downloadHandler(
            filename = function() paste0("monocle3_trajectory_genes_", Sys.Date(), ".csv"),
            content = function(file) {
                write.csv(monocle_gene_result()$table, file, row.names = FALSE)
            }
        )

        output$download_monocle_modules <- downloadHandler(
            filename = function() paste0("monocle3_module_genes_", Sys.Date(), ".csv"),
            content = function(file) {
                req(monocle_module_result())
                modules <- monocle_module_result()$modules
                gene_col <- intersect(
                    c("gene_id", "id", "gene_short_name", "gene"),
                    colnames(modules)
                )
                if (length(gene_col)) {
                    modules$gene <- as.character(modules[[gene_col[[1]]]])
                } else {
                    modules$gene <- rownames(modules)
                }
                if (!"module" %in% colnames(modules)) {
                    modules$module <- "Module 1"
                }
                modules <- modules[, unique(c("gene", "module", colnames(modules))),
                                   drop = FALSE]
                write.csv(modules, file, row.names = FALSE)
            }
        )

        traj_de_result <- eventReactive(input$run_traj_de, {
            req(ti_result())
            res <- ti_result()

            if (!identical(res$method, "slingshot")) {
                showNotification(
                    "Trajectory DE is currently available for Slingshot results only.",
                    type = "warning",
                    duration = 10
                )
                return(NULL)
            }

            engine <- input$traj_de_method %||% "pseudotimede"
            if (identical(engine, "pseudotimede") && !pseudotimede_available()) {
                showNotification(
                    pseudotimede_install_message(),
                    type = "warning",
                    duration = 15
                )
                return(NULL)
            }

            progress_message <- switch(
                engine,
                pseudotimede = "Running PseudotimeDE-fast...",
                tradeseq = "Running tradeSeq. This can take a long time on large datasets...",
                "Running scMaSigPro..."
            )

            withProgress(message = progress_message, value = 0.3, {
                setProgress(0.5, message = "Preparing variable genes and lineage metadata...")
                variable_features <- Seurat::VariableFeatures(res$obj)
                if (!length(variable_features)) {
                    variable_features <- rownames(res$obj)
                }
                result <- tryCatch(
                    if (identical(engine, "tradeseq")) {
                        run_tradeseq(
                            res$sce,
                            res$assignment,
                            gene_subset = variable_features,
                            nknots = input$tradeseq_nknots %||% 6
                        )
                    } else if (identical(engine, "pseudotimede")) {
                        run_pseudotimede_fast(
                            res$sce,
                            gene_subset = variable_features,
                            cores = 2
                        )
                    } else {
                        run_scmasigpro(
                            res$sce,
                            gene_subset = variable_features
                        )
                    },
                    error = function(e) {
                        showNotification(
                            paste0("Trajectory DE failed: ", conditionMessage(e)),
                            type = "error",
                            duration = 15
                        )
                        NULL
                    }
                )
                if (!is.null(result)) {
                    if (!is.null(result$message)) {
                        showNotification(result$message, type = "message", duration = 10)
                    }
                    setProgress(1, message = paste0(nrow(result$table), " trajectory genes found"))
                }
                result
            })
        })

        output$traj_genes_table <- renderReactable({
            req(traj_de_result())
            df <- utils::head(traj_de_result()$table, input$n_top_genes)
            if (nrow(df) == 0) {
                return(reactable(data.frame(Message = "No significant trajectory genes found")))
            }

            num_cols <- sapply(df, is.numeric)
            df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))
            if ("lineages" %in% colnames(df)) {
                df$lineages <- gsub("vs", " vs ", df$lineages, fixed = TRUE)
                df$lineages <- gsub(",\\s*", ", ", df$lineages)
                df$lineages <- gsub("\\s+", " ", df$lineages)
            }

            columns <- list()
            if ("gene" %in% colnames(df)) {
                columns$gene <- reactable::colDef(name = "Gene", minWidth = 95)
            }
            if ("p_value" %in% colnames(df)) {
                columns$p_value <- reactable::colDef(name = "P value", width = 90)
            }
            if ("rsquared" %in% colnames(df)) {
                columns$rsquared <- reactable::colDef(name = "R squared", width = 95)
            }
            if ("n_lineages" %in% colnames(df)) {
                columns$n_lineages <- reactable::colDef(name = "Lineages", width = 80)
            }
            if ("lineages" %in% colnames(df)) {
                columns$lineages <- reactable::colDef(
                    name = "Lineage terms",
                    minWidth = 280,
                    style = list(whiteSpace = "normal", wordBreak = "break-word"),
                    headerStyle = list(whiteSpace = "normal")
                )
            }

            reactable(
                df,
                striped = TRUE,
                highlight = TRUE,
                defaultPageSize = 25,
                defaultColDef = reactable::colDef(minWidth = 90),
                columns = columns
            )
        })

        output$traj_heatmap <- renderPlot({
            req(traj_de_result(), ti_result())
            obj <- ti_result()$obj
            top_genes <- utils::head(traj_de_result()$table$gene, input$n_top_genes)
            top_genes <- intersect(top_genes, rownames(obj))
            req(length(top_genes) > 0)

            trajectory_heatmap_plot(obj, top_genes)
        })

        output$traj_feature_plot <- renderPlot({
            req(traj_de_result(), ti_result())
            trajectory_feature_plot(
                ti_result()$obj,
                traj_de_result()$table$gene,
                max_genes = min(input$n_top_genes %||% 10, 12)
            )
        })

        output$download_traj_genes <- downloadHandler(
            filename = function() paste0("trajectory_genes_", Sys.Date(), ".csv"),
            content = function(file) {
                write.csv(traj_de_result()$table, file, row.names = FALSE)
            }
        )

        list(
            ti_result = ti_result,
            traj_de_result = traj_de_result,
            paga_gene_result = paga_gene_result,
            monocle_gene_result = monocle_gene_result,
            monocle_module_result = monocle_module_result
        )
    })
}
