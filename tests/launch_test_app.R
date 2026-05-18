# ============================================================================
# Asc-Seurat v3 — Standalone Launcher
# ============================================================================

cat("Asc-Seurat v3 launcher\n")
cat("Loading dependencies...\n")

suppressPackageStartupMessages({
    library(shiny)
    library(bslib)
    library(Seurat)
    library(SeuratObject)
    library(scCustomize)
    library(slingshot)
    library(harmony)
    library(ComplexHeatmap)
    library(SingleCellExperiment)
    library(ggplot2)
    library(dplyr)
    library(patchwork)
    library(reactable)
    library(DT)
    library(shinyWidgets)
    library(shinyFeedback)
    library(shinycssloaders)
    library(sctransform)
    library(hdf5r)
    library(future)
    library(metap)
    library(circlize)
    library(reticulate)
    library(sass)
})

enableBookmarking("server")

if (interactive()) {
    base_dir <- rstudioapi::getActiveProject()
    if (is.null(base_dir)) {
        base_dir <- getwd()
    }
} else {
    args <- commandArgs(trailingOnly = FALSE)
    script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
    base_dir <- normalizePath(file.path(dirname(script_path), ".."))
}

cat("Project root:", base_dir, "\n")

r_dir <- file.path(base_dir, "R")
source_files <- c(
    "ascseurat-package.R",
    "utils_compat.R",
    "utils_doublets.R",
    "utils_seurat.R",
    "utils_plots.R",
    "utils_ui.R",
    "mod_data_loading.R",
    "mod_qc.R",
    "mod_normalization.R",
    "mod_clustering.R",
    "mod_de.R",
    "mod_visualization.R",
    "mod_integration.R",
    "mod_trajectory.R",
    "mod_celltype_annotation.R",
    "mod_advanced_plots.R",
    "run_app.R"
)
for (f in source_files) {
    source(file.path(r_dir, f), local = FALSE)
}

dir.create(file.path(base_dir, "data"), showWarnings = FALSE)
dir.create(file.path(base_dir, "RDS_files"), showWarnings = FALSE)

css_path <- file.path(base_dir, "inst", "app", "www", "custom.css")

asc_theme <- bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#6366f1",
    secondary = "#64748b",
    success = "#10b981",
    info = "#06b6d4",
    warning = "#f59e0b",
    danger = "#ef4444",
    "font-size-base" = "0.95rem",
    "enable-rounded" = TRUE,
    "enable-shadows" = TRUE
)
if (file.exists(css_path)) {
    asc_theme <- asc_theme |> bs_add_rules(sass::sass_file(css_path))
}

ui <- function(request) {
    tagList(
        shinyFeedback::useShinyFeedback(),
        page_navbar(
            title = tags$span(
                icon("project-diagram", class = "me-2"),
                "Asc-Seurat v3"
            ),
            id = "main_nav",
            theme = asc_theme,
            fillable = FALSE,
            header = tags$head(
                tags$meta(charset = "UTF-8"),
                tags$meta(
                    name = "description",
                    content = "Asc-Seurat: Interactive single-cell RNA-seq analysis"
                ),
                tags$meta(
                    name = "viewport",
                    content = "width=device-width, initial-scale=1"
                )
            ),
            nav_spacer(),
            nav_item(actionButton(
                "session_report",
                "Session report",
                icon = icon("file-lines"),
                class = "btn btn-sm btn-outline-secondary"
            )),
            nav_item(bookmarkButton(
                label = NULL,
                icon = icon("bookmark"),
                title = "Bookmark this analysis state",
                class = "btn btn-sm btn-outline-primary"
            )),
            nav_item(input_dark_mode(id = "dark_mode", mode = "light")),
            nav_panel(
                title = "Home",
                icon = icon("house"),
                value = "home",
                div(
                    class = "container-fluid py-4",
                    card(
                        card_header(
                            tags$h2("Welcome to Asc-Seurat v3"),
                            class = "bg-primary text-white text-center"
                        ),
                        card_body(
                            class = "text-center py-4",
                            tags$p(
                                class = "lead",
                                "Analytical Single Cell Seurat-Based Web Application"
                            ),
                            tags$p(
                                "An interactive platform for single-cell RNA-seq analysis, ",
                                "powered by Seurat v5."
                            ),
                            tags$hr(),
                            layout_columns(
                                col_widths = c(4, 4, 4),
                                value_box(
                                    title = "Single Sample",
                                    value = "Analysis",
                                    showcase = icon("braille"),
                                    theme = "primary",
                                    p("QC, normalization, clustering, and DE analysis.")
                                ),
                                value_box(
                                    title = "Multi-Sample",
                                    value = "Integration",
                                    showcase = icon("project-diagram"),
                                    theme = "info",
                                    p("Integrate multiple samples using RPCA or Harmony.")
                                ),
                                value_box(
                                    title = "Trajectory",
                                    value = "Inference",
                                    showcase = icon("chart-line"),
                                    theme = "success",
                                    p("Infer trajectories with Slingshot, PAGA, and Monocle 3.")
                                )
                            ),
                            tags$hr(),
                            tags$p(
                                class = "text-muted",
                                "Asc-Seurat v3.0.0 — ",
                                tags$a(
                                    href = "https://github.com/Pereira-Lab-UF/asc_seurat",
                                    "GitHub",
                                    target = "_blank"
                                ),
                                " | ",
                                tags$a(
                                    href = "https://doi.org/10.1186/s12859-021-04472-2",
                                    "Paper",
                                    target = "_blank"
                                )
                            )
                        )
                    )
                )
            ),
            nav_panel(
                title = "Single Sample",
                icon = icon("braille"),
                value = "single_sample",
                div(
                    class = "container-fluid py-3",
                    mod_data_loading_ui("single_data"),
                    mod_qc_ui("single_qc"),
                    mod_normalization_ui("single_norm"),
                    mod_clustering_ui("single_clust"),
                    mod_de_ui("single_de"),
                    mod_visualization_ui("single_viz")
                )
            ),
            nav_panel(
                title = "Integration",
                icon = icon("project-diagram"),
                value = "integration",
                div(
                    class = "container-fluid py-3",
                    mod_integration_ui("multi_int"),
                    mod_clustering_ui("multi_clust"),
                    mod_de_ui("multi_de", allow_conserved = TRUE),
                    mod_visualization_ui("multi_viz")
                )
            ),
            nav_panel(
                title = "Trajectory Inference",
                icon = icon("chart-line"),
                value = "trajectory",
                div(
                    class = "container-fluid py-3",
                    mod_trajectory_ui("traj")
                )
            ),
            nav_panel(
                title = "Cell-type Annotation",
                icon = icon("tags"),
                value = "celltype_annotation",
                div(
                    class = "container-fluid py-3",
                    mod_celltype_annotation_ui("celltype")
                )
            ),
            nav_panel(
                title = "Advanced Plots",
                icon = icon("chart-area"),
                value = "advanced_plots",
                div(
                    class = "container-fluid py-3",
                    mod_advanced_plots_ui("adv_plots")
                )
            )
        )
    )
}

server <- function(input, output, session) {
    onBookmarked(function(url) {
        showModal(modalDialog(
            title = "Bookmark created",
            tags$p("Use this URL to restore your analysis parameters."),
            tags$p("Uploaded files are not restored and must be selected again."),
            tags$input(
                type = "text",
                value = url,
                readonly = TRUE,
                class = "form-control"
            ),
            easyClose = TRUE
        ))
    })

    session_log <- reactiveValues(events = list())

    extract_scoped_inputs <- function(prefix) {
        values <- reactiveValuesToList(input)
        scoped_names <- grep(paste0("^", prefix), names(values), value = TRUE)
        scoped <- values[scoped_names]
        names(scoped) <- sub(paste0("^", prefix), "", names(scoped))

        lapply(scoped, function(value) {
            if (is.null(value)) {
                return(NULL)
            }
            if (is.list(value) && !is.null(value$name)) {
                return(paste(value$name, collapse = ", "))
            }
            if (length(value) > 10) {
                value <- value[seq_len(10)]
            }
            paste(as.character(value), collapse = ", ")
        })
    }

    append_log <- function(action, prefix) {
        session_log$events <- c(session_log$events, list(list(
            timestamp = Sys.time(),
            action = action,
            params = extract_scoped_inputs(prefix)
        )))
    }

    log_map <- list(
        list(id = "single_data-load_data", action = "Load single-sample data", prefix = "single_data-"),
        list(id = "single_data-load_demo", action = "Load demo data", prefix = "single_data-"),
        list(id = "single_data-apply_metadata", action = "Add single-sample metadata", prefix = "single_data-"),
        list(id = "single_qc-apply_filter", action = "QC filter (single sample)", prefix = "single_qc-"),
        list(id = "single_norm-run_norm", action = "Normalize single sample", prefix = "single_norm-"),
        list(id = "single_clust-run_clustering", action = "Cluster single sample", prefix = "single_clust-"),
        list(id = "single_clust-apply_rename", action = "Rename single-sample clusters", prefix = "single_clust-"),
        list(id = "single_clust-reanalyze", action = "Reanalyze single-sample subset", prefix = "single_clust-"),
        list(id = "single_de-run_de", action = "Run single-sample DE", prefix = "single_de-"),
        list(id = "multi_int-run_integration", action = "Run integration", prefix = "multi_int-"),
        list(id = "multi_clust-run_clustering", action = "Cluster integrated data", prefix = "multi_clust-"),
        list(id = "multi_clust-apply_rename", action = "Rename integrated clusters", prefix = "multi_clust-"),
        list(id = "multi_clust-reanalyze", action = "Reanalyze integrated subset", prefix = "multi_clust-"),
        list(id = "multi_de-run_de", action = "Run integrated DE", prefix = "multi_de-"),
        list(id = "traj-run_ti", action = "Run trajectory inference", prefix = "traj-"),
        list(id = "traj-run_traj_de", action = "Run trajectory DE", prefix = "traj-"),
        list(id = "celltype-run_singler", action = "Run SingleR annotation", prefix = "celltype-"),
        list(id = "adv_plots-make_plot", action = "Generate advanced plot", prefix = "adv_plots-")
    )
    for (entry in log_map) {
        local({
            current <- entry
            observeEvent(input[[current$id]], {
                append_log(current$action, current$prefix)
            }, ignoreInit = TRUE)
        })
    }

    build_session_report <- function() {
        lines <- c(
            "# Asc-Seurat session report",
            "",
            paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
            paste("ascseurat version:", as.character(packageVersion("ascseurat"))),
            paste("Seurat version:", as.character(packageVersion("Seurat"))),
            "",
            "## Steps run",
            ""
        )
        if (!length(session_log$events)) {
            return(c(lines, "No actions have been logged yet."))
        }

        for (i in seq_along(session_log$events)) {
            event <- session_log$events[[i]]
            lines <- c(
                lines,
                paste0(i, ". **", event$action, "** (", format(event$timestamp, "%H:%M:%S"), ")")
            )
            if (length(event$params)) {
                lines <- c(lines, vapply(names(event$params), function(name) {
                    paste0("   - ", name, ": ", event$params[[name]])
                }, character(1)))
            }
            lines <- c(lines, "")
        }
        lines
    }

    observeEvent(input$session_report, {
        showModal(modalDialog(
            title = "Session report",
            tags$p("Download a Markdown summary of the steps and parameters used in this session."),
            downloadButton("download_session_report", "Download Markdown",
                           class = "btn btn-outline-secondary"),
            easyClose = TRUE
        ))
    })

    output$download_session_report <- downloadHandler(
        filename = function() paste0("ascseurat_session_", Sys.Date(), ".md"),
        content = function(file) {
            writeLines(build_session_report(), file)
        }
    )

    single_data <- mod_data_loading_server("single_data")
    single_qc <- mod_qc_server("single_qc", seurat_obj = single_data)
    single_norm <- mod_normalization_server("single_norm", seurat_obj = single_qc)
    single_clust <- mod_clustering_server("single_clust", norm_result = single_norm)
    mod_de_server("single_de", seurat_obj = single_clust)
    mod_visualization_server("single_viz", seurat_obj = single_clust)

    multi_int <- mod_integration_server("multi_int")
    multi_clust <- mod_clustering_server("multi_clust", norm_result = multi_int)
    mod_de_server("multi_de", seurat_obj = multi_clust, allow_conserved = TRUE)
    mod_visualization_server("multi_viz", seurat_obj = multi_clust)

    mod_trajectory_server("traj")
    mod_celltype_annotation_server(
        "celltype",
        seurat_obj_single = single_clust,
        seurat_obj_multi = multi_clust
    )
    mod_advanced_plots_server(
        "adv_plots",
        seurat_obj_single = single_clust,
        seurat_obj_multi = multi_clust
    )
}

cat("\nLaunching Asc-Seurat v3 on port 7775...\n")
cat("Open: http://127.0.0.1:7775\n\n")

addResourcePath("www", file.path(base_dir, "inst", "app", "www"))

shinyApp(
    ui = ui,
    server = server,
    options = list(host = "127.0.0.1", port = 7775, launch.browser = FALSE)
)
