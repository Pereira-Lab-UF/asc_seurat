# Asc-Seurat v3 — Main Shiny Application

library(shiny)
library(bslib)

# In development (devtools::load_all), skip re-loading the installed package so
# source changes are picked up immediately without reinstalling.
if (!requireNamespace("pkgload", quietly = TRUE) || !pkgload::is_dev_package("ascseurat")) {
    library(ascseurat)
}

options(shiny.maxRequestSize = getOption("shiny.maxRequestSize", 5120 * 1024^2))

enableBookmarking("server")

asc_theme <- bs_theme(
    version = 5,
    # No bootswatch — we want to start from a clean Bootstrap base and apply
    # the "Editorial Lab Notebook" design system in custom.css.
    primary    = "#0c6e4a",   # deep emerald — single confident accent
    secondary  = "#52606d",   # warm slate
    success    = "#0c6e4a",
    info       = "#1f6f8b",
    warning    = "#a36b1c",
    danger     = "#a4341d",
    "body-bg"          = "#f7f5ee",   # parchment off-white
    "body-color"       = "#1c1b18",   # ink
    "border-color"     = "#e6e3d6",
    "font-family-base" = "'Inter Tight', system-ui, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif",
    "headings-font-family" = "'Fraunces', 'Iowan Old Style', Palatino, 'Times New Roman', serif",
    "headings-font-weight" = "600",
    "headings-color"   = "#1c1b18",
    "font-size-base"   = "0.95rem",
    "line-height-base" = "1.55",
    "border-radius"    = "0.4rem",
    "enable-rounded"   = TRUE,
    "enable-shadows"   = FALSE
) |>
    bs_add_rules(paste(readLines(
        system.file("app", "www", "custom.css", package = "ascseurat")
    ), collapse = "\n"))

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
                "Session Report",
                icon = icon("file-lines"),
                class = "btn btn-sm btn-light asc-session-report"
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

                    # ── Welcome banner ────────────────────────────────────────
                    card(
                        card_header(
                            class = "bg-primary text-white text-center py-3",
                            tags$h2(class = "mb-0", "Welcome to Asc-Seurat v3")
                        ),
                        card_body(
                            class = "text-center py-3",
                            tags$p(
                                class = "lead mb-0",
                                "An interactive platform for single-cell RNA-seq analysis, ",
                                "powered by Seurat v5. Designed for wet-lab biologists — ",
                                "no programming required."
                            )
                        )
                    ),

                    tags$hr(),

                    # ── 4-step Quick Start ────────────────────────────────────
                    tags$h5(
                        class = "mb-3",
                        icon("rocket", class = "me-1 text-primary"),
                        "How to analyze your data in 4 steps"
                    ),
                    layout_columns(
                        col_widths = c(3, 3, 3, 3),
                        card(
                            class = "h-100 border-primary",
                            card_header(
                                class = "bg-primary text-white",
                                tags$span(class = "badge bg-white text-primary me-2", "1"),
                                "Load Data"
                            ),
                            card_body(
                                tags$p(
                                    icon("folder-open", class = "text-primary me-1"),
                                    "Upload your data: 10X Genomics directory, HDF5 (.h5), ",
                                    "AnnData (.h5ad), CSV/TSV count matrix, or a Seurat RDS."
                                ),
                                tags$p(
                                    class = "text-muted small mb-0",
                                    icon("arrow-right", class = "me-1"),
                                    tags$em("Single Sample tab \u2192 Step 1")
                                )
                            )
                        ),
                        card(
                            class = "h-100 border-success",
                            card_header(
                                class = "bg-success text-white",
                                tags$span(class = "badge bg-white text-success me-2", "2"),
                                "QC & Filter"
                            ),
                            card_body(
                                tags$p(
                                    icon("filter", class = "text-success me-1"),
                                    "Remove low-quality cells. Set thresholds for gene counts ",
                                    "and mitochondrial content. Optionally detect doublets."
                                ),
                                tags$p(
                                    class = "text-muted small mb-0",
                                    icon("arrow-right", class = "me-1"),
                                    tags$em("Single Sample tab \u2192 Step 2")
                                )
                            )
                        ),
                        card(
                            class = "h-100 border-info",
                            card_header(
                                class = "bg-info text-white",
                                tags$span(class = "badge bg-white text-info me-2", "3"),
                                "Normalize & Cluster"
                            ),
                            card_body(
                                tags$p(
                                    icon("chart-bar", class = "text-info me-1"),
                                    "Normalize expression data, run PCA, cluster cells, ",
                                    "and generate UMAP visualizations."
                                ),
                                tags$p(
                                    class = "text-muted small mb-0",
                                    icon("arrow-right", class = "me-1"),
                                    tags$em("Single Sample tab \u2192 Steps 3 & 4")
                                )
                            )
                        ),
                        card(
                            class = "h-100 border-warning",
                            card_header(
                                class = "bg-warning",
                                tags$span(class = "badge bg-white text-warning me-2", "4"),
                                "Explore Results"
                            ),
                            card_body(
                                tags$p(
                                    icon("table", class = "text-warning me-1"),
                                    "Run differential expression, visualize marker genes, ",
                                    "and build publication-ready expression plots."
                                ),
                                tags$p(
                                    class = "text-muted small mb-0",
                                    icon("arrow-right", class = "me-1"),
                                    tags$em("Single Sample tab \u2192 Steps 5 & 6")
                                )
                            )
                        )
                    ),

                    tags$hr(),

                    # ── Which pipeline + supported formats ────────────────────
                    layout_columns(
                        col_widths = c(6, 6),
                        card(
                            card_header(
                                class = "bg-light",
                                tags$h5(
                                    class = "mb-0",
                                    icon("question-circle", class = "me-1"),
                                    "Which pipeline should I use?"
                                )
                            ),
                            card_body(
                                tags$dl(
                                    class = "row mb-0",
                                    tags$dt(
                                        class = "col-sm-4",
                                        icon("braille", class = "me-1 text-primary"),
                                        "Single Sample"
                                    ),
                                    tags$dd(
                                        class = "col-sm-8",
                                        "One sample or dataset. Start here for ",
                                        "QC \u2192 cluster \u2192 DE analysis."
                                    ),
                                    tags$dt(
                                        class = "col-sm-4",
                                        icon("project-diagram", class = "me-1 text-info"),
                                        "Integration"
                                    ),
                                    tags$dd(
                                        class = "col-sm-8",
                                        "Multiple samples (e.g., control + treatment). ",
                                        "Batch-corrected with Harmony or RPCA."
                                    ),
                                    tags$dt(
                                        class = "col-sm-4",
                                        icon("chart-line", class = "me-1 text-success"),
                                        "Trajectory"
                                    ),
                                    tags$dd(
                                        class = "col-sm-8",
                                        "Studying differentiation or a developmental ",
                                        "process. Pseudotime inference with Slingshot."
                                    ),
                                    tags$dt(
                                        class = "col-sm-4",
                                        icon("tags", class = "me-1 text-danger"),
                                        "Cell-type Annotation"
                                    ),
                                    tags$dd(
                                        class = "col-sm-8",
                                        "Automatically label clusters using a reference ",
                                        "dataset (SingleR)."
                                    )
                                )
                            )
                        ),
                        card(
                            card_header(
                                class = "bg-light",
                                tags$h5(
                                    class = "mb-0",
                                    icon("file", class = "me-1"),
                                    "Supported Input Formats"
                                )
                            ),
                            card_body(
                                tags$table(
                                    class = "table table-sm table-borderless mb-1",
                                    tags$tbody(
                                        tags$tr(
                                            tags$td(tags$code("directory")),
                                            tags$td("10X Genomics"),
                                            tags$td(
                                                class = "text-muted small",
                                                "barcodes / features / matrix files"
                                            )
                                        ),
                                        tags$tr(
                                            tags$td(tags$code(".h5")),
                                            tags$td("10X HDF5"),
                                            tags$td(
                                                class = "text-muted small",
                                                "Filtered feature-barcode matrix"
                                            )
                                        ),
                                        tags$tr(
                                            tags$td(tags$code(".h5ad")),
                                            tags$td("AnnData"),
                                            tags$td(
                                                class = "text-muted small",
                                                "Python / Scanpy output"
                                            )
                                        ),
                                        tags$tr(
                                            tags$td(tags$code(".csv / .tsv")),
                                            tags$td("Count matrix"),
                                            tags$td(
                                                class = "text-muted small",
                                                "Genes as rows, cells as columns"
                                            )
                                        ),
                                        tags$tr(
                                            tags$td(tags$code(".rds")),
                                            tags$td("Seurat RDS"),
                                            tags$td(
                                                class = "text-muted small",
                                                "Resume from a previously saved object"
                                            )
                                        )
                                    )
                                ),
                                tags$p(
                                    class = "text-muted small mb-0",
                                    icon("info-circle", class = "me-1"),
                                    "For 10X Genomics data, specify the folder containing ",
                                    "the three matrix files — not an individual file."
                                )
                            )
                        )
                    ),

                    div(
                        class = "home-cta-stack",
                        tags$a(
                            class = "btn btn-outline-secondary btn-lg",
                            href = "https://asc-seurat.readthedocs.io",
                            target = "_blank",
                            icon("book"), " Documentation"
                        ),
                        actionButton(
                            "go_to_demo",
                            label = "Try it with Demo Data",
                            icon = icon("play"),
                            class = "btn btn-success btn-lg"
                        )
                    ),

                    tags$hr(),

                    tags$p(
                        class = "text-center text-muted small",
                        "Asc-Seurat v3.0.0 \u2014 ",
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
            ),
            nav_panel(
                title = "Single Sample",
                icon = icon("braille"),
                value = "single_sample",
                div(
                    class = "container-fluid py-3",
                    tags$div(
                        class = "visually-hidden",
                        textOutput("single_data-data_ready"),
                        textOutput("single_norm-norm_ready")
                    ),
                    mod_data_loading_ui("single_data"),
                    conditionalPanel(
                        condition = "output['single_data-data_ready'] == 'true'",
                        mod_qc_ui("single_qc"),
                        mod_normalization_ui("single_norm")
                    ),
                    conditionalPanel(
                        condition = "output['single_norm-norm_ready'] == 'true'",
                        mod_clustering_ui("single_clust", step_number = 4)
                    ),
                    conditionalPanel(
                        condition = "output['single_clust-clustering_ready'] == 'true'",
                        mod_de_ui("single_de", step_number = 5),
                        mod_visualization_ui("single_viz", step_number = 6)
                    ),
                    uiOutput("single_clust-rename_section")
                )
            ),
            nav_panel(
                title = "Integration",
                icon = icon("project-diagram"),
                value = "integration",
                div(
                    class = "container-fluid py-3",
                    tags$div(
                        class = "visually-hidden",
                        textOutput("multi_int-integration_ready")
                    ),
                    mod_integration_ui("multi_int"),
                    conditionalPanel(
                        condition = "output['multi_int-integration_ready'] == 'true'",
                        mod_clustering_ui(
                            "multi_clust",
                            step_number = 2,
                            download_label = "Download Integrated and Clustered Object"
                        )
                    ),
                    conditionalPanel(
                        condition = "output['multi_clust-clustering_ready'] == 'true'",
                        mod_de_ui("multi_de", allow_conserved = TRUE, step_number = 3),
                        mod_visualization_ui("multi_viz", step_number = 4)
                    ),
                    uiOutput("multi_clust-rename_section")
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
                title = "Tools",
                icon = icon("toolbox"),
                value = "tools",
                div(
                    class = "container-fluid py-3",
                    bslib::navset_card_tab(
                        id = "tools_nav",
                        bslib::nav_panel(
                            title = tagList(icon("tags"), "Cell-type Annotation"),
                            value = "celltype_annotation",
                            mod_celltype_annotation_ui("celltype")
                        ),
                        bslib::nav_panel(
                            title = tagList(icon("chart-area"), "Advanced Plots"),
                            value = "advanced_plots",
                            mod_advanced_plots_ui("adv_plots")
                        )
                    )
                )
            ),
            nav_panel(
                title = "Demo",
                icon = icon("play-circle"),
                value = "demo",
                div(
                    class = "container-fluid py-3",
                    tags$div(
                        class = "visually-hidden",
                        textOutput("demo_data-data_ready"),
                        textOutput("demo_norm-norm_ready")
                    ),
                    card(
                        class = "mb-3 border-success",
                        card_header(
                            class = "bg-success text-white d-flex align-items-center gap-2",
                            icon("play-circle"),
                            tags$strong("Demo Mode — Example Dataset (2,000 cells)")
                        ),
                        card_body(
                            class = "py-2",
                            tags$p(
                                class = "mb-0",
                                "The demo dataset has been loaded automatically. ",
                                "Walk through each step below to explore a complete analysis. ",
                                "When you are ready to use your own data, switch to the ",
                                tags$strong("Single Sample"), " tab."
                            )
                        )
                    ),
                    mod_data_loading_ui("demo_data", demo_mode = TRUE),
                    conditionalPanel(
                        condition = "output['demo_data-data_ready'] == 'true'",
                        mod_qc_ui("demo_qc"),
                        mod_normalization_ui("demo_norm")
                    ),
                    conditionalPanel(
                        condition = "output['demo_norm-norm_ready'] == 'true'",
                        mod_clustering_ui("demo_clust", step_number = 4)
                    ),
                    conditionalPanel(
                        condition = "output['demo_clust-clustering_ready'] == 'true'",
                        mod_de_ui("demo_de", step_number = 5),
                        mod_visualization_ui("demo_viz", step_number = 6)
                    ),
                    uiOutput("demo_clust-rename_section")
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

    demo_trigger <- reactiveVal(0)

    observeEvent(input$go_to_demo, {
        bslib::nav_select("main_nav", "demo")
        demo_trigger(demo_trigger() + 1)
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
        list(id = "demo_data-load_demo", action = "Load demo data (demo tab)", prefix = "demo_data-"),
        list(id = "demo_qc-apply_filter", action = "QC filter (demo)", prefix = "demo_qc-"),
        list(id = "demo_norm-run_norm", action = "Normalize (demo)", prefix = "demo_norm-"),
        list(id = "demo_clust-run_clustering", action = "Cluster (demo)", prefix = "demo_clust-"),
        list(id = "demo_de-run_de", action = "Run DE (demo)", prefix = "demo_de-"),
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
        pkg_version <- function(pkg) {
            if (requireNamespace(pkg, quietly = TRUE)) {
                as.character(packageVersion(pkg))
            } else {
                "not installed"
            }
        }

        lines <- c(
            "# Asc-Seurat session report",
            "",
            paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
            paste("ascseurat version:", as.character(packageVersion("ascseurat"))),
            paste("R version:", R.version.string),
            paste("Seurat version:", pkg_version("Seurat")),
            paste("SeuratObject version:", pkg_version("SeuratObject")),
            paste("SingleR version:", pkg_version("SingleR")),
            paste("slingshot version:", pkg_version("slingshot")),
            "",
            "## Citation",
            "",
            "Pereira WJ, Almeida FM, Balmant KM, Rodriguez DC, Triozzi PM, Schmidt HW, Dervinis C, Pappas Jr. GJ, Kirst M. Asc-Seurat: analytical single-cell Seurat-based web application. BMC Bioinformatics 22, 556 (2021). https://doi.org/10.1186/s12859-021-04472-2",
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
                param_lines <- vapply(names(event$params), function(name) {
                    paste0("   - ", name, ": ", event$params[[name]])
                }, character(1))
                lines <- c(lines, param_lines)
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

    demo_data <- mod_data_loading_server("demo_data", trigger_demo = demo_trigger)
    demo_qc <- mod_qc_server("demo_qc", seurat_obj = demo_data)
    demo_norm <- mod_normalization_server("demo_norm", seurat_obj = demo_qc)
    demo_clust <- mod_clustering_server("demo_clust", norm_result = demo_norm)
    demo_de <- mod_de_server("demo_de", seurat_obj = demo_clust)
    mod_visualization_server("demo_viz", seurat_obj = demo_clust, markers_de = demo_de, de_step_number = 5)

    single_data <- mod_data_loading_server("single_data")
    single_qc <- mod_qc_server("single_qc", seurat_obj = single_data)
    single_norm <- mod_normalization_server("single_norm", seurat_obj = single_qc)
    single_clust <- mod_clustering_server("single_clust", norm_result = single_norm)
    single_de <- mod_de_server("single_de", seurat_obj = single_clust)
    mod_visualization_server("single_viz", seurat_obj = single_clust, markers_de = single_de, de_step_number = 5)

    multi_int <- mod_integration_server("multi_int")
    multi_clust <- mod_clustering_server("multi_clust", norm_result = multi_int)
    multi_de <- mod_de_server("multi_de", seurat_obj = multi_clust, allow_conserved = TRUE)
    mod_visualization_server("multi_viz", seurat_obj = multi_clust, markers_de = multi_de, de_step_number = 3)

    observeEvent(input$main_nav, {
        if (!identical(input$main_nav, "demo")) {
            return(invisible(NULL))
        }

        demo_missing <- tryCatch(
            is.null(demo_data()),
            error = function(e) TRUE
        )
        if (isTRUE(demo_missing)) {
            demo_trigger(demo_trigger() + 1)
        }
    }, ignoreInit = TRUE)

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

shinyApp(ui = ui, server = server)
