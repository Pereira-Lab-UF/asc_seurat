#' Normalization Module - UI
#'
#' @param id Module namespace ID.
#' @export
#' @keywords internal
mod_normalization_ui <- function(id) {
    ns <- NS(id)

    step_card(
        title = "Step 3: Normalization & PCA",
        next_step = "Review the elbow plot, choose the number of PCs, then proceed to Step 4: Clustering",

        layout_columns(
            col_widths = c(4, 4, 4),
            card(
                card_header("Normalization"),
                card_body(
                    input_norm_method(ns("norm_method")),
                    numericInput(
                        ns("n_variable_features"),
                        "Number of variable features",
                        value = 3000, min = 500, max = 10000, step = 500
                    ),
                    conditionalPanel(
                        condition = sprintf("input['%s'] == 'lognorm'", ns("norm_method")),
                        numericInput(
                            ns("scale_factor"),
                            "Scale factor",
                            value = 10000
                        )
                    ),
                    conditionalPanel(
                        condition = sprintf("input['%s'] == 'sctransform'", ns("norm_method")),
                        selectizeInput(
                            ns("vars_to_regress"),
                            "Variables to regress out (optional):",
                            choices = NULL,
                            multiple = TRUE,
                            options = list(
                                placeholder = "e.g. percent.mt"
                            )
                        )
                    )
                )
            ),
            card(
                card_header("Cell-cycle scoring (optional)"),
                card_body(
                    checkboxInput(ns("score_cc"),
                                  "Score cell cycle",
                                  value = FALSE),
                    conditionalPanel(
                        condition = sprintf("input['%s'] == true", ns("score_cc")),
                        selectInput(
                            ns("cc_organism"),
                            "Organism:",
                            choices = c("Human" = "human", "Mouse" = "mouse"),
                            selected = "human"
                        ),
                        checkboxInput(
                            ns("regress_cc"),
                            "Regress out cell-cycle scores during scaling",
                            value = FALSE
                        )
                    )
                )
            ),
            card(
                card_header("Run"),
                card_body(
                    tags$p("Select the normalization method and click Run."),
                    tags$p(
                        tags$strong("LogNormalize:"),
                        "Standard log-normalization."
                    ),
                    tags$p(
                        tags$strong("SCTransform v2:"),
                        "Variance-stabilizing normalization."
                    ),
                    tags$p(
                        class = "text-warning small fw-semibold mb-0",
                        "Apply filters in Step 2, choose No filtering, or apply doublet removal if that option is selected before running Normalization + PCA."
                    ),
                    tags$hr(),
                    action_btn(
                        ns("run_norm"),
                        "Run Normalization + PCA",
                        icon = icon("play")
                    )
                )
            )
        ),

        conditionalPanel(
            condition = sprintf("output['%s'] == 'true'", ns("norm_ready")),
            tags$hr(),
            layout_columns(
                col_widths = c(8, 4),
                card(
                    class = "normalization-elbow-card",
                    card_header("Elbow Plot - Select Number of PCs"),
                    card_body(
                        class = "plot-container plot-container-wide",
                        with_spinner(plotOutput(ns("elbow_plot"), height = "440px", width = "100%"))
                    )
                ),
                card(
                    class = "h-100",
                    card_header("PCA Dimensions"),
                    card_body(
                        numericInput(
                            ns("n_pcs"),
                            "Number of PCs to use for downstream analysis",
                            value = 30, min = 5, max = 100
                        ),
                        tags$p(
                            class = "text-muted",
                            "Choose the number of PCs that capture most ",
                            "of the variance (the 'elbow' in the plot)."
                        ),
                        card_plot_download("pca", ns("download_elbow"), ns = ns)
                    )
                )
            )
        )
    )
}


#' Normalization Module - Server
#'
#' @param id Module namespace ID.
#' @param seurat_obj Reactive Seurat object from QC module.
#' @return A list of reactives with the normalized Seurat object and chosen PCs.
#' @export
#' @keywords internal
mod_normalization_server <- function(id, seurat_obj) {
    moduleServer(id, function(input, output, session) {
        normalized_obj <- reactiveVal(NULL)
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

        observeEvent(seurat_obj(), {
            normalized_obj(NULL)
            metadata_choices <- grep(
                "^(orig\\.ident|seurat_clusters|RNA_snn)",
                colnames(seurat_obj()@meta.data),
                invert = TRUE,
                value = TRUE
            )
            updateSelectizeInput(
                session,
                "vars_to_regress",
                choices = metadata_choices,
                server = FALSE
            )
        }, ignoreNULL = TRUE)

        output$norm_ready <- shiny::renderText({
            if (is.null(normalized_obj())) "false" else "true"
        })
        shiny::outputOptions(output, "norm_ready", suspendWhenHidden = FALSE)

        observeEvent(input$run_norm, {
            obj_in <- tryCatch(
                isolate(seurat_obj()),
                error = function(e) NULL
            )
            if (is.null(obj_in)) {
                showNotification(
                    ui = tagList(
                        tags$strong("QC choice needed"),
                        tags$div("Apply filters in Step 2, choose No filtering, or apply doublet removal if that option is selected before running Normalization + PCA.")
                    ),
                    type = "warning",
                    duration = 6
                )
                return(invisible(NULL))
            }
            withProgress(message = "Normalizing data...", value = 0.2, {
                obj <- run_normalization(
                    obj = obj_in,
                    method = input$norm_method,
                    nfeatures = input$n_variable_features,
                    scale_factor = if (!is.null(input$scale_factor)) input$scale_factor else 10000,
                    vars_to_regress = input$vars_to_regress,
                    score_cell_cycle = isTRUE(input$score_cc),
                    cc_organism = if (!is.null(input$cc_organism)) input$cc_organism else "human",
                    regress_cc = isTRUE(input$regress_cc)
                )
                setProgress(1, message = "Normalization complete!")
                normalized_obj(obj)
            })
        })

        output$elbow_plot <- renderPlot({
            req(normalized_obj())
            old_par <- graphics::par(no.readonly = TRUE)
            on.exit(graphics::par(old_par), add = TRUE)
            graphics::par(mar = c(4.5, 4.5, 1.5, 1))
            ElbowPlot(normalized_obj(), ndims = 50)
        }, height = function() 440, width = function() plot_width("elbow_plot"))

        output$download_elbow <- plot_download_handler(
            plot_reactive = reactive({ ElbowPlot(normalized_obj(), ndims = 50) }),
            filename_base = "elbow_plot",
            input = input,
            prefix = "pca"
        )

        list(
            obj = reactive({
                req(normalized_obj())
                normalized_obj()
            }),
            n_pcs = reactive({
                req(normalized_obj(), input$n_pcs)
                input$n_pcs
            })
        )
    })
}
