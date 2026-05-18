#' ascseurat: Analytical Single Cell Seurat-Based Web Application
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom shiny NS moduleServer reactive reactiveVal eventReactive req observe observeEvent isolate showNotification withProgress setProgress downloadHandler downloadButton fileInput textInput numericInput radioButtons selectInput selectizeInput checkboxGroupInput checkboxInput actionButton conditionalPanel uiOutput renderUI plotOutput renderPlot verbatimTextOutput renderPrint tagList tags icon HTML div addResourcePath updateNumericInput bookmarkButton enableBookmarking updateSelectizeInput updateRadioButtons setBookmarkExclude onBookmarked showModal modalDialog
#' @importFrom bslib page_navbar nav_panel nav_menu nav_spacer nav_item card card_header card_body card_footer layout_columns value_box bs_theme bs_add_rules input_dark_mode
#' @importFrom Seurat Read10X Read10X_h5 CreateSeuratObject NormalizeData FindVariableFeatures ScaleData RunPCA FindNeighbors FindClusters RunUMAP RunTSNE ElbowPlot PercentageFeatureSet FindAllMarkers FindMarkers FindConservedMarkers DimPlot FeaturePlot VlnPlot DotPlot DoHeatmap WhichCells SCTransform GetAssayData UpdateSeuratObject VariableFeatures IntegrateLayers RPCAIntegration as.SingleCellExperiment Embeddings CellCycleScoring
#' @importFrom SeuratObject Idents "Idents<-" DefaultAssay LayerData JoinLayers as.sparse
#' @importFrom dplyr filter mutate arrange group_by summarise "%>%" slice_max pull
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom reactable reactable reactableOutput renderReactable
#' @importFrom DT DTOutput renderDT datatable
#' @importFrom shinyWidgets pickerInput
#' @importFrom shinyFeedback useShinyFeedback showFeedbackDanger hideFeedback showToast
#' @importFrom shinycssloaders withSpinner
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom sass sass
#' @importFrom ggplot2 ggsave theme_void theme element_blank scale_fill_viridis_c ggtitle
#' @importFrom methods is slotNames
#' @importFrom stats p.adjust median
#' @importFrom utils read.csv write.csv packageVersion read.delim zip head tail
#' @importFrom grDevices png pdf svg tiff jpeg dev.off hcl.colors
NULL

utils::globalVariables(c(
    ".data",
    "Count",
    "avg_exp",
    "bin",
    "cluster",
    "gene",
    "nFeature_RNA",
    "pct_exp",
    "percent.mt",
    "scDblFinder_class",
    "value"
))

.onLoad <- function(libname, pkgname) {
    declare_paga_python_requirements()
}
