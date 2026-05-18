#' Launch the Asc-Seurat Shiny Application
#'
#' Starts the Asc-Seurat interactive web application for single-cell
#' RNA-seq analysis. Creates necessary working directories if they
#' don't exist.
#'
#' @param host Character. The host address to serve the app on.
#'   Defaults to "127.0.0.1" (localhost). Use "0.0.0.0" for Docker.
#' @param port Integer. The port to serve the app on. Defaults to 3838.
#' @param launch.browser Logical. Whether to open the app in a browser.
#'   Defaults to TRUE.
#' @param workdir Character. Directory used for user data, generated files, and
#'   Shiny runtime state. Defaults to the current working directory.
#' @param demo Character. Bundled demo dataset identifier. Defaults to "pbmc".
#' @param max_upload_size Numeric. Maximum browser upload size in megabytes.
#'   Defaults to 5120 MB. Use file paths inside the app for very large RDS files.
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @export
#' @examples
#' \dontrun{
#' ascseurat::run_app()
#' }
run_app <- function(host = "127.0.0.1", port = 3838,
                    launch.browser = TRUE, workdir = getwd(),
                    demo = "pbmc", max_upload_size = 5120, ...) {

    demo <- match.arg(demo, c("pbmc"))
    if (!dir.exists(workdir)) {
        dir.create(workdir, recursive = TRUE)
    }
    workdir <- normalizePath(workdir, mustWork = TRUE)
    old_workdir <- setwd(workdir)
    on.exit(setwd(old_workdir), add = TRUE)
    options(
        ascseurat.demo = demo,
        ascseurat.workdir = workdir,
        shiny.maxRequestSize = max_upload_size * 1024^2
    )

    # Create working directories if they don't exist
    if (!dir.exists("data")) {
        dir.create("data", recursive = TRUE)
        message("Created 'data/' directory in: ", workdir)
    }
    if (!dir.exists("RDS_files")) {
        dir.create("RDS_files", recursive = TRUE)
        message("Created 'RDS_files/' directory in: ", workdir)
    }

    # Locate app directory
    app_dir <- system.file("app", package = "ascseurat")
    if (app_dir == "") {
        stop(
            "Could not find the app directory. ",
            "Try re-installing the `ascseurat` package.",
            call. = FALSE
        )
    }

    # Launch
    message("Starting Asc-Seurat v3...")
    message("App directory: ", app_dir)
    message("Working directory: ", workdir)
    shiny::runApp(
        appDir = app_dir,
        host = host,
        port = port,
        launch.browser = launch.browser,
        ...
    )
}


paga_python_packages <- function() {
    c("scanpy", "anndata", "numpy", "scipy", "pandas", "leidenalg", "igraph")
}


declare_paga_python_requirements <- function() {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        return(invisible(FALSE))
    }

    if (utils::packageVersion("reticulate") >= "1.41.0") {
        reticulate::py_require(paga_python_packages())
    }

    invisible(TRUE)
}


#' Set Up Python Environment for PAGA
#'
#' Installs and verifies the Python packages used by the supported PAGA
#' trajectory workflow. Docker images already include this stack.
#'
#' @param method Character. Installation method for Python.
#'   Defaults to "auto" which uses conda if available, otherwise virtualenv.
#' @param envname Character. Name of the Python environment.
#'   Defaults to "ascseurat-paga".
#'
#' @export
#' @examples
#' \dontrun{
#' ascseurat::setup_paga()
#' }
setup_paga <- function(method = "auto", envname = "ascseurat-paga") {

    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package 'reticulate' is required. Install with: install.packages('reticulate')",
             call. = FALSE)
    }

    message("Setting up Python environment for PAGA...")

    # Install miniconda if no Python available
    if (!reticulate::py_available(initialize = FALSE)) {
        message("No Python found. Installing Miniconda...")
        tryCatch(
            reticulate::install_miniconda(),
            error = function(e) {
                stop(
                    "Miniconda installation failed: ", conditionMessage(e),
                    "\n\nThis usually means no internet connection or insufficient disk space. ",
                    "You can install Python manually and re-run setup_paga().",
                    call. = FALSE
                )
            }
        )
    }

    # Create environment and install the Scanpy stack
    message("Installing the Scanpy/PAGA stack in environment '", envname, "'...")
    tryCatch(
        reticulate::py_install(
            packages = paga_python_packages(),
            envname = envname,
            method = method,
            pip = TRUE
        ),
        error = function(e) {
            stop(
                "PAGA Python dependency installation failed: ", conditionMessage(e),
                "\n\nCheck your internet connection and conda/pip configuration.",
                call. = FALSE
            )
        }
    )

    declare_paga_python_requirements()

    message("\nPAGA setup complete!")
    message("To use PAGA, the environment will be activated automatically when possible.")
    message("You can verify with: reticulate::py_module_available('scanpy')")

    invisible(TRUE)
}


#' Check if PAGA is Available
#'
#' @return Logical. TRUE if Python and scanpy are available.
#' @keywords internal
paga_available <- function() {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
    declare_paga_python_requirements()
    tryCatch({
        reticulate::py_available(initialize = TRUE) &&
            all(vapply(
                paga_python_packages(),
                reticulate::py_module_available,
                logical(1)
            ))
    }, error = function(e) FALSE)
}
