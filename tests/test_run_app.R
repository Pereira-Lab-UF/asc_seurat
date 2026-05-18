args <- commandArgs(trailingOnly = TRUE)
port <- if (length(args)) as.integer(args[[1]]) else 17775L
host <- "127.0.0.1"

if (!requireNamespace("callr", quietly = TRUE)) {
    stop("Package 'callr' is required to run tests/test_run_app.R", call. = FALSE)
}

app_proc <- callr::r_bg(
    func = function(port) {
        ascseurat::run_app(
            host = "127.0.0.1",
            port = port,
            launch.browser = FALSE
        )
    },
    args = list(port = port),
    supervise = TRUE,
    stdout = "|",
    stderr = "|"
)

on.exit({
    if (app_proc$is_alive()) {
        app_proc$kill()
    }
}, add = TRUE)

deadline <- Sys.time() + 30
html <- NULL

while (Sys.time() < deadline) {
    if (!app_proc$is_alive()) {
        stop(
            "ascseurat::run_app() exited before the app became reachable.\n",
            paste(app_proc$read_all_error(), collapse = "\n"),
            call. = FALSE
        )
    }

    html <- tryCatch(
        paste(readLines(sprintf("http://%s:%s", host, port), warn = FALSE), collapse = "\n"),
        error = function(e) NULL
    )

    if (!is.null(html) && nzchar(html)) {
        break
    }

    Sys.sleep(1)
}

if (is.null(html) || !nzchar(html)) {
    stop("Timed out waiting for Asc-Seurat to respond on the configured port.", call. = FALSE)
}

required_strings <- c(
    "Asc-Seurat v3",
    "Home",
    "Single Sample",
    "Integration",
    "Trajectory Inference",
    "Cell-type Annotation",
    "Advanced Plots",
    "Annotation"
)

missing_strings <- required_strings[!vapply(required_strings, grepl, logical(1), x = html, fixed = TRUE)]
if (length(missing_strings)) {
    stop(
        "Root page is missing expected content: ",
        paste(missing_strings, collapse = ", "),
        call. = FALSE
    )
}

message("Asc-Seurat runtime smoke test passed on port ", port, ".")
