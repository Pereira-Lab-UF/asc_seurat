#' Open a Graphics Device
#'
#' Opens a graphics device with cm-based sizing.
#'
#' @param format Character. One of "png", "tiff", "jpeg", "pdf", "svg".
#' @param file Character. Output file path.
#' @param width Numeric. Width in cm.
#' @param height Numeric. Height in cm.
#' @param dpi Integer. Resolution in DPI.
#' @keywords internal
open_plot_device <- function(file, format = "png",
                             width = 15, height = 10, dpi = 300) {
    switch(
        format,
        png = grDevices::png(
            filename = file,
            width = width,
            height = height,
            units = "cm",
            res = dpi
        ),
        tiff = grDevices::tiff(
            filename = file,
            width = width,
            height = height,
            units = "cm",
            res = dpi,
            compression = "lzw"
        ),
        jpeg = grDevices::jpeg(
            filename = file,
            width = width,
            height = height,
            units = "cm",
            res = dpi,
            quality = 100
        ),
        pdf = grDevices::pdf(
            file = file,
            width = width / 2.54,
            height = height / 2.54
        ),
        svg = grDevices::svg(
            filename = file,
            width = width / 2.54,
            height = height / 2.54
        ),
        stop("Unsupported plot format: ", format, call. = FALSE)
    )
}


#' Write a Plot to Disk
#'
#' Saves either a ggplot/patchwork object or a drawing function.
#'
#' @param plot The plot object to save, or a zero-argument drawing function.
#' @param file Character. Output file path.
#' @param format Character. One of "png", "tiff", "jpeg", "pdf", "svg".
#' @param width Numeric. Width in cm.
#' @param height Numeric. Height in cm.
#' @param dpi Integer. Resolution in DPI.
#' @keywords internal
write_plot_file <- function(plot, file, format = "png",
                            width = 15, height = 10, dpi = 300) {
    if (is.function(plot)) {
        open_plot_device(
            file = file,
            format = format,
            width = width,
            height = height,
            dpi = dpi
        )
        on.exit(grDevices::dev.off(), add = TRUE)
        plot()
    } else {
        ggplot2::ggsave(
            filename = file,
            plot = plot,
            width = width,
            height = height,
            units = "cm",
            dpi = dpi,
            device = format
        )
    }

    invisible(file)
}


#' Save Plot with Configurable Options
#'
#' Wrapper around ggsave/device-based saving with support for
#' multiple formats and configurable resolution.
#'
#' @param plot The plot object to save, or a zero-argument drawing function.
#' @param filename Character. Output filename (without extension).
#' @param format Character. One of "png", "tiff", "jpeg", "pdf", "svg".
#' @param width Numeric. Width in cm.
#' @param height Numeric. Height in cm.
#' @param dpi Integer. Resolution in DPI.
#' @param path Character. Directory to save in. Defaults to current directory.
#' @return The file path of the saved plot (invisibly).
#' @keywords internal
save_plot <- function(plot, filename, format = "png",
                      width = 15, height = 10, dpi = 300,
                      path = ".") {

    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }

    full_path <- file.path(path, paste0(filename, ".", format))
    write_plot_file(
        plot = plot,
        file = full_path,
        format = format,
        width = width,
        height = height,
        dpi = dpi
    )
}


#' Create a Download Handler for Plots
#'
#' Factory function that creates a Shiny downloadHandler for plots
#' with user-configurable format, size, and resolution.
#'
#' @param plot_reactive A reactive expression that returns the plot.
#' @param filename_base Character. Base filename (no extension).
#' @param input The Shiny input object.
#' @param prefix Character. Prefix for input IDs (e.g., "p1" matches
#'   "p1_height", "p1_width", "p1_res", "p1_format").
#' @return A \code{\link[shiny]{downloadHandler}}.
#' @keywords internal
plot_download_handler <- function(plot_reactive, filename_base, input, prefix) {

    shiny::downloadHandler(
        filename = function() {
            fmt <- input[[paste0(prefix, "_format")]]
            paste0(filename_base, "_", Sys.Date(), ".", fmt)
        },
        content = function(file) {
            plt <- plot_reactive()
            fmt <- input[[paste0(prefix, "_format")]]
            w   <- input[[paste0(prefix, "_width")]]
            h   <- input[[paste0(prefix, "_height")]]
            res <- as.numeric(input[[paste0(prefix, "_res")]])

            write_plot_file(
                plot = plt,
                file = file,
                format = fmt,
                width = w,
                height = h,
                dpi = res
            )
        }
    )
}


#' Custom Spinner Wrapper
#'
#' Wraps a UI element with a loading spinner.
#'
#' @param ui_element A Shiny UI element.
#' @param hide_ui Logical. Whether to hide the UI while loading.
#' @return The UI element wrapped with a spinner.
#' @keywords internal
with_spinner <- function(ui_element, hide_ui = FALSE) {
    shinycssloaders::withSpinner(
        ui_element,
        type = 4,
        color = "#6366f1",
        size = 0.75,
        hide.ui = hide_ui
    )
}
