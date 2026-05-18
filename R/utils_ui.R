#' Numeric Input for Plot Height
#' @param id Input ID.
#' @param value Default value in cm.
#' @keywords internal
input_plot_height <- function(id, value = 15) {
    numericInput(id, "Height (cm)", step = 0.5, value = value)
}

#' Numeric Input for Plot Width
#' @param id Input ID.
#' @param value Default value in cm.
#' @keywords internal
input_plot_width <- function(id, value = 15) {
    numericInput(id, "Width (cm)", step = 0.5, value = value)
}

#' Select Input for Plot Resolution
#' @param id Input ID.
#' @keywords internal
input_plot_res <- function(id) {
    selectInput(id, "Resolution (DPI)",
                choices = c("100", "200", "300", "400", "500", "600"),
                selected = "300")
}

#' Select Input for Plot Format
#' @param id Input ID.
#' @keywords internal
input_plot_format <- function(id) {
    selectInput(id, "File type",
                choices = list("png" = "png", "tiff" = "tiff",
                               "jpeg" = "jpeg", "pdf" = "pdf", "svg" = "svg"),
                selected = "png")
}

#' Plot Download Options Card
#'
#' Creates a bslib card containing all download options for a plot.
#'
#' @param prefix Character. Prefix for all input IDs.
#' @param download_id Character. ID for the download button.
#' @param height_default Numeric. Default height in cm.
#' @param width_default Numeric. Default width in cm.
#' @param ns Namespace function for module-scoped inputs.
#' @param plot_choices Optional named character vector/list of plot choices.
#' @param plot_selected Optional default selected plot choice.
#' @param plot_label Character. Label for the optional plot-choice selector.
#' @return A bslib card UI element.
#' @keywords internal
card_plot_download <- function(prefix, download_id,
                               height_default = 10,
                               width_default = 15,
                               ns = identity,
                               plot_choices = NULL,
                               plot_selected = NULL,
                               plot_label = "Plot to download") {
    plot_selector <- NULL
    if (!is.null(plot_choices)) {
        if (is.null(plot_selected)) {
            plot_selected <- unname(unlist(plot_choices, use.names = FALSE))[[1]]
        }
        plot_selector <- selectInput(
            ns(paste0(prefix, "_plot_choice")),
            plot_label,
            choices = plot_choices,
            selected = plot_selected
        )
    }

    bslib::card(
        class = "asc-download-options-card",
        bslib::card_body(
            tags$details(
                class = "asc-collapse-panel asc-download-options",
                tags$summary(
                    class = "asc-collapse-summary",
                    icon("download", class = "me-2"),
                    "Download options"
                ),
                tags$div(
                    class = "asc-download-options-body",
                    plot_selector,
                    input_plot_height(ns(paste0(prefix, "_height")), value = height_default),
                    input_plot_width(ns(paste0(prefix, "_width")), value = width_default),
                    input_plot_res(ns(paste0(prefix, "_res"))),
                    input_plot_format(ns(paste0(prefix, "_format"))),
                    downloadButton(download_id, "Download Plot", class = "btn-download-soft w-100")
                )
            )
        )
    )
}

#' Download Button with Busy Feedback
#'
#' Prevents repeated clicks while a long-running download is being prepared.
#'
#' @param output_id Download output ID.
#' @param label Button label.
#' @param class CSS class.
#' @param busy_label Label shown while the file is being prepared.
#' @param timeout_ms Fallback time before the button is re-enabled.
#' @return A download button UI element.
#' @keywords internal
busy_download_button <- function(output_id, label,
                                 class = "btn-download-soft w-100",
                                 busy_label = "Preparing download...",
                                 timeout_ms = 120000) {
    tagList(
        downloadButton(
            output_id,
            label,
            class = paste(class, "asc-busy-download"),
            `data-busy-label` = busy_label,
            `data-busy-timeout` = timeout_ms
        ),
        tags$script(HTML("
(function() {
  if (window.ascDownloadBusyInstalled) return;
  window.ascDownloadBusyInstalled = true;

  function restoreDownloadButton(el) {
    if (!el) return;
    if (el.ascDownloadTimer) {
      clearTimeout(el.ascDownloadTimer);
      el.ascDownloadTimer = null;
    }
    el.classList.remove('disabled', 'asc-download-busy');
    el.removeAttribute('aria-disabled');
    el.removeAttribute('aria-busy');
    if (el.dataset.originalHtml) {
      el.innerHTML = el.dataset.originalHtml;
    }
  }

  function markDownloadBusy(el) {
    if (!el || el.classList.contains('asc-download-busy')) return;
    el.dataset.originalHtml = el.innerHTML;
    var busyLabel = el.dataset.busyLabel || 'Preparing download...';
    el.classList.add('disabled', 'asc-download-busy');
    el.setAttribute('aria-disabled', 'true');
    el.setAttribute('aria-busy', 'true');
    el.innerHTML = '<span class=\"spinner-border spinner-border-sm me-2\" role=\"status\" aria-hidden=\"true\"></span>' + busyLabel;
    var timeout = parseInt(el.dataset.busyTimeout || '120000', 10);
    el.ascDownloadTimer = setTimeout(function() {
      restoreDownloadButton(el);
    }, timeout);
  }

  document.addEventListener('click', function(event) {
    var el = event.target.closest('a.asc-busy-download');
    if (!el) return;
    if (el.classList.contains('asc-download-busy')) {
      event.preventDefault();
      event.stopImmediatePropagation();
      return false;
    }
    markDownloadBusy(el);
  }, true);

  if (window.Shiny && Shiny.addCustomMessageHandler) {
    Shiny.addCustomMessageHandler('ascDownloadReady', function(message) {
      restoreDownloadButton(document.getElementById(message.id));
    });
  }
})();
"))
    )
}

#' Step Card
#'
#' Creates a styled card for an analysis step with a header, body, and an
#' optional "next step" hint footer.
#'
#' @param title Character. Title for the step.
#' @param ... UI elements for the card body.
#' @param id Optional card ID.
#' @param full_screen Logical. Allow full-screen viewing.
#' @param next_step Optional character. If supplied, a muted footer is appended
#'   with a right-pointing arrow and this text to guide the user toward the
#'   next action (e.g. "Proceed to Step 3: Normalization").
#' @return A bslib card.
#' @keywords internal
step_card <- function(title, ..., id = NULL, full_screen = FALSE,
                      next_step = NULL) {
    footer_el <- if (!is.null(next_step)) {
        bslib::card_footer(
            class = "text-muted small text-end",
            icon("arrow-right", class = "me-1"),
            next_step
        )
    } else {
        NULL
    }

    bslib::card(
        id = id,
        full_screen = full_screen,
        bslib::card_header(tags$h4(title), class = "bg-primary text-white"),
        bslib::card_body(
            class = "p-3",
            ...
        ),
        footer_el
    )
}

#' Action Button (Styled)
#'
#' A styled action button consistent with the v3 theme.
#'
#' @param id Input ID.
#' @param label Button label.
#' @param icon An optional icon.
#' @param class CSS class. Defaults to "btn-primary".
#' @keywords internal
action_btn <- function(id, label, icon = NULL, class = "btn-primary") {
    actionButton(id, label = label, icon = icon,
                 class = paste("btn-block", class))
}

#' Select Normalization Method
#' @param id Input ID.
#' @keywords internal
input_norm_method <- function(id) {
    radioButtons(id,
                 "Normalization method",
                 choices = list("LogNormalize" = "lognorm",
                                "SCTransform v2" = "sctransform"),
                 selected = "lognorm")
}

#' Select Statistical Test for DE
#' @param id Input ID.
#' @keywords internal
input_de_test <- function(id) {
    shinyWidgets::pickerInput(
        inputId = id,
        label = "Statistical test",
        choices = c("wilcox", "bimod", "roc", "t", "MAST"),
        selected = "wilcox",
        multiple = FALSE
    )
}

#' File Input for Markers List
#' @param id Input ID.
#' @keywords internal
input_markers_file <- function(id) {
    tagList(
        fileInput(id,
                  label = "Upload marker gene list",
                  accept = c("text/csv", ".csv", "text/tsv", ".tsv",
                             "text/comma-separated-values",
                             "text/tab-separated-values")),
        preserve_file_input_scroll_js(id)
    )
}


#' Preserve Scroll Position Around File Input Dialogs
#'
#' @param id Input ID.
#' @keywords internal
preserve_file_input_scroll_js <- function(id) {
    tags$script(HTML(sprintf("
        (function() {
            var fileId = %s;
            var key = 'ascseurat-file-scroll-y-' + fileId;
            var activeKey = 'ascseurat-file-scroll-active';

            function fileContainer() {
                var input = document.getElementById(fileId);
                return input && input.closest ? input.closest('.shiny-input-container') : null;
            }

            function isFileInputClick(target) {
                var container = fileContainer();
                return !!(container && target && container.contains(target));
            }

            function rememberScroll() {
                var y = window.scrollY || document.documentElement.scrollTop || 0;
                window.sessionStorage.setItem(key, String(y));
                window.sessionStorage.setItem(activeKey, key);
            }

            function restoreScroll() {
                if (window.sessionStorage.getItem(activeKey) !== key) return;
                var y = window.sessionStorage.getItem(key);
                if (y === null) return;
                y = parseInt(y, 10);
                if (!Number.isFinite(y)) return;
                window.scrollTo({ top: y, left: 0, behavior: 'auto' });
            }

            function restoreSoon() {
                window.setTimeout(restoreScroll, 0);
                window.setTimeout(restoreScroll, 75);
                window.setTimeout(restoreScroll, 250);
            }

            document.addEventListener('mousedown', function(event) {
                if (isFileInputClick(event.target)) rememberScroll();
            }, true);

            document.addEventListener('click', function(event) {
                if (isFileInputClick(event.target)) rememberScroll();
            }, true);

            document.addEventListener('change', function(event) {
                if (event.target && event.target.id === fileId) {
                    rememberScroll();
                    restoreSoon();
                }
            }, true);
        })();
    ", shQuote(id))))
}

#' Mitochondrial Gene Regex Input
#' @param id Input ID.
#' @param value Default regex pattern.
#' @keywords internal
input_mito_regex <- function(id, value = "") {
    textInput(
        id,
        label = "Mitochondrial gene pattern (regex)",
        value = value,
        placeholder = "Optional, e.g. ^MT-"
    )
}
