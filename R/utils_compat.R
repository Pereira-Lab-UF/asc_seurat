#' Load a Seurat RDS File with Backward Compatibility
#'
#' Loads an RDS file and automatically updates old Seurat objects
#' (v2/v4) to the current Seurat v5 format.
#'
#' @param path Character. Path to the RDS file.
#' @return A Seurat object (v5 format).
#' @keywords internal
load_seurat_rds <- function(path) {
    if (!file.exists(path)) {
        stop("File not found: ", path, call. = FALSE)
    }

    message("Loading RDS file: ", basename(path))
    obj <- readRDS(path)

    # Check if it's a Seurat object
    if (!inherits(obj, "Seurat")) {
        stop("The loaded RDS file does not contain a Seurat object.", call. = FALSE)
    }

    # Update old Seurat objects to v5
    tryCatch({
        obj_version <- obj@version
        current_version <- packageVersion("SeuratObject")
        if (obj_version < current_version) {
            message("Updating Seurat object from v", obj_version, " to v", current_version, "...")
            obj <- UpdateSeuratObject(obj)
            message("Update complete.")
        }
    }, error = function(e) {
        message("Attempting object update due to version mismatch...")
        obj <<- UpdateSeuratObject(obj)
    })

    return(obj)
}


# Null-coalescing helper.
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}


#' Calculate Robust QC Limits
#'
#' Uses median +/- n MADs for numeric QC metrics.
#'
#' @param values Numeric vector.
#' @param nmads Number of MADs from the median.
#' @param lower_floor Minimum allowed lower bound.
#' @param upper_ceiling Maximum allowed upper bound.
#' @return Named numeric vector with lower and upper values.
#' @keywords internal
qc_mad_limits <- function(values,
                          nmads = 3,
                          lower_floor = 0,
                          upper_ceiling = Inf) {
    values <- suppressWarnings(as.numeric(values))
    values <- values[is.finite(values)]
    if (!length(values)) {
        return(c(lower = NA_real_, upper = NA_real_))
    }

    center <- stats::median(values, na.rm = TRUE)
    spread <- stats::mad(values, center = center, na.rm = TRUE)
    if (!is.finite(spread) || spread <= 0) {
        spread <- stats::IQR(values, na.rm = TRUE) / 1.349
    }
    if (!is.finite(spread) || spread <= 0) {
        spread <- stats::sd(values, na.rm = TRUE)
    }
    if (!is.finite(spread) || spread <= 0) {
        spread <- 0
    }

    lower <- center - (nmads * spread)
    upper <- center + (nmads * spread)
    lower <- max(lower_floor, lower, na.rm = TRUE)
    upper <- min(upper_ceiling, upper, na.rm = TRUE)

    if (!is.finite(lower)) lower <- NA_real_
    if (!is.finite(upper)) upper <- NA_real_
    c(lower = lower, upper = upper)
}


#' Build QC Filter Defaults from Metadata
#'
#' @param meta Seurat metadata data frame.
#' @param nmads Number of MADs from the median.
#' @param defaults Named fallback values.
#' @return Named list of QC default values.
#' @keywords internal
qc_filter_defaults <- function(meta,
                               nmads = 3,
                               defaults = list(
                                   min_genes = 200,
                                   max_genes = 5000,
                                   max_mito = 5,
                                   min_counts = 0,
                                   max_counts = NA_real_
                               )) {
    result <- defaults

    if ("nFeature_RNA" %in% colnames(meta)) {
        limits <- qc_mad_limits(meta$nFeature_RNA, nmads = nmads, lower_floor = 0)
        if (is.finite(limits[["lower"]])) {
            result$min_genes <- max(0, floor(limits[["lower"]]))
        }
        if (is.finite(limits[["upper"]])) {
            result$max_genes <- max(result$min_genes + 1, ceiling(limits[["upper"]]))
        }
    }

    if ("nCount_RNA" %in% colnames(meta)) {
        limits <- qc_mad_limits(meta$nCount_RNA, nmads = nmads, lower_floor = 0)
        if (is.finite(limits[["lower"]])) {
            result$min_counts <- max(0, floor(limits[["lower"]]))
        }
        if (is.finite(limits[["upper"]])) {
            result$max_counts <- max(result$min_counts + 1, ceiling(limits[["upper"]]))
        }
    }

    if ("percent.mt" %in% colnames(meta)) {
        mito_values <- suppressWarnings(as.numeric(meta$percent.mt))
        limits <- qc_mad_limits(
            mito_values,
            nmads = nmads,
            lower_floor = 0,
            upper_ceiling = 100
        )
        if (isTRUE(any(mito_values > 0, na.rm = TRUE)) &&
            is.finite(limits[["upper"]]) &&
            limits[["upper"]] > 0) {
            result$max_mito <- min(100, max(5, round(limits[["upper"]], 1)))
        }
    }

    result
}


#' Read a Delimited Text Table
#'
#' Reads CSV, TSV, or TXT input based on file extension.
#'
#' @param path Character. Path to the table.
#' @param header Logical. Whether the file has a header row.
#' @param row.names Optional row names column.
#' @param check.names Logical. Passed through to the reader.
#' @param stringsAsFactors Logical. Passed through to the reader.
#' @return A data frame.
#' @keywords internal
read_uploaded_table <- function(path,
                                header = TRUE,
                                row.names = NULL,
                                check.names = FALSE,
                                stringsAsFactors = FALSE) {
    ext <- tolower(tools::file_ext(path))

    if (ext %in% c("tsv", "txt")) {
        utils::read.delim(
            path,
            header = header,
            row.names = row.names,
            check.names = check.names,
            stringsAsFactors = stringsAsFactors
        )
    } else {
        utils::read.csv(
            path,
            header = header,
            row.names = row.names,
            check.names = check.names,
            stringsAsFactors = stringsAsFactors
        )
    }
}


#' Detect Available Memory in Bytes
#'
#' Attempts to detect currently available memory, respecting container limits
#' when possible.
#'
#' @return Numeric scalar with available bytes, or NA_real_ if unavailable.
#' @keywords internal
detect_available_memory_bytes <- function() {
    read_numeric_file <- function(path) {
        if (!file.exists(path)) {
            return(NA_real_)
        }

        value <- trimws(readLines(path, warn = FALSE, n = 1))
        if (identical(value, "max")) {
            return(NA_real_)
        }

        suppressWarnings(as.numeric(value))
    }

    parse_vm_stat_value <- function(lines, label) {
        line <- grep(paste0("^", label), lines, value = TRUE)
        if (!length(line)) {
            return(0)
        }
        suppressWarnings(as.numeric(gsub("[^0-9]", "", line[[1]])))
    }

    cgroup_limit <- NA_real_
    cgroup_usage <- NA_real_

    # cgroup v2
    if (file.exists("/sys/fs/cgroup/memory.max")) {
        cgroup_limit <- read_numeric_file("/sys/fs/cgroup/memory.max")
        cgroup_usage <- read_numeric_file("/sys/fs/cgroup/memory.current")
    }

    # cgroup v1
    if (!is.finite(cgroup_limit) &&
        file.exists("/sys/fs/cgroup/memory/memory.limit_in_bytes")) {
        cgroup_limit <- read_numeric_file("/sys/fs/cgroup/memory/memory.limit_in_bytes")
        cgroup_usage <- read_numeric_file("/sys/fs/cgroup/memory/memory.usage_in_bytes")
    }

    if (is.finite(cgroup_limit) && is.finite(cgroup_usage) &&
        cgroup_limit > cgroup_usage && cgroup_limit < 1e18) {
        return(cgroup_limit - cgroup_usage)
    }

    if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] == "Linux" &&
        file.exists("/proc/meminfo")) {
        meminfo <- readLines("/proc/meminfo", warn = FALSE)
        line <- grep("^MemAvailable:", meminfo, value = TRUE)
        if (!length(line)) {
            line <- grep("^MemTotal:", meminfo, value = TRUE)
        }
        if (length(line)) {
            kb <- suppressWarnings(as.numeric(gsub("[^0-9]", "", line[[1]])))
            if (is.finite(kb)) {
                return(kb * 1024)
            }
        }
    }

    if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] == "Darwin") {
        vm_stat <- tryCatch(
            system2("vm_stat", stdout = TRUE, stderr = FALSE),
            error = function(e) character()
        )

        if (length(vm_stat)) {
            page_size_match <- regexec("page size of ([0-9]+) bytes", vm_stat[[1]])
            page_size_tokens <- regmatches(vm_stat[[1]], page_size_match)[[1]]
            page_size <- if (length(page_size_tokens) >= 2) {
                suppressWarnings(as.numeric(page_size_tokens[[2]]))
            } else {
                4096
            }

            free_pages <- parse_vm_stat_value(vm_stat, "Pages free")
            speculative_pages <- parse_vm_stat_value(vm_stat, "Pages speculative")
            inactive_pages <- parse_vm_stat_value(vm_stat, "Pages inactive")
            purgeable_pages <- parse_vm_stat_value(vm_stat, "Pages purgeable")

            bytes <- (free_pages + speculative_pages + inactive_pages + purgeable_pages) * page_size
            if (is.finite(bytes) && bytes > 0) {
                return(bytes)
            }
        }

        memsize <- tryCatch(
            system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE, stderr = FALSE),
            error = function(e) character()
        )
        if (length(memsize)) {
            bytes <- suppressWarnings(as.numeric(memsize[[1]]))
            if (is.finite(bytes)) {
                return(bytes)
            }
        }
    }

    NA_real_
}


#' Format Bytes for Human Readability
#'
#' @param bytes Numeric scalar in bytes.
#' @return Character string.
#' @keywords internal
format_bytes <- function(bytes) {
    if (!is.finite(bytes) || is.na(bytes)) {
        return("unknown")
    }

    units <- c("B", "KiB", "MiB", "GiB", "TiB")
    idx <- 1L
    value <- as.numeric(bytes)

    while (value >= 1024 && idx < length(units)) {
        value <- value / 1024
        idx <- idx + 1L
    }

    sprintf("%.2f %s", value, units[[idx]])
}


#' Get the Asc-Seurat Working Directory
#'
#' Shiny may evaluate the app from \code{inst/app}; this helper keeps user
#' paths anchored to the directory supplied to \code{run_app(workdir = ...)}.
#'
#' @return Normalized working directory path.
#' @keywords internal
ascseurat_workdir <- function() {
    normalizePath(
        getOption("ascseurat.workdir", getwd()),
        winslash = "/",
        mustWork = FALSE
    )
}


#' Temporarily Increase future.globals.maxSize
#'
#' Computes a memory-aware limit and applies it for the current scope.
#'
#' @param fraction Numeric fraction of detected available memory to permit.
#' @param reserve_bytes Numeric bytes to leave unallocated.
#' @param minimum_bytes Numeric minimum limit.
#' @param fallback_bytes Numeric fallback when memory cannot be detected.
#' @return A list containing the previous option value and the applied limit.
#' @keywords internal
local_future_globals_max_size <- function(fraction = 0.95,
                                          reserve_bytes = 512 * 1024^2,
                                          minimum_bytes = 2 * 1024^3,
                                          fallback_bytes = 8 * 1024^3) {
    env_override <- suppressWarnings(as.numeric(Sys.getenv(
        "ASCSEURAT_FUTURE_GLOBALS_MAXSIZE",
        unset = NA_character_
    )))

    detected_available <- detect_available_memory_bytes()

    target_bytes <- if (is.finite(env_override) && env_override > 0) {
        env_override
    } else if (is.finite(detected_available) && detected_available > 0) {
        usable_bytes <- max(detected_available - reserve_bytes, minimum_bytes)
        max(minimum_bytes, floor(min(usable_bytes, detected_available * fraction)))
    } else {
        fallback_bytes
    }

    current_value <- getOption("future.globals.maxSize", NA_real_)
    if (is.finite(current_value) && current_value > 0) {
        target_bytes <- max(target_bytes, current_value)
    }

    old_options <- options(future.globals.maxSize = target_bytes)

    list(
        old = old_options,
        target = target_bytes,
        available = detected_available
    )
}


#' Detect Input Data Format
#'
#' Determines the format of input data (10X directory, h5, h5ad, csv, rds)
#' and returns a standardized format identifier.
#'
#' @param path Character. Path to the input file or directory.
#' @return Character. One of: "10x_dir", "10x_h5", "h5ad", "csv", "tsv", "rds", "unknown".
#' @keywords internal
detect_input_format <- function(path) {
    if (dir.exists(path)) {
        # Check if it's a 10X directory (contains matrix.mtx, barcodes.tsv, genes/features.tsv)
        files <- list.files(path)
        has_mtx <- any(grepl("matrix\\.mtx", files))
        has_barcodes <- any(grepl("barcodes\\.tsv", files))
        has_features <- any(grepl("(genes|features)\\.tsv", files))
        if (has_mtx && has_barcodes && has_features) {
            return("10x_dir")
        }
        return("unknown")
    }

    if (!file.exists(path)) {
        return("unknown")
    }

    ext <- tolower(tools::file_ext(path))

    switch(ext,
        "h5"   = "10x_h5",
        "h5ad" = "h5ad",
        "csv"  = "csv",
        "tsv"  = "tsv",
        "rds"  = "rds",
        "unknown"
    )
}


#' Resolve a User-Supplied Input Path
#'
#' Relative paths are resolved against the app working directory.
#'
#' @param path Character. User-supplied file or directory path.
#' @return Normalized path when possible.
#' @keywords internal
resolve_input_path <- function(path) {
    path <- trimws(path.expand(path))
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
        path <- file.path(ascseurat_workdir(), path)
    }

    normalizePath(path, winslash = "/", mustWork = FALSE)
}


#' Load Single-Cell Data from Multiple Formats
#'
#' Loads scRNA-seq data from various file formats and returns a Seurat object.
#'
#' @param path Character. Path to the input file or directory.
#' @param project Character. Project name for the Seurat object.
#' @param min.cells Integer. Minimum number of cells expressing a gene (for filtering).
#' @param min.features Integer. Minimum number of features (genes) per cell.
#' @return A Seurat object.
#' @keywords internal
load_input_data <- function(path, project = "AscSeurat",
                            min.cells = 3, min.features = 200) {

    format <- detect_input_format(path)

    obj <- switch(format,
        "10x_dir" = {
            message("Loading 10X Genomics data from directory...")
            read_dir <- path
            if (file.exists(file.path(path, "features.tsv")) &&
                !file.exists(file.path(path, "features.tsv.gz"))) {
                read_dir <- tempfile("ascseurat_10x_")
                dir.create(read_dir, recursive = TRUE, showWarnings = FALSE)
                on.exit(unlink(read_dir, recursive = TRUE), add = TRUE)

                gzip_copy <- function(from, to) {
                    input <- file(from, open = "rb")
                    output <- gzfile(to, open = "wb")
                    on.exit(close(input), add = TRUE)
                    on.exit(close(output), add = TRUE)

                    repeat {
                        chunk <- readBin(input, what = "raw", n = 1024 * 1024)
                        if (!length(chunk)) {
                            break
                        }
                        writeBin(chunk, output)
                    }
                }

                for (filename in c("matrix.mtx", "barcodes.tsv", "features.tsv")) {
                    source_file <- file.path(path, filename)
                    if (!file.exists(source_file)) {
                        stop("Incomplete 10X directory; missing ", filename, call. = FALSE)
                    }
                    gzip_copy(source_file, file.path(read_dir, paste0(filename, ".gz")))
                }
            }
            data <- Seurat::Read10X(data.dir = read_dir)
            Seurat::CreateSeuratObject(
                counts = data,
                project = project,
                min.cells = min.cells,
                min.features = min.features
            )
        },

        "10x_h5" = {
            message("Loading 10X Genomics h5 file...")
            data <- Seurat::Read10X_h5(filename = path)
            Seurat::CreateSeuratObject(
                counts = data,
                project = project,
                min.cells = min.cells,
                min.features = min.features
            )
        },

        "h5ad" = {
            message("Loading h5ad (AnnData) file...")
            if (!requireNamespace("anndataR", quietly = TRUE)) {
                stop("Reading .h5ad files requires the 'anndataR' package. ",
                     "Install it with: install.packages('anndataR')",
                     call. = FALSE)
            }
            coerce_h5ad_seurat <- function(seurat_obj) {
                layers <- SeuratObject::Layers(seurat_obj)
                source_layer <- intersect(c("counts", "data", "X"), layers)
                if (!length(source_layer)) {
                    source_layer <- layers[[1]]
                } else {
                    source_layer <- source_layer[[1]]
                }

                counts <- Seurat::GetAssayData(seurat_obj, layer = source_layer)
                if (!inherits(counts, "Matrix")) {
                    warning(
                        "The AnnData X matrix is not sparse; converting it may use substantial memory.",
                        call. = FALSE
                    )
                    counts <- SeuratObject::as.sparse(as.matrix(counts))
                }

                metadata <- seurat_obj@meta.data
                metadata <- metadata[colnames(counts), , drop = FALSE]

                Seurat::CreateSeuratObject(
                    counts = counts,
                    project = project,
                    min.cells = min.cells,
                    min.features = min.features,
                    meta.data = metadata
                )
            }

            adata <- anndataR::read_h5ad(path, as = "Seurat")
            coerce_h5ad_seurat(adata)
        },

        "csv" = {
            message("Loading CSV count matrix...")
            data <- read_uploaded_table(path, row.names = 1, check.names = FALSE)
            Seurat::CreateSeuratObject(
                counts = SeuratObject::as.sparse(as.matrix(data)),
                project = project,
                min.cells = min.cells,
                min.features = min.features
            )
        },

        "tsv" = {
            message("Loading TSV count matrix...")
            data <- read_uploaded_table(path, row.names = 1, check.names = FALSE)
            Seurat::CreateSeuratObject(
                counts = SeuratObject::as.sparse(as.matrix(data)),
                project = project,
                min.cells = min.cells,
                min.features = min.features
            )
        },

        "rds" = {
            message("Loading pre-existing Seurat RDS file...")
            load_seurat_rds(path)
        },

        stop("Unsupported input format. Supported: 10X directory, .h5, .h5ad, .csv, .tsv, .rds",
             call. = FALSE)
    )

    message("Loaded ", ncol(obj), " cells and ", nrow(obj), " features.")
    return(obj)
}
