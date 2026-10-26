test_that("detect_input_format recognizes common formats", {
    csv_path <- tempfile(fileext = ".csv")
    tsv_path <- tempfile(fileext = ".tsv")
    rds_path <- tempfile(fileext = ".rds")

    write.csv(data.frame(gene = "A", cell = 1), csv_path, row.names = FALSE)
    write.table(data.frame(gene = "A", cell = 1), tsv_path,
                sep = "\t", row.names = FALSE)
    saveRDS(demo_fixture(), rds_path)

    expect_equal(ascseurat:::detect_input_format(csv_path), "csv")
    expect_equal(ascseurat:::detect_input_format(tsv_path), "tsv")
    expect_equal(ascseurat:::detect_input_format(rds_path), "rds")
    expect_equal(ascseurat:::detect_input_format(tempfile(fileext = ".foo")), "unknown")
})

test_that("resolve_input_path tolerates pasted shell quotes", {
    old_workdir <- getOption("ascseurat.workdir", NULL)
    on.exit({
        if (is.null(old_workdir)) {
            options(ascseurat.workdir = NULL)
        } else {
            options(ascseurat.workdir = old_workdir)
        }
    }, add = TRUE)

    workdir <- tempdir()
    options(ascseurat.workdir = workdir)

    absolute_path <- file.path(tempdir(), "sample folder")
    expect_equal(
        ascseurat:::resolve_input_path(paste0('"', absolute_path, '"')),
        normalizePath(absolute_path, winslash = "/", mustWork = FALSE)
    )
    expect_equal(
        ascseurat:::resolve_input_path(paste0("'", absolute_path)),
        normalizePath(absolute_path, winslash = "/", mustWork = FALSE)
    )
    expect_equal(
        ascseurat:::resolve_input_path("'relative folder'"),
        normalizePath(file.path(ascseurat:::ascseurat_workdir(), "relative folder"),
                      winslash = "/", mustWork = FALSE)
    )
})

test_that("format_bytes and local_future_globals_max_size return usable values", {
    expect_match(ascseurat:::format_bytes(1024), "KiB")
    info <- ascseurat:::local_future_globals_max_size()
    expect_type(info$target, "double")
    expect_true(info$target > 0)
})

test_that("QC mitochondrial default never starts at zero", {
    low_mito <- data.frame(percent.mt = c(0, 0.001, 0.002, 0.003, 0.004))
    zero_mito <- data.frame(percent.mt = rep(0, 5))

    expect_equal(ascseurat:::qc_filter_defaults(low_mito)$max_mito, 5)
    expect_equal(ascseurat:::qc_filter_defaults(zero_mito)$max_mito, 5)
})

h5ad_test_counts <- function() {
    set.seed(42)
    counts <- matrix(
        as.numeric(stats::rpois(30 * 25, 2)), nrow = 30, ncol = 25,
        dimnames = list(paste0("cell", 1:30), paste0("gene", 1:25))
    )
    counts[, 1] <- counts[, 1] + 1
    counts
}

# Scanpy-style normalize_total(target_sum = 1e4) followed by log1p.
h5ad_test_lognorm <- function(counts) {
    log1p(counts / rowSums(counts) * 1e4)
}

# Writes a cells x genes AnnData file; raw_X is added as .raw following the
# AnnData on-disk spec, since anndataR cannot write .raw itself.
write_test_h5ad <- function(X, layers = NULL, raw_X = NULL) {
    path <- tempfile(fileext = ".h5ad")
    obs <- data.frame(row.names = rownames(X), sample = rep("s1", nrow(X)))
    anndataR::AnnData(
        X = X,
        obs = obs,
        var = data.frame(row.names = colnames(X)),
        layers = layers
    )$write_h5ad(path)

    if (!is.null(raw_X)) {
        raw_path <- tempfile(fileext = ".h5ad")
        anndataR::AnnData(
            X = raw_X,
            obs = data.frame(row.names = rownames(raw_X)),
            var = data.frame(row.names = colnames(raw_X))
        )$write_h5ad(raw_path)
        h5 <- hdf5r::H5File$new(path, mode = "r+")
        raw_h5 <- hdf5r::H5File$new(raw_path, mode = "r")
        raw <- h5$create_group("raw")
        raw$obj_copy_from(raw_h5, "X", "X")
        raw$obj_copy_from(raw_h5, "var", "var")
        hdf5r::h5attr(raw, "encoding-type") <- "raw"
        hdf5r::h5attr(raw, "encoding-version") <- "0.1.0"
        raw_h5$close_all()
        h5$close_all()
    }
    path
}

load_test_h5ad <- function(path) {
    suppressWarnings(
        ascseurat:::load_input_data(path, project = "h5ad", min.cells = 0, min.features = 0)
    )
}

skip_without_h5ad_writer <- function() {
    skip_if_not_installed("anndataR")
    skip_if_not_installed("rhdf5")
    skip_if_not_installed("hdf5r")
}

test_that("h5ad loading reads raw counts from .raw when X is normalized", {
    skip_without_h5ad_writer()
    counts <- h5ad_test_counts()
    # X keeps only a gene subset, as after highly variable gene selection.
    normalized_x <- methods::as(h5ad_test_lognorm(counts)[, 1:15], "CsparseMatrix")
    raw_inputs <- list(
        csr = methods::as(counts, "RsparseMatrix"),
        csc = methods::as(counts, "CsparseMatrix"),
        dense = counts
    )

    for (encoding in names(raw_inputs)) {
        obj <- load_test_h5ad(write_test_h5ad(normalized_x, raw_X = raw_inputs[[encoding]]))
        expect_equal(
            as.matrix(SeuratObject::LayerData(obj, layer = "counts")),
            t(counts),
            info = encoding
        )
        expect_equal(obj@meta.data$sample, rep("s1", nrow(counts)), info = encoding)
        expect_equal(unname(obj$nCount_RNA), unname(rowSums(counts)), info = encoding)
        expect_equal(unname(obj$nFeature_RNA), unname(rowSums(counts > 0)), info = encoding)
    }
})

test_that("h5ad loading prefers a raw-counts layer over a normalized X", {
    skip_without_h5ad_writer()
    counts <- h5ad_test_counts()
    normalized_x <- methods::as(h5ad_test_lognorm(counts), "CsparseMatrix")

    for (layer in c("counts", "raw_counts")) {
        layers <- stats::setNames(list(methods::as(counts, "CsparseMatrix")), layer)
        obj <- load_test_h5ad(write_test_h5ad(normalized_x, layers = layers))
        expect_equal(
            as.matrix(SeuratObject::LayerData(obj, layer = "counts")),
            t(counts),
            info = layer
        )
    }
})

test_that("h5ad loading keeps integer X and rejects normalized X without counts", {
    skip_without_h5ad_writer()
    counts <- h5ad_test_counts()
    normalized_x <- methods::as(h5ad_test_lognorm(counts), "CsparseMatrix")

    obj <- load_test_h5ad(write_test_h5ad(methods::as(counts, "CsparseMatrix")))
    expect_equal(as.matrix(SeuratObject::LayerData(obj, layer = "counts")), t(counts))

    expect_error(
        load_test_h5ad(write_test_h5ad(normalized_x)),
        "normalized .*adata\\.layers\\['counts'\\]"
    )
    # Scanpy's `adata.raw = adata` after log1p leaves normalized values in .raw.
    expect_error(
        load_test_h5ad(write_test_h5ad(normalized_x, raw_X = normalized_x)),
        "normalized"
    )
})

test_that("h5ad .raw counts must align with the cells in obs", {
    skip_without_h5ad_writer()
    counts <- h5ad_test_counts()
    normalized_x <- methods::as(h5ad_test_lognorm(counts), "CsparseMatrix")

    expect_error(
        load_test_h5ad(write_test_h5ad(normalized_x, raw_X = counts[1:20, ])),
        "does not match the 30 cells in obs"
    )
})

# Rewrites the sparse .raw/X so each row (CSR) or column (CSC) lists its
# entries twice, at half value and in reverse index order. scipy allows
# both unsorted and duplicate indices; summing restores the counts.
scramble_raw_sparse <- function(path) {
    h5 <- hdf5r::H5File$new(path, mode = "r+")
    on.exit(h5$close_all(), add = TRUE)
    raw_x <- h5[["raw/X"]]
    indptr <- raw_x[["indptr"]]$read()
    indices <- raw_x[["indices"]]$read()
    values <- raw_x[["data"]]$read()

    new_indices <- new_values <- numeric()
    for (k in seq_len(length(indptr) - 1L)) {
        segment <- rev(indptr[k] + seq_len(indptr[k + 1L] - indptr[k]))
        new_indices <- c(new_indices, indices[segment], indices[segment])
        new_values <- c(new_values, values[segment] / 2, values[segment] / 2)
    }
    for (name in c("indptr", "indices", "data")) {
        raw_x$link_delete(name)
    }
    raw_x$create_dataset("indptr", robj = as.integer(indptr * 2))
    raw_x$create_dataset("indices", robj = as.integer(new_indices))
    raw_x$create_dataset("data", robj = new_values)
    path
}

test_that("h5ad .raw sparse counts may have unsorted and duplicate indices", {
    skip_without_h5ad_writer()
    counts <- h5ad_test_counts()
    normalized_x <- methods::as(h5ad_test_lognorm(counts)[, 1:15], "CsparseMatrix")
    raw_inputs <- list(
        csr = methods::as(counts, "RsparseMatrix"),
        csc = methods::as(counts, "CsparseMatrix")
    )

    for (encoding in names(raw_inputs)) {
        path <- scramble_raw_sparse(
            write_test_h5ad(normalized_x, raw_X = raw_inputs[[encoding]])
        )
        obj <- load_test_h5ad(path)
        expect_equal(
            as.matrix(SeuratObject::LayerData(obj, layer = "counts")),
            t(counts),
            info = encoding
        )
    }
})

test_that("h5ad loading accepts gene names with underscores", {
    skip_without_h5ad_writer()
    counts <- h5ad_test_counts()
    colnames(counts) <- paste0("gene_", seq_len(ncol(counts)))
    normalized_x <- methods::as(h5ad_test_lognorm(counts), "CsparseMatrix")
    expected <- t(counts)
    # Seurat replaces underscores in feature names with dashes.
    rownames(expected) <- sub("_", "-", rownames(expected), fixed = TRUE)

    files <- list(
        integer_x = write_test_h5ad(methods::as(counts, "CsparseMatrix")),
        counts_layer = write_test_h5ad(
            normalized_x,
            layers = list(counts = methods::as(counts, "CsparseMatrix"))
        ),
        raw = write_test_h5ad(normalized_x, raw_X = counts)
    )
    for (layout in names(files)) {
        obj <- load_test_h5ad(files[[layout]])
        expect_equal(
            as.matrix(SeuratObject::LayerData(obj, layer = "counts")),
            expected,
            info = layout
        )
        expect_equal(obj@meta.data$sample, rep("s1", nrow(counts)), info = layout)
    }
})
