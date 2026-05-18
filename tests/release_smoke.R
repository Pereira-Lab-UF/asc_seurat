#!/usr/bin/env Rscript

options(warn = 1)
options(repos = c(CRAN = "https://cloud.r-project.org"))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
    sub("^--file=", "", script_arg[[1]])
} else {
    file.path("tests", "release_smoke.R")
}
base_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!file.exists(file.path(base_dir, "DESCRIPTION"))) {
    base_dir <- normalizePath(getwd())
}
setwd(base_dir)

suppressPackageStartupMessages({
    pkgload::load_all(".", quiet = TRUE)
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
})
source(file.path("tests", "testthat", "helper-fixtures.R"))

out_dir <- file.path("tests", "output", "release_smoke")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

results <- data.frame(
    area = character(),
    check = character(),
    status = character(),
    seconds = numeric(),
    detail = character(),
    stringsAsFactors = FALSE
)

record <- function(area, name, status, seconds, detail = "") {
    results <<- rbind(
        results,
        data.frame(
            area = area,
            check = name,
            status = status,
            seconds = round(seconds, 2),
            detail = detail,
            stringsAsFactors = FALSE
        )
    )
}

check <- function(area, name, expr) {
    cat(sprintf("• %-18s %s ... ", paste0("[", area, "]"), name))
    started <- Sys.time()
    detail <- ""
    ok <- tryCatch(
        {
            value <- force(expr)
            if (!is.null(value)) {
                detail <- paste(utils::capture.output(print(value)), collapse = " ")
            }
            TRUE
        },
        error = function(e) {
            detail <<- conditionMessage(e)
            FALSE
        }
    )
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (ok) {
        cat("PASS\n")
        record(area, name, "PASS", elapsed, detail)
    } else {
        cat("FAIL\n")
        cat("  ", detail, "\n", sep = "")
        record(area, name, "FAIL", elapsed, detail)
    }
    invisible(ok)
}

assert <- function(condition, message) {
    if (!isTRUE(condition)) {
        stop(message, call. = FALSE)
    }
}

make_counts <- function(n_genes = 200, n_cells = 80) {
    counts <- Seurat::GetAssayData(demo_fixture(), layer = "counts")
    counts <- counts[seq_len(min(n_genes, nrow(counts))),
                     seq_len(min(n_cells, ncol(counts))), drop = FALSE]
    counts <- counts[Matrix::rowSums(counts) > 0, Matrix::colSums(counts) > 0, drop = FALSE]
    counts
}

write_10x_dir <- function(counts, dir) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    Matrix::writeMM(counts, file.path(dir, "matrix.mtx"))
    writeLines(colnames(counts), file.path(dir, "barcodes.tsv"))
    features <- data.frame(
        id = rownames(counts),
        name = rownames(counts),
        type = "Gene Expression"
    )
    utils::write.table(
        features,
        file.path(dir, "features.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
}

write_h5ad <- function(counts, path) {
    scanpy_ok <- ascseurat:::paga_available()
    assert(scanpy_ok, "Scanpy/AnnData Python stack is unavailable")
    anndata <- reticulate::import("anndata", convert = FALSE)
    pandas <- reticulate::import("pandas", convert = FALSE)
    scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
    adata <- anndata$AnnData(
        X = scipy_sparse$csr_matrix(t(as.matrix(counts))),
        obs = pandas$DataFrame(index = reticulate::r_to_py(colnames(counts))),
        var = pandas$DataFrame(index = reticulate::r_to_py(rownames(counts)))
    )
    adata$write_h5ad(path)
    invisible(path)
}

save_plot_check <- function(plot, filename, width = 8, height = 6) {
    ggplot2::ggsave(filename, plot, width = width, height = height, dpi = 100)
    assert(file.exists(filename) && file.info(filename)$size > 0,
           paste("Plot was not written:", filename))
}

cat("\nAsc-Seurat release smoke test\n")
cat("Output:", normalizePath(out_dir, mustWork = FALSE), "\n\n")

required_packages <- c(
    "shiny", "bslib", "Seurat", "SeuratObject", "scCustomize", "slingshot",
    "harmony", "SingleCellExperiment", "DelayedMatrixStats",
    "PseudotimeDE", "tradeSeq", "scMaSigPro", "monocle3", "reticulate",
    "anndataR", "DropletUtils", "SingleR", "celldex",
    "scDblFinder", "metap"
)
for (pkg in required_packages) {
    check("Dependencies", pkg, {
        assert(requireNamespace(pkg, quietly = TRUE), paste("Missing package:", pkg))
    })
}

check("Dependencies", "PAGA Python stack", {
    assert(ascseurat:::paga_available(), "PAGA Python stack is unavailable")
    modules <- ascseurat:::paga_python_packages()
    missing <- modules[!vapply(modules, reticulate::py_module_available, logical(1))]
    assert(!length(missing), paste("Missing Python modules:", paste(missing, collapse = ", ")))
})

counts <- make_counts()
tmp_formats <- tempfile("ascseurat_formats_")
dir.create(tmp_formats)

check("Loading", "CSV count matrix", {
    path <- file.path(tmp_formats, "counts.csv")
    utils::write.csv(as.matrix(counts), path)
    obj <- ascseurat:::load_input_data(path, project = "csv", min.cells = 0, min.features = 0)
    assert(inherits(obj, "Seurat"), "CSV did not load as Seurat")
    assert(ncol(obj) == ncol(counts), "CSV cell count changed")
})

check("Loading", "TSV count matrix", {
    path <- file.path(tmp_formats, "counts.tsv")
    utils::write.table(as.matrix(counts), path, sep = "\t", quote = FALSE)
    obj <- ascseurat:::load_input_data(path, project = "tsv", min.cells = 0, min.features = 0)
    assert(inherits(obj, "Seurat"), "TSV did not load as Seurat")
    assert(nrow(obj) == nrow(counts), "TSV gene count changed")
})

check("Loading", "10X directory", {
    path <- file.path(tmp_formats, "tenx_dir")
    write_10x_dir(counts, path)
    assert(identical(ascseurat:::detect_input_format(path), "10x_dir"), "10X directory not detected")
    obj <- ascseurat:::load_input_data(path, project = "tenx", min.cells = 0, min.features = 0)
    assert(inherits(obj, "Seurat"), "10X directory did not load")
})

check("Loading", "10X H5", {
    path <- file.path(tmp_formats, "counts.h5")
    DropletUtils::write10xCounts(path, counts, overwrite = TRUE, type = "HDF5", version = "3")
    obj <- ascseurat:::load_input_data(path, project = "h5", min.cells = 0, min.features = 0)
    assert(inherits(obj, "Seurat"), "10X H5 did not load")
    assert("counts" %in% SeuratObject::Layers(obj), "H5 counts layer missing")
})

check("Loading", "H5AD", {
    path <- file.path(tmp_formats, "counts.h5ad")
    write_h5ad(counts, path)
    obj <- suppressWarnings(ascseurat:::load_input_data(path, project = "h5ad", min.cells = 0, min.features = 0))
    assert(inherits(obj, "Seurat"), "H5AD did not load")
    assert("counts" %in% SeuratObject::Layers(obj), "H5AD counts layer missing")
})

check("Loading", "RDS round trip", {
    path <- file.path(tmp_formats, "demo.rds")
    saveRDS(demo_fixture(), path)
    obj <- ascseurat:::load_input_data(path, project = "rds", min.cells = 0, min.features = 0)
    assert(inherits(obj, "Seurat"), "RDS did not load")
    assert(ncol(obj) == ncol(demo_fixture()), "RDS cell count changed")
})

check("QC", "MAD defaults", {
    defaults <- ascseurat:::qc_filter_defaults(demo_fixture()@meta.data)
    assert(defaults$min_genes >= 0, "Invalid min gene default")
    assert(defaults$max_genes > defaults$min_genes, "Invalid max gene default")
    assert(defaults$max_mito >= 0 && defaults$max_mito <= 100, "Invalid mitochondrial default")
})

check("QC", "No filtering path", {
    obj <- demo_fixture()
    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(session$setInputs(filter_mode = "none"))
            out <- qc_result()
            assert(ncol(out) == ncol(obj), "No filtering changed cell count")
            assert(identical(attr(out, "ascseurat_filter_mode", exact = TRUE), "none"),
                   "No filtering mode attr missing")
        }
    )
})

check("QC", "Threshold filtering path", {
    obj <- demo_fixture()
    shiny::testServer(
        ascseurat::mod_qc_server,
        args = list(seurat_obj = shiny::reactive(obj)),
        {
            suppressWarnings(session$setInputs(
                filter_mode = "filter",
                min_genes = 0,
                max_genes = max(obj$nFeature_RNA, na.rm = TRUE) + 1,
                max_mito = 100,
                detect_doublets = FALSE,
                apply_filter = 1
            ))
            assert(inherits(qc_result(), "Seurat"), "Filtering did not return a Seurat object")
        }
    )
})

check("QC", "Doublet detection path", {
    obj <- demo_fixture()[, seq_len(min(300, ncol(demo_fixture())))]
    obj <- ascseurat:::detect_doublets(obj)
    assert("scDblFinder_class" %in% colnames(obj@meta.data), "Doublet class missing")
})

check("Normalization", "LogNormalize", {
    obj <- ascseurat:::run_normalization(demo_fixture(), method = "lognorm", nfeatures = 500)
    assert("pca" %in% names(obj@reductions), "PCA missing after LogNormalize")
})

check("Normalization", "SCTransform v2", {
    obj <- ascseurat:::run_normalization(demo_fixture()[, 1:400], method = "sctransform", nfeatures = 500)
    assert("SCT" %in% names(obj@assays), "SCT assay missing")
    assert("pca" %in% names(obj@reductions), "PCA missing after SCTransform")
})

check("Normalization", "Human cell-cycle scoring/regression", {
    obj <- ascseurat:::run_normalization(
        demo_fixture()[, 1:400],
        method = "lognorm",
        nfeatures = 500,
        score_cell_cycle = TRUE,
        cc_organism = "human",
        regress_cc = TRUE
    )
    assert(all(c("S.Score", "G2M.Score", "Phase") %in% colnames(obj@meta.data)),
           "Cell-cycle metadata missing")
})

check("Normalization", "Mouse cell-cycle scoring", {
    genes <- unique(c(
        paste0(toupper(substr(tolower(Seurat::cc.genes.updated.2019$s.genes), 1, 1)),
               substr(tolower(Seurat::cc.genes.updated.2019$s.genes), 2, 100)),
        paste0(toupper(substr(tolower(Seurat::cc.genes.updated.2019$g2m.genes), 1, 1)),
               substr(tolower(Seurat::cc.genes.updated.2019$g2m.genes), 2, 100))
    ))
    genes <- head(genes, 80)
    mat <- Matrix::Matrix(sample(0:6, length(genes) * 80, replace = TRUE), nrow = length(genes), sparse = TRUE)
    rownames(mat) <- genes
    colnames(mat) <- paste0("cell", seq_len(ncol(mat)))
    obj <- Seurat::CreateSeuratObject(mat, min.cells = 0, min.features = 0)
    obj <- ascseurat:::run_normalization(obj, method = "lognorm", nfeatures = 50,
                                         score_cell_cycle = TRUE, cc_organism = "mouse")
    assert("Phase" %in% colnames(obj@meta.data), "Mouse cell-cycle phase missing")
})

clustered <- NULL
check("Clustering", "UMAP and tSNE", {
    obj <- ascseurat:::run_normalization(demo_fixture(), method = "lognorm", nfeatures = 500)
    clustered <<- ascseurat:::run_clustering(obj, dims = 1:10, resolution = 0.4, run_tsne = TRUE)
    assert(all(c("umap", "tsne") %in% names(clustered@reductions)), "UMAP/tSNE missing")
    assert(length(levels(SeuratObject::Idents(clustered))) > 1, "Only one cluster found")
})

check("Clustering", "Module uses selected upstream PCs", {
    shiny::testServer(
        ascseurat::mod_clustering_server,
        args = list(norm_result = list(
            obj = shiny::reactive(ascseurat:::run_normalization(demo_fixture(), method = "lognorm", nfeatures = 500)),
            n_pcs = shiny::reactive(10)
        )),
        {
            session$setInputs(resolution = 0.4, run_clustering = 1)
            assert(identical(attr(active_obj(), "ascseurat_clustering_n_pcs", exact = TRUE), 10L),
                   "Clustering did not use upstream PC count")
            session$setInputs(rename_1 = "Renamed cluster", apply_rename = 1)
            assert(any(grepl("Renamed cluster", levels(SeuratObject::Idents(active_obj())))),
                   "Cluster rename failed")
        }
    )
})

make_merged <- function() {
    obj <- demo_fixture()[1:800, 1:240]
    obj1 <- obj[, 1:120]
    obj2 <- obj[, 121:240]
    obj1$orig.ident <- "sample_1"
    obj2$orig.ident <- "sample_2"
    merge(obj1, obj2, add.cell.ids = c("s1", "s2"))
}

for (method in c("rpca", "harmony")) {
    for (norm_method in c("lognorm", "sctransform")) {
        check("Integration", paste(method, norm_method), {
            obj <- ascseurat:::run_integration(
                make_merged(),
                method = method,
                normalization_method = norm_method,
                dims = 1:10,
                nfeatures = 500
            )
            reduction <- attr(obj, "ascseurat_default_reduction", exact = TRUE)
            assert(reduction %in% names(obj@reductions), "Integrated reduction missing")
            assert("umap" %in% names(obj@reductions), "Integrated UMAP missing")
        })
    }
}

for (de_type in c("all", "one", "compare", "conditions", "conserved")) {
    check("DE", de_type, {
        obj <- clustered
        clusters <- levels(SeuratObject::Idents(obj))
        obj$samples <- rep(c("A", "B"), length.out = ncol(obj))
        cluster <- clusters[[1]]
        cluster_cells <- names(SeuratObject::Idents(obj))[as.character(SeuratObject::Idents(obj)) == cluster]
        obj$condition <- "other"
        obj$condition[cluster_cells] <- rep(c("A", "B"), length.out = length(cluster_cells))
        shiny::testServer(
            ascseurat::mod_de_server,
            args = list(seurat_obj = shiny::reactive(obj), allow_conserved = TRUE),
            {
                common <- list(
                    test_use = "wilcox",
                    pval_cutoff = 0.999999,
                    logfc_threshold = 0,
                    min_pct = 0,
                    only_pos = FALSE,
                    run_de = 1
                )
                inputs <- switch(
                    de_type,
                    all = c(list(de_type = "all"), common),
                    one = c(list(de_type = "one", cluster_1 = cluster), common),
                    compare = c(list(de_type = "compare", cluster_1 = cluster, cluster_2 = clusters[[2]]), common),
                    conditions = c(list(
                        de_type = "conditions",
                        cluster_1 = cluster,
                        condition_column = "condition",
                        condition_value_1 = "A",
                        condition_value_2 = "B"
                    ), common),
                    conserved = c(list(
                        de_type = "conserved",
                        cluster_1 = cluster,
                        conserved_grouping = "samples"
                    ), common)
                )
                do.call(session$setInputs, inputs)
                assert(is.data.frame(markers()), "DE result is not a data frame")
                assert("gene" %in% colnames(markers()), "DE gene column missing")
            }
        )
    })
}

check("Visualization", "Upload markers, grouped genes, all expression slots", {
    genes <- head(rownames(clustered), 4)
    marker_file <- file.path(out_dir, "markers_grouped.csv")
    utils::write.csv(data.frame(gene = genes, group = c("A", "A", "B", "B")), marker_file, row.names = FALSE)
    shiny::testServer(
        ascseurat::mod_visualization_server,
        args = list(
            seurat_obj = shiny::reactive(clustered),
            markers_de = shiny::reactive(data.frame(
                gene = genes,
                cluster = c("0", "0", "1", "1"),
                p_val_adj = c(0.01, 0.02, 0.03, 0.04)
            ))
        ),
        {
            session$setInputs(
                marker_source = "upload",
                markers_file = list(datapath = marker_file),
                header_opt = "Yes",
                load_markers = 1,
                selected_group = "A",
                upload_gene_count = 2,
                selected_genes = genes[1:2],
                show_expression = 1
            )
            for (slot in c("counts", "data", "scale.data")) {
                session$setInputs(slot_selection = slot)
                assert(!is.null(clustered_dot_plot()), paste("Clustered dot plot failed for", slot))
            }
            assert(!is.null(feature_plot()), "Feature plot failed")
            assert(!is.null(violin_plot()), "Violin plot failed")
            session$setInputs(marker_source = "de_results", de_cluster = "0", top_n_genes = 2)
            assert(length(active_genes()) == 2, "DE result marker source failed")
        }
    )
})

for (plot_type in c("stacked_vln", "dotplot")) {
    check("Advanced Plots", plot_type, {
        genes <- head(rownames(clustered), 4)
        marker_file <- file.path(out_dir, paste0("advanced_", plot_type, ".csv"))
        utils::write.csv(data.frame(gene = genes), marker_file, row.names = FALSE)
        shiny::testServer(
            ascseurat::mod_advanced_plots_server,
            args = list(
                seurat_obj_single = shiny::reactive(clustered),
                seurat_obj_multi = shiny::reactive(clustered)
            ),
            {
                session$setInputs(
                    data_source = "single",
                    plot_type = plot_type,
                    gene_list = list(datapath = marker_file),
                    make_plot = 1
                )
                assert(!is.null(adv_plot()), "Advanced plot did not render")
            }
        )
    })
}

for (method in c("slingshot", "paga", "monocle3")) {
    check("Trajectory", method, {
        res <- switch(
            method,
            slingshot = ascseurat:::run_slingshot_trajectory(clustered, start_cluster = "0", end_cluster = NULL),
            paga = ascseurat:::run_paga_trajectory(clustered, start_cluster = "0"),
            monocle3 = ascseurat:::run_monocle3_trajectory(clustered, start_cluster = "0")
        )
        assert(all(c("trajectory_method", "pseudotime", "trajectory_path") %in%
            colnames(res$obj@meta.data)), "Standard trajectory metadata missing")
        assert(any(is.finite(res$obj$pseudotime)), "Pseudotime contains no finite values")
    })
}

check("Trajectory DE", "PseudotimeDE-fast", {
    res <- ascseurat:::run_slingshot_trajectory(clustered, start_cluster = "0")
    de <- ascseurat:::run_pseudotimede_fast(res$sce, gene_subset = head(rownames(clustered), 8), cores = 1)
    assert(nrow(de$table) > 0, "PseudotimeDE returned no rows")
})

check("Trajectory DE", "tradeSeq", {
    res <- ascseurat:::run_slingshot_trajectory(clustered, start_cluster = "0")
    de <- ascseurat:::run_tradeseq(res$sce, res$assignment, gene_subset = head(rownames(clustered), 8), nknots = 3)
    assert(nrow(de$table) > 0, "tradeSeq returned no rows")
})

check("Trajectory DE", "scMaSigPro", {
    res <- ascseurat:::run_slingshot_trajectory(clustered, start_cluster = "0")
    de <- ascseurat:::run_scmasigpro(
        res$sce,
        gene_subset = head(rownames(clustered), 8),
        p_value = 0.99,
        rsq = 0,
        poly_degree = 2
    )
    assert(nrow(de$table) > 0, "scMaSigPro returned no rows")
})

check("SingleR", "All references and granularities", {
    obj <- clustered[, seq_len(min(120, ncol(clustered)))]
    expr <- Seurat::GetAssayData(obj, layer = "data")
    refs <- c("hpca", "blueprint", "monaco", "mousernaseq", "immgen")
    cache_dir <- file.path(tools::R_user_dir("ascseurat", "cache"), "celldex")
    for (ref_name in refs) {
        ref <- readRDS(file.path(cache_dir, paste0(ref_name, ".rds")))
        for (granularity in c("label.main", "label.fine")) {
            labels <- ref[[granularity]]
            if (all(is.na(labels))) {
                labels <- ref$label.main
            }
            pred <- SingleR::SingleR(test = expr, ref = ref, labels = labels)
            assert(length(pred$labels) == ncol(obj),
                   paste("SingleR label count mismatch:", ref_name, granularity))
        }
    }
})

check("UI", "Collapsed download options", {
    html <- htmltools::renderTags(ascseurat:::card_plot_download("x", "download_x"))$html
    assert(grepl("<details", html, fixed = TRUE), "Download options are not collapsible details")
    assert(!grepl("<details open", html, fixed = TRUE), "Download options are expanded by default")
})

check("UI", "Busy RDS download button", {
    html <- htmltools::renderTags(ascseurat:::busy_download_button("download_rds", "Save Seurat Object (RDS)"))$html
    assert(grepl("asc-busy-download", html, fixed = TRUE), "Busy download class missing")
    assert(grepl("Preparing download", html, fixed = TRUE), "Busy download label missing")
})

write.csv(results, file.path(out_dir, "release_smoke_results.csv"), row.names = FALSE)

cat("\n")
print(table(results$status))

failed <- results[results$status != "PASS", , drop = FALSE]
if (nrow(failed)) {
    cat("\nFailed checks:\n")
    print(failed[, c("area", "check", "detail")], row.names = FALSE)
    quit(save = "no", status = 1)
}

cat("\nAll release smoke checks passed.\n")
