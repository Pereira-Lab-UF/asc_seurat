#!/usr/bin/env Rscript
# ============================================================================
# Asc-Seurat v3 — Non-Interactive Pipeline Test
# ============================================================================
# Tests all core utility functions (not Shiny modules) using the Arabidopsis
# dataset from Ryu et al. 2019.
#
# Usage: Rscript tests/test_v3_pipeline.R
# ============================================================================

cat("\n====================================================================\n")
cat("  ASC-SEURAT v3 — PIPELINE TEST\n")
cat("====================================================================\n\n")

if (!identical(Sys.getenv("ASCSEURAT_RUN_INTEGRATION_TESTS"), "1")) {
    cat("Skipping slow integration test. Set ASCSEURAT_RUN_INTEGRATION_TESTS=1 to run it.\n")
    quit(save = "no", status = 0)
}

# ─── Setup ──────────────────────────────────────────────────────────────────
cat("PHASE 0: Loading packages...\n")

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(ggplot2)
    library(dplyr)
    library(slingshot)
    library(SingleCellExperiment)
})

# Check for scCustomize
has_scCustomize <- requireNamespace("scCustomize", quietly = TRUE)
if (has_scCustomize) {
    suppressPackageStartupMessages(library(scCustomize))
    cat("  ✓ scCustomize loaded\n")
} else {
    cat("  ⚠ scCustomize not available — will use base Seurat plots\n")
}

# Check for harmony
has_harmony <- requireNamespace("harmony", quietly = TRUE)
if (has_harmony) cat("  ✓ harmony loaded\n")

# Source utility files
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
    sub("^--file=", "", script_arg[[1]])
} else if (!is.null(sys.frame(1)$ofile)) {
    sys.frame(1)$ofile
} else {
    file.path(getwd(), "tests", "test_v3_pipeline.R")
}

base_dir <- normalizePath(file.path(dirname(script_path), ".."))
cat("  Base directory:", base_dir, "\n")
source(file.path(base_dir, "R", "utils_compat.R"))
source(file.path(base_dir, "R", "utils_seurat.R"))
source(file.path(base_dir, "R", "utils_plots.R"))

# Output directory
out_dir <- file.path(base_dir, "tests", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Data paths
wt_path <- file.path(base_dir, "arabidopsis_mutant_scRNAseq",
                     "Ryu_et_al_2019_PMID_30718350_Arabidopsis_WT",
                     "filtered_gene_bc_matrices", "TAIR10")
rhd6_path <- file.path(base_dir, "arabidopsis_mutant_scRNAseq",
                       "Ryu_et_al_2019_PMID_30718350_Arabidopsis_rhd6",
                       "filtered_gene_bc_matrices", "TAIR10")

# Track test results
passed <- 0
failed <- 0
errors <- c()

test <- function(name, expr) {
    cat(sprintf("\n  TEST: %s...", name))
    tryCatch({
        eval(expr)
        cat(" ✓ PASS\n")
        passed <<- passed + 1
    }, error = function(e) {
        cat(sprintf(" ✗ FAIL: %s\n", conditionMessage(e)))
        failed <<- failed + 1
        errors <<- c(errors, paste0(name, ": ", conditionMessage(e)))
    })
}


# ─── PHASE 1: Data Loading ──────────────────────────────────────────────────
cat("\n--------------------------------------------------------------------\n")
cat("PHASE 1: Data Loading & Format Detection\n")
cat("--------------------------------------------------------------------\n")

test("detect_input_format() — 10X directory", {
    fmt <- detect_input_format(wt_path)
    stopifnot(fmt == "10x_dir")
})

test("load_input_data() — WT sample (10X)", {
    wt <- load_input_data(wt_path, project = "WT", min.cells = 3, min.features = 200)
    stopifnot(inherits(wt, "Seurat"))
    stopifnot(ncol(wt) > 0)
    cat(sprintf(" [%d cells, %d genes]", ncol(wt), nrow(wt)))
})

wt <- load_input_data(wt_path, project = "WT", min.cells = 3, min.features = 200)

test("load_input_data() — rhd6 sample (10X)", {
    rhd6 <- load_input_data(rhd6_path, project = "rhd6", min.cells = 3, min.features = 200)
    stopifnot(inherits(rhd6, "Seurat"))
    cat(sprintf(" [%d cells, %d genes]", ncol(rhd6), nrow(rhd6)))
})

rhd6 <- load_input_data(rhd6_path, project = "rhd6", min.cells = 3, min.features = 200)

# Test RDS save/load cycle
test("RDS save/load round-trip", {
    tmp_rds <- file.path(out_dir, "test_wt.rds")
    saveRDS(wt, tmp_rds)
    fmt <- detect_input_format(tmp_rds)
    stopifnot(fmt == "rds")
    wt_reloaded <- load_seurat_rds(tmp_rds)
    stopifnot(inherits(wt_reloaded, "Seurat"))
    stopifnot(ncol(wt_reloaded) == ncol(wt))
    unlink(tmp_rds)
})


# ─── PHASE 2: QC & Filtering ────────────────────────────────────────────────
cat("\n--------------------------------------------------------------------\n")
cat("PHASE 2: QC & Filtering\n")
cat("--------------------------------------------------------------------\n")

# Arabidopsis mitochondria: genes starting with "ATMG" or on the Mt chromosome
# Also try standard "^MT-" which won't match — test graceful handling
test("PercentageFeatureSet() — mito detection", {
    wt[["percent.mt"]] <- PercentageFeatureSet(wt, pattern = "^ATMG")
    stopifnot("percent.mt" %in% colnames(wt@meta.data))
    cat(sprintf(" [median mito: %.2f%%]", median(wt$percent.mt)))
})

test("run_qc_filter()", {
    n_before <- ncol(wt)
    wt_filtered <- run_qc_filter(wt, mito_pattern = "^ATMG",
                                  min_features = 200, max_features = 5000,
                                  max_mito_pct = 20)
    n_after <- ncol(wt_filtered)
    stopifnot(n_after > 0)
    stopifnot(n_after <= n_before)
    cat(sprintf(" [%d → %d cells]", n_before, n_after))
})

wt_filtered <- run_qc_filter(wt, mito_pattern = "^ATMG",
                              min_features = 200, max_features = 5000,
                              max_mito_pct = 20)

if (has_scCustomize) {
    test("QC Violin plots (scCustomize)", {
        p <- VlnPlot_scCustom(wt_filtered,
                              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                              pt.size = 0.1)
        ggsave(file.path(out_dir, "qc_violin.png"), p, width = 12, height = 5, dpi = 150)
    })
}


# ─── PHASE 3: Normalization + PCA ───────────────────────────────────────────
cat("\n--------------------------------------------------------------------\n")
cat("PHASE 3: Normalization & PCA\n")
cat("--------------------------------------------------------------------\n")

test("run_normalization() — LogNormalize", {
    wt_norm <- run_normalization(wt_filtered, method = "lognorm", nfeatures = 2000)
    stopifnot(inherits(wt_norm, "Seurat"))
    stopifnot("pca" %in% names(wt_norm@reductions))
    cat(sprintf(" [%d variable features]", length(VariableFeatures(wt_norm))))
})

wt_norm <- run_normalization(wt_filtered, method = "lognorm", nfeatures = 2000)

test("ElbowPlot", {
    p <- ElbowPlot(wt_norm, ndims = 50)
    ggsave(file.path(out_dir, "elbow_plot.png"), p, width = 8, height = 5, dpi = 150)
})

test("run_normalization() — SCTransform v2", {
    wt_sct <- run_normalization(wt_filtered, method = "sctransform", nfeatures = 2000)
    stopifnot(inherits(wt_sct, "Seurat"))
    stopifnot("pca" %in% names(wt_sct@reductions))
})


# ─── PHASE 4: Clustering + UMAP ─────────────────────────────────────────────
cat("\n--------------------------------------------------------------------\n")
cat("PHASE 4: Clustering & UMAP\n")
cat("--------------------------------------------------------------------\n")

test("run_clustering()", {
    wt_clust <- run_clustering(wt_norm, dims = 1:20, resolution = 0.6)
    stopifnot("umap" %in% names(wt_clust@reductions))
    n_clust <- length(levels(Idents(wt_clust)))
    stopifnot(n_clust > 1)
    cat(sprintf(" [%d clusters]", n_clust))
})

wt_clust <- run_clustering(wt_norm, dims = 1:20, resolution = 0.6)

if (has_scCustomize) {
    test("DimPlot_scCustom (figure_plot=TRUE)", {
        p <- DimPlot_scCustom(wt_clust, figure_plot = TRUE, label = TRUE)
        ggsave(file.path(out_dir, "umap_clusters.png"), p, width = 8, height = 7, dpi = 150)
    })

    test("FeaturePlot_scCustom", {
        # Pick a gene that's likely expressed
        genes <- head(VariableFeatures(wt_clust), 4)
        p <- FeaturePlot_scCustom(wt_clust, features = genes, num_columns = 2)
        ggsave(file.path(out_dir, "feature_plots.png"), p, width = 10, height = 10, dpi = 150)
    })
}


# ─── PHASE 5: Differential Expression ───────────────────────────────────────
cat("\n--------------------------------------------------------------------\n")
cat("PHASE 5: Differential Expression\n")
cat("--------------------------------------------------------------------\n")

test("FindAllMarkers()", {
    markers <- FindAllMarkers(wt_clust, only.pos = TRUE, logfc.threshold = 0.25,
                              min.pct = 0.1, verbose = FALSE)
    stopifnot(nrow(markers) > 0)
    stopifnot(all(c("gene", "cluster", "p_val_adj") %in% colnames(markers)))
    cat(sprintf(" [%d markers across %d clusters]",
        nrow(markers), length(unique(markers$cluster))))
})

markers <- FindAllMarkers(wt_clust, only.pos = TRUE, logfc.threshold = 0.25,
                          min.pct = 0.1, verbose = FALSE)

test("FindMarkers() — cluster 0 vs rest", {
    m <- FindMarkers(wt_clust, ident.1 = "0", verbose = FALSE)
    stopifnot(nrow(m) > 0)
    cat(sprintf(" [%d DE genes]", nrow(m)))
})

test("Heatmap (DoHeatmap)", {
    top_genes <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5) %>%
        pull(gene) %>% unique()
    p <- DoHeatmap(wt_clust, features = top_genes) + scale_fill_viridis_c()
    ggsave(file.path(out_dir, "heatmap.png"), p, width = 12, height = 10, dpi = 150)
})


# ─── PHASE 6: Trajectory (Slingshot) ────────────────────────────────────────
cat("\n--------------------------------------------------------------------\n")
cat("PHASE 6: Trajectory Inference (Slingshot)\n")
cat("--------------------------------------------------------------------\n")

test("Seurat → SingleCellExperiment conversion", {
    sce <- as.SingleCellExperiment(wt_clust)
    stopifnot(inherits(sce, "SingleCellExperiment"))
})

test("slingshot()", {
    sce <- as.SingleCellExperiment(wt_clust)
    sce <- slingshot(sce, clusterLabels = "ident", reducedDim = "UMAP")
    pt <- slingPseudotime(sce)
    stopifnot(!is.null(pt))
    stopifnot(nrow(pt) == ncol(wt_clust))
    n_lineages <- ncol(pt)
    cat(sprintf(" [%d lineages]", n_lineages))
})

sce <- as.SingleCellExperiment(wt_clust)
sce <- slingshot(sce, clusterLabels = "ident", reducedDim = "UMAP")
pt <- slingPseudotime(sce)
wt_clust$pseudotime <- pt[, 1]

test("Pseudotime UMAP", {
    if (has_scCustomize) {
        p <- FeaturePlot_scCustom(wt_clust, features = "pseudotime", order = TRUE)
    } else {
        p <- FeaturePlot(wt_clust, features = "pseudotime", order = TRUE)
    }
    ggsave(file.path(out_dir, "pseudotime_umap.png"), p, width = 8, height = 7, dpi = 150)
})

test("Slingshot lineage curves plot", {
    emb <- Embeddings(wt_clust, "umap")
    sds <- SlingshotDataSet(sce)
    png(file.path(out_dir, "slingshot_curves.png"), width = 800, height = 700, res = 150)
    plot(emb, col = as.numeric(as.factor(Idents(wt_clust))),
         pch = 16, cex = 0.4, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
         main = "Slingshot Lineages")
    lines(sds, lwd = 2, col = "black")
    dev.off()
})

test("Trajectory DE — correlation-based", {
    expr <- LayerData(wt_clust, layer = "data")
    test_genes <- head(VariableFeatures(wt_clust), 200)
    cors <- sapply(test_genes, function(g) {
        tryCatch(
            cor(as.numeric(expr[g, ]), wt_clust$pseudotime,
                use = "complete.obs", method = "spearman"),
            error = function(e) NA
        )
    })
    result <- data.frame(gene = names(cors), correlation = cors, abs_cor = abs(cors))
    result <- result[!is.na(result$correlation), ]
    result <- result[order(-result$abs_cor), ]
    stopifnot(nrow(result) > 0)
    cat(sprintf(" [top gene: %s (r=%.3f)]", result$gene[1], result$correlation[1]))
    write.csv(result, file.path(out_dir, "trajectory_de_genes.csv"), row.names = FALSE)
})


# ─── PHASE 7: Integration ───────────────────────────────────────────────────
cat("\n--------------------------------------------------------------------\n")
cat("PHASE 7: Multi-Sample Integration\n")
cat("--------------------------------------------------------------------\n")

test("Load + merge WT + rhd6", {
    rhd6 <- load_input_data(rhd6_path, project = "rhd6", min.cells = 3, min.features = 200)
    rhd6[["percent.mt"]] <- PercentageFeatureSet(rhd6, pattern = "^ATMG")
    rhd6 <- subset(rhd6, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

    merged <- merge(wt_filtered, y = rhd6, add.cell.ids = c("WT", "rhd6"))
    stopifnot(ncol(merged) == ncol(wt_filtered) + ncol(rhd6))
    cat(sprintf(" [%d total cells]", ncol(merged)))
})

rhd6 <- load_input_data(rhd6_path, project = "rhd6", min.cells = 3, min.features = 200)
rhd6[["percent.mt"]] <- PercentageFeatureSet(rhd6, pattern = "^ATMG")
rhd6 <- subset(rhd6, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
merged <- merge(wt_filtered, y = rhd6, add.cell.ids = c("WT", "rhd6"))

test("run_normalization() on merged object", {
    merged_norm <- run_normalization(merged, method = "lognorm", nfeatures = 2000)
    stopifnot("pca" %in% names(merged_norm@reductions))
})

merged_norm <- run_normalization(merged, method = "lognorm", nfeatures = 2000)

if (has_harmony) {
    test("run_integration() — Harmony", {
        merged_int <- run_integration(merged, method = "harmony", dims = 1:20)
        stopifnot("umap" %in% names(merged_int@reductions))
        cat(sprintf(" [%d cells integrated]", ncol(merged_int)))
    })

    merged_int <- run_integration(merged, method = "harmony", dims = 1:20)

    test("Clustering on integrated data", {
        merged_clust <- FindClusters(merged_int, resolution = 0.6, verbose = FALSE)
        n_clust <- length(levels(Idents(merged_clust)))
        cat(sprintf(" [%d clusters]", n_clust))
    })

    merged_clust <- FindClusters(merged_int, resolution = 0.6, verbose = FALSE)

    if (has_scCustomize) {
        test("Integration UMAP — by cluster", {
            p <- DimPlot_scCustom(merged_clust, figure_plot = TRUE, label = TRUE)
            ggsave(file.path(out_dir, "integrated_umap_clusters.png"), p,
                   width = 8, height = 7, dpi = 150)
        })

        test("Integration UMAP — by sample (split)", {
            p <- DimPlot_scCustom(merged_clust, group.by = "orig.ident",
                                  figure_plot = TRUE)
            ggsave(file.path(out_dir, "integrated_umap_samples.png"), p,
                   width = 8, height = 7, dpi = 150)
        })
    }
} else {
    cat("  ⚠ Harmony not available — skipping integration tests\n")
}


# ─── PHASE 8: Advanced Plots ────────────────────────────────────────────────
cat("\n--------------------------------------------------------------------\n")
cat("PHASE 9: Advanced Plots\n")
cat("--------------------------------------------------------------------\n")

if (has_scCustomize) {
    test("Stacked_VlnPlot", {
        top_genes <- head(unique(markers$gene), 8)
        p <- Stacked_VlnPlot(wt_clust, features = top_genes, x_lab_rotate = TRUE)
        ggsave(file.path(out_dir, "stacked_violin.png"), p,
               width = 10, height = 12, dpi = 150)
    })

    test("DotPlot_scCustom", {
        top_genes <- head(unique(markers$gene), 10)
        p <- DotPlot_scCustom(wt_clust, features = top_genes, flip_axes = TRUE)
        ggsave(file.path(out_dir, "dotplot.png"), p, width = 10, height = 6, dpi = 150)
    })

    test("VlnPlot_scCustom — single gene", {
        gene <- VariableFeatures(wt_clust)[1]
        p <- VlnPlot_scCustom(wt_clust, features = gene, pt.size = 0.1)
        ggsave(file.path(out_dir, "violin_single.png"), p, width = 8, height = 5, dpi = 150)
    })
}


# ─── PHASE 10: save_plot() utility ──────────────────────────────────────────
cat("\n--------------------------------------------------------------------\n")
cat("PHASE 10: Utility Functions\n")
cat("--------------------------------------------------------------------\n")

test("save_plot() — multiple formats", {
    p <- ElbowPlot(wt_norm, ndims = 30)
    for (fmt in c("png", "pdf")) {
        save_plot(p, "test_elbow", format = fmt, path = out_dir)
        f <- file.path(out_dir, paste0("test_elbow.", fmt))
        stopifnot(file.exists(f))
    }
})

test("RDS backward compatibility (save + reload)", {
    rds_path <- file.path(out_dir, "test_clustered.rds")
    saveRDS(wt_clust, rds_path)
    reloaded <- load_seurat_rds(rds_path)
    stopifnot(ncol(reloaded) == ncol(wt_clust))
    stopifnot(length(levels(Idents(reloaded))) == length(levels(Idents(wt_clust))))
    unlink(rds_path)
})


# ─── RESULTS ─────────────────────────────────────────────────────────────────
cat("\n====================================================================\n")
cat(sprintf("  RESULTS: %d PASSED, %d FAILED\n", passed, failed))
cat("====================================================================\n")

if (failed > 0) {
    cat("\nFailed tests:\n")
    for (e in errors) {
        cat("  ✗", e, "\n")
    }
}

# List output files
out_files <- list.files(out_dir, full.names = FALSE)
if (length(out_files) > 0) {
    cat("\nOutput files in tests/output/:\n")
    for (f in out_files) {
        size <- file.info(file.path(out_dir, f))$size
        cat(sprintf("  → %s (%s)\n", f, format(size, big.mark = ",")))
    }
}

cat("\nDone!\n")
