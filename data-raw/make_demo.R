# data-raw/make_demo.R
# Generates inst/extdata/pbmc_demo.rds, the bundled teaching demo dataset.
#
# The demo is a reproducible 2,000-cell subset of the 10X Genomics PBMC 3k
# filtered feature-barcode matrix used by many Seurat tutorials. It is large
# enough to exercise the workflow while still loading quickly on a laptop.
#
# Run from the repository root:
#   Rscript data-raw/make_demo.R

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
})

set.seed(42)

demo_url <- paste0(
    "https://cf.10xgenomics.com/samples/cell/pbmc3k/",
    "pbmc3k_filtered_gene_bc_matrices.tar.gz"
)

tmpdir <- tempfile("ascseurat_pbmc_demo_")
dir.create(tmpdir, recursive = TRUE)
on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

archive <- file.path(tmpdir, "pbmc3k_filtered_gene_bc_matrices.tar.gz")
message("Downloading PBMC 3k matrix from 10X Genomics...")
utils::download.file(demo_url, archive, mode = "wb", quiet = FALSE)
utils::untar(archive, exdir = tmpdir)

data_dir <- file.path(tmpdir, "filtered_gene_bc_matrices", "hg19")
if (!dir.exists(data_dir)) {
    stop("PBMC data directory not found after extraction: ", data_dir, call. = FALSE)
}

counts <- Read10X(data.dir = data_dir)
selected_cells <- sample(colnames(counts), min(2000, ncol(counts)))

demo <- CreateSeuratObject(
    counts = counts[, selected_cells, drop = FALSE],
    project = "PBMC_Demo",
    min.cells = 3,
    min.features = 200
)
demo[["percent.mt"]] <- PercentageFeatureSet(demo, pattern = "^MT-")

message("Demo dataset created: ", ncol(demo), " cells, ", nrow(demo), " genes")
message(
    "percent.mt range: ",
    round(min(demo$percent.mt), 3), " - ",
    round(max(demo$percent.mt), 3)
)

out_path <- file.path("inst", "extdata", "pbmc_demo.rds")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(demo, out_path, compress = "xz")
message("Saved to ", out_path)
