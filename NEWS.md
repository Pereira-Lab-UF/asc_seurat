# ascseurat 3.0.9

## Install and distribution

* The Docker image builds again: Monocle 3's dependencies are installed
  explicitly before it is built (the 3.0.8 image was never published).

## Bug fixes: the session no longer closes

* Blank or invalid inputs no longer close the session: a cleared Max
  mitochondrial %, p-value, minimum fraction, or Log2FC field, an empty or
  negative Top N, an invalid mitochondrial regex, metadata with duplicate or
  missing barcodes, and hidden normalization options (cell-cycle regression
  after unticking Score cell cycle) now show a message instead.
* Reanalyzing or re-clustering a very small subset (fewer than 10 cells, or
  every cluster excluded) is refused with a message.
* A 0-byte marker file shows a warning instead of closing the session, and
  one-gene dot plots now render.

## Bug fixes: results

* Conserved markers are filtered on the largest per-group adjusted p-value
  (new column `max_p_val_adj`), use the selected test (`roc` is refused), and
  explain when fewer than two groups can be tested. DE between conditions can
  use any suitable metadata column.
* Harmony uses "Number of PCs for integration", and the clustering PC input
  is capped to the integrated reduction.
* Recomputing a subset of an integrated object repeats the batch correction
  (RPCA or Harmony) as far as the subset allows; samples too small to be
  their own batch are grouped with the largest one, and the subset falls back
  to plain normalization with a warning when integration is not possible.
* Recomputing a subset repeats the Step 3 normalization settings (method,
  variable features, scale factor, regressed variables, cell-cycle scoring),
  and keeps no duplicate RNA layers.
* Scaled expression plots keep the values of genes that were already scaled,
  regress the same variables as normalization, and rescale all genes after a
  subset that kept the current clusters. After SCTransform they no longer
  show RNA scaled values carried in from an uploaded object.
* Advanced Plots use normalized values for uploaded SCTransform objects whose
  RNA assay was never normalized.
* h5ad files: raw counts are read from `layers['counts']`,
  `layers['raw_counts']`, `.raw`, or an integer `X`; normalized-only files are
  rejected with an explanation, and gene names with underscores load.
* Trajectory: Monocle 3 orders cells from a single root node, PAGA uses the
  clustering reduction (Harmony or RPCA after integration), tradeSeq uses
  whole-transcriptome library sizes and a fixed seed, PseudotimeDE-fast tests
  each Slingshot lineage separately, gene tables are ranked consistently,
  `pseudotime_PathN` is `NA` off-lineage, scMaSigPro failures are reported as
  errors, and the Monocle 3 graph works on Integration-tab objects.
* Re-normalizing an integrated `.rds` drops its stale Harmony/RPCA reduction.
* tSNE runs only when its box is ticked and the object has one sample.

## Bug fixes: inputs, stale results, and downloads

* Results that no longer match the current object are cleared: QC,
  normalization, clustering, DE tables and downloads, expression plots,
  Advanced Plots, SingleR results, and trajectory DE, heatmaps and gene
  tables. Re-runs that return the same object keep them.
* A failed Load & Integrate keeps the previous integration; the Integration
  Summary reports the method and PCs actually used, and Number of PCs for
  integration must be a whole number from 5 to 100.
* Integration-tab sample names, paths, and QC thresholds keep typed values
  when samples are renamed or added, QC suggestions use each sample's Mito
  pattern, and restored bookmarks keep their thresholds. The clustering PC
  input is updated after the first integration.
* The Session Report lists only steps that completed, and now records demo
  cluster renaming and subset reanalysis and the post-integration filter.

# ascseurat 3.0.8

## Bug fixes

* Downloading the integrated, clustered, annotated, or trajectory object is
  about 4x faster: objects are saved with fast gzip compression (files are
  ~10% larger and still open with `readRDS()`). The download button now stays
  busy until the file is ready instead of resetting after two minutes, and its
  label says that large datasets can take several minutes.
* Doublet detection no longer fails with `could not find function "assay<-"`
  after Monocle 3 has been used in the same R session.
* Feature plots split by sample no longer fail when more genes than samples
  are selected, and their height now fits one row of sample panels per gene.
* The Docker path setup snippet now shows one command per line.
* The home page shows the installed Asc-Seurat version.
* Server-side bookmarks are stored in the app work directory (the installed
  package folder is read-only in Docker), and restoring a bookmark no longer
  re-runs the analysis steps whose buttons had been clicked.
* Differential expression and SingleR show a clear message when a required
  selection is missing instead of an empty error.
* The "Generating expression plots" alert no longer stays on screen after the
  plots have rendered.
* The session report records demo-data loading.
* Trajectory gene-table and object downloads show a message instead of saving
  an empty file when the analysis has not been run.

## Performance

* The Docker image includes `presto`, which speeds up Wilcoxon marker
  detection. MAST is offered as a DE test only when it is installed.
* Integration no longer runs PCA twice or builds neighbors and a UMAP that
  clustering recomputes; the "not clustered" integrated download therefore no
  longer contains a UMAP.
* tSNE is optional in single-sample clustering.
* Integration reuses the counts read for the per-sample QC preview instead of
  reading each 10X, h5, CSV, or TSV sample from disk again.
* SCTransform workflows no longer store a second dense scaled matrix for the
  RNA assay; scaled expression plots compute it for the requested genes.
* tradeSeq fits genes in parallel (up to four workers).

# ascseurat 3.0.7

## New features

* Gene visualization plots can split clusters by sample for integrated data.
* Data-path help explains how Docker bind mounts map host folders into the app.

# ascseurat 3.0.6

## Install and distribution

* ARM64 Docker builds now switch Ubuntu ports package sources to an HTTPS
  mirror before `apt-get update`, avoiding release failures when
  `ports.ubuntu.com` is unreachable over plain HTTP on GitHub runners.

# ascseurat 3.0.5

## Install and distribution

* Docker builds now retry transient Ubuntu apt mirror and CRAN/RSPM download
  failures before failing release builds.
* macOS CI now retries CRAN bootstrap installs for `remotes` and `pak`, with a
  longer timeout for large binary downloads.

# ascseurat 3.0.4

## Install and distribution

* GitHub Actions now pre-install `BPCells` with `remotes` before handing the
  rest of the dependency graph to `pak`, avoiding the GitHub sub-directory
  archive unpacking failure that broke pkgdown and R CMD check.
* macOS CI now installs and exports OpenSSL/HDF5 build paths before source
  package compilation, fixing Bioconductor packages that link against
  `libssl`.
* Docker release publishing now downloads per-architecture digest artifacts by
  name before creating the multi-arch manifest.

# ascseurat 3.0.3

## Install and distribution

* R packages used by app features are now declared as runtime dependencies.
  Local R installs therefore install `PseudotimeDE`, Monocle 3, SingleR,
  scDblFinder, and other app packages automatically, while docs now list the
  system compiler/HDF5 prerequisites that users still need to install.
* README and ReadTheDocs installation instructions now direct users to either
  Docker or the R-session install commands, with the BPCells HDF5 prerequisite
  and macOS GNU Fortran requirement called out explicitly.
* Docker usage docs now explain bind mounts for local data paths, and the app
  strips pasted path quotes before resolving local files.
* macOS package-check smoke tests install the official R GNU Fortran runtime
  before building source packages that link Fortran libraries.

# ascseurat 3.0.2

## Install and distribution

* Docker runtime now includes `libuv1t64` for Shiny/httpuv, `libwebpmux3`
  for ragg graphics devices, and `libnode109` for V8/node-linked packages.
* Package GitHub Actions now install only hard/runtime dependencies, avoiding
  pak's optional GitHub sub-directory remote path for `BPCells`; the Docker
  release workflow remains the full optional-stack runtime gate.
* Docker smoke tests now preserve container logs on failure.

# ascseurat 3.0.1

## Install and distribution

* README install instructions now pre-install `BPCells` via `remotes` before
  running `pak::pkg_install`. Works around a `pak` 0.9.5 issue with GitHub
  sub-directory remotes that blocked the install of `monocle3` and its
  downstream packages.
* Docker image runtime stage slimmed: dropped compile-time `*-dev` packages
  and now relies on the runtime shared libraries already shipped in the
  `rocker/r-ver:4.5.3` (Ubuntu 24.04) base. Local image is ~4.8 GB.
* Docker build CI split into native `linux/amd64` and `linux/arm64` jobs
  (no more QEMU emulation) and now triggers only on `v*` release tags or
  manual workflow dispatch.
* Docker Hub multi-arch manifest: `pereiralabbio/asc-seurat:3` and `:latest`
  now resolve transparently to the native variant on Apple Silicon and x86.
* GitHub releases auto-archive to Zenodo.

## Docs

* Read the Docs site rewritten to match the v3 editorial UI and to drop
  the v2 migration-narrative framing throughout.

# ascseurat 3.0.0

## Major changes

* Refactored the codebase from a flat Shiny app into a proper R package.
  The app is now installable via `remotes::install_github("Pereira-Lab-UF/asc_seurat")`
  and launched with `ascseurat::run_app()`.
* Updated to Seurat v5 throughout. Old Seurat v3/v4 RDS files are
  automatically updated on load via `UpdateSeuratObject()`.
* New modular Shiny architecture (`mod_*_ui` / `mod_*_server`).

## Removed (intentional)

* **Dynverse / dyno trajectory inference**. Removed because Dynverse required
  Docker-in-Docker, was slow, and frequently broke between releases.
  Replaced by supported Slingshot, PAGA, and Monocle 3 workflows.
* **biomaRt online annotation**. Removed because it required a working
  internet connection and the Ensembl BioMart server was frequently down.
  Replaced by offline GO enrichment via `clusterProfiler` + organism-specific
  Bioconductor `org.*.db` packages.
* **`rclipboard` copy buttons**. Not used in practice and added a JS
  dependency.

## Added

* H5AD (AnnData) input via `anndataR`.
* RPCA integration (Seurat v5 `IntegrateLayers`).
* PAGA and Monocle 3 trajectory inference methods in the Shiny workflow.
* `PseudotimeDE-fast` as the recommended default trajectory differential
  expression engine for Slingshot outputs, with `scMaSigPro` and optional
  `tradeSeq` alternatives.
* Cell-type renaming UI in the clustering tab.
* Per-gene plot-bundle (zip) download.

# ascseurat 2.2.1 and earlier

See: <https://github.com/Pereira-Lab-UF/asc_seurat/releases>
