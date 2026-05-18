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
