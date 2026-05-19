# ascseurat 3.0.2

## Install and distribution

* Docker runtime now includes `libuv1t64` for Shiny/httpuv, `libwebpmux3`
  for ragg graphics devices, and `libnode109` for V8/node-linked packages.
* Package GitHub Actions now install only hard/runtime dependencies, avoiding
  pak's optional GitHub sub-directory remote path for `BPCells`; the Docker
  release workflow remains the full optional-stack runtime gate.
* R packages used by app features are now declared as runtime dependencies.
  Local R installs therefore install `PseudotimeDE`, Monocle 3, SingleR,
  scDblFinder, and other app packages automatically, while docs now list the
  system compiler/HDF5 prerequisites that users still need to install.
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
