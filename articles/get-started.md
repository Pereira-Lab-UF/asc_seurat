# Get started

Asc-Seurat v3 can be launched either from Docker or from the GitHub R
package. For the R package install, follow the upstream [BPCells R
installation
instructions](https://github.com/bnprks/BPCells#r-installation) first so
the HDF5 system dependency is available.

``` r
# install.packages(c("pak", "remotes"))
# remotes::install_github("bnprks/BPCells/r", upgrade = "never")
# pak::pkg_install("Pereira-Lab-UF/asc_seurat")
# ascseurat::run_app()
```

For a quick first run, use the built-in 300-cell PBMC demo dataset.
