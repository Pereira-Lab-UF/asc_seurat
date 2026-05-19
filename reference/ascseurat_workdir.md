# Get the Asc-Seurat Working Directory

Shiny may evaluate the app from `inst/app`; this helper keeps user paths
anchored to the directory supplied to `run_app(workdir = ...)`.

## Usage

``` r
ascseurat_workdir()
```

## Value

Normalized working directory path.
