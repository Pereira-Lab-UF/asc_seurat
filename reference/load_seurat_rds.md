# Load a Seurat RDS File with Backward Compatibility

Loads an RDS file and automatically updates old Seurat objects (v2/v4)
to the current Seurat v5 format.

## Usage

``` r
load_seurat_rds(path)
```

## Arguments

- path:

  Character. Path to the RDS file.

## Value

A Seurat object (v5 format).
