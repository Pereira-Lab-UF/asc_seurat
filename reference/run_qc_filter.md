# Run QC on a Seurat Object

Adds percentage mitochondrial content and filters cells.

## Usage

``` r
run_qc_filter(
  obj,
  mito_pattern = "^MT-",
  min_features = 200,
  max_features = Inf,
  max_mito_pct = 5
)
```

## Arguments

- obj:

  A Seurat object.

- mito_pattern:

  Character. Regex to identify mitochondrial genes.

- min_features:

  Integer. Min genes per cell.

- max_features:

  Integer. Max genes per cell.

- max_mito_pct:

  Numeric. Max mitochondrial percentage.

## Value

A filtered Seurat object.
