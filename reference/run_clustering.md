# Run Clustering Pipeline

Performs neighbor finding, clustering, and dimensionality reduction.

## Usage

``` r
run_clustering(
  obj,
  dims = 1:30,
  resolution = 0.6,
  run_tsne = FALSE,
  reduction = "pca"
)
```

## Arguments

- obj:

  A Seurat object with PCA computed.

- dims:

  Integer vector. PCA dimensions to use.

- resolution:

  Numeric. Clustering resolution.

- run_tsne:

  Logical. Whether to also run tSNE.

- reduction:

  Character. Dimensionality reduction to use for neighbors and UMAP.
  Defaults to "pca".

## Value

A clustered Seurat object with UMAP (and optionally tSNE).
