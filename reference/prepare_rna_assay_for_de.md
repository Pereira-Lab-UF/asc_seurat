# Prepare RNA Assay for Differential Expression and Visualization

SCTransform is used through PCA/clustering, while RNA-normalized values
are kept available for DE and expression plots.

## Usage

``` r
prepare_rna_assay_for_de(obj, nfeatures = 3000, scale_factor = 10000)
```

## Arguments

- obj:

  A Seurat object.

- nfeatures:

  Integer. Number of variable features.

- scale_factor:

  Numeric. Scale factor for LogNormalize.

## Value

A Seurat object with RNA data and scale.data prepared when possible.
