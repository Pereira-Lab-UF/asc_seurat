# Run Integration with RPCA or Harmony

Integrates multiple samples using either RPCA or Harmony.

## Usage

``` r
run_integration(
  obj,
  method = "rpca",
  dims = 1:30,
  sample_col = "orig.ident",
  normalization_method = "lognorm",
  nfeatures = 3000,
  scale_factor = 10000
)
```

## Arguments

- obj:

  A merged Seurat object with multiple samples.

- method:

  Character. "rpca" or "harmony".

- dims:

  Integer vector. Dimensions for integration.

- sample_col:

  Character. Column in metadata identifying samples.

- normalization_method:

  Character. "lognorm" or "sctransform".

- nfeatures:

  Integer. Number of variable features.

- scale_factor:

  Numeric. Scale factor for LogNormalize.

## Value

An integrated Seurat object.
