# Run Normalization Pipeline

Normalizes data using either LogNormalize or SCTransform v2.

## Usage

``` r
run_normalization(
  obj,
  method = "lognorm",
  nfeatures = 3000,
  scale_factor = 10000,
  vars_to_regress = NULL,
  score_cell_cycle = FALSE,
  cc_organism = "human",
  regress_cc = FALSE
)
```

## Arguments

- obj:

  A Seurat object.

- method:

  Character. "lognorm" or "sctransform".

- nfeatures:

  Integer. Number of variable features.

- scale_factor:

  Numeric. Scale factor for LogNormalize.

- vars_to_regress:

  Optional character vector of metadata columns to regress.

- score_cell_cycle:

  Logical. Whether to score cell-cycle phase.

- cc_organism:

  Character. "human" or "mouse" for cell-cycle genes.

- regress_cc:

  Logical. Whether to regress out cell-cycle scores.

## Value

A normalized Seurat object with PCA computed.
