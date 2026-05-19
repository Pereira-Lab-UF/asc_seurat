# Build QC Filter Defaults from Metadata

Build QC Filter Defaults from Metadata

## Usage

``` r
qc_filter_defaults(
  meta,
  nmads = 3,
  defaults = list(min_genes = 200, max_genes = 5000, max_mito = 5, min_counts = 0,
    max_counts = NA_real_)
)
```

## Arguments

- meta:

  Seurat metadata data frame.

- nmads:

  Number of MADs from the median.

- defaults:

  Named fallback values.

## Value

Named list of QC default values.
