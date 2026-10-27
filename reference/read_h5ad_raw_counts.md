# Read Raw Counts from the .raw Slot of an h5ad File

anndataR does not import `adata.raw`, where Scanpy and CELLxGENE files
often keep the raw counts behind a normalized `X`.

## Usage

``` r
read_h5ad_raw_counts(path)
```

## Arguments

- path:

  Character. Path to the h5ad file.

## Value

A sparse genes x cells matrix, or NULL when the file has no `.raw`
matrix.
