# Load Single-Cell Data from Multiple Formats

Loads scRNA-seq data from various file formats and returns a Seurat
object.

## Usage

``` r
load_input_data(path, project = "AscSeurat", min.cells = 3, min.features = 200)
```

## Arguments

- path:

  Character. Path to the input file or directory.

- project:

  Character. Project name for the Seurat object.

- min.cells:

  Integer. Minimum number of cells expressing a gene (for filtering).

- min.features:

  Integer. Minimum number of features (genes) per cell.

## Value

A Seurat object.
