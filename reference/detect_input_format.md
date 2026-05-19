# Detect Input Data Format

Determines the format of input data (10X directory, h5, h5ad, csv, rds)
and returns a standardized format identifier.

## Usage

``` r
detect_input_format(path)
```

## Arguments

- path:

  Character. Path to the input file or directory.

## Value

Character. One of: "10x_dir", "10x_h5", "h5ad", "csv", "tsv", "rds",
"unknown".
