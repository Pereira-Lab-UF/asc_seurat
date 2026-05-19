# Read a Delimited Text Table

Reads CSV, TSV, or TXT input based on file extension.

## Usage

``` r
read_uploaded_table(
  path,
  header = TRUE,
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
```

## Arguments

- path:

  Character. Path to the table.

- header:

  Logical. Whether the file has a header row.

- row.names:

  Optional row names column.

- check.names:

  Logical. Passed through to the reader.

- stringsAsFactors:

  Logical. Passed through to the reader.

## Value

A data frame.
