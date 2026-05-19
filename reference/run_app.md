# Launch the Asc-Seurat Shiny Application

Starts the Asc-Seurat interactive web application for single-cell
RNA-seq analysis. Creates necessary working directories if they don't
exist.

## Usage

``` r
run_app(
  host = "127.0.0.1",
  port = 3838,
  launch.browser = TRUE,
  workdir = getwd(),
  demo = "pbmc",
  max_upload_size = 5120,
  ...
)
```

## Arguments

- host:

  Character. The host address to serve the app on. Defaults to
  "127.0.0.1" (localhost). Use "0.0.0.0" for Docker.

- port:

  Integer. The port to serve the app on. Defaults to 3838.

- launch.browser:

  Logical. Whether to open the app in a browser. Defaults to TRUE.

- workdir:

  Character. Directory used for user data, generated files, and Shiny
  runtime state. Defaults to the current working directory.

- demo:

  Character. Bundled demo dataset identifier. Defaults to "pbmc".

- max_upload_size:

  Numeric. Maximum browser upload size in megabytes. Defaults to 5120
  MB. Use file paths inside the app for very large RDS files.

- ...:

  Additional arguments passed to
  [`runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Examples

``` r
if (FALSE) { # \dontrun{
ascseurat::run_app()
} # }
```
