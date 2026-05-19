# Create a Download Handler for Plots

Factory function that creates a Shiny downloadHandler for plots with
user-configurable format, size, and resolution.

## Usage

``` r
plot_download_handler(plot_reactive, filename_base, input, prefix)
```

## Arguments

- plot_reactive:

  A reactive expression that returns the plot.

- filename_base:

  Character. Base filename (no extension).

- input:

  The Shiny input object.

- prefix:

  Character. Prefix for input IDs (e.g., "p1" matches "p1_height",
  "p1_width", "p1_res", "p1_format").

## Value

A
[`downloadHandler`](https://rdrr.io/pkg/shiny/man/downloadHandler.html).
