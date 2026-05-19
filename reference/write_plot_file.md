# Write a Plot to Disk

Saves either a ggplot/patchwork object or a drawing function.

## Usage

``` r
write_plot_file(plot, file, format = "png", width = 15, height = 10, dpi = 300)
```

## Arguments

- plot:

  The plot object to save, or a zero-argument drawing function.

- file:

  Character. Output file path.

- format:

  Character. One of "png", "tiff", "jpeg", "pdf", "svg".

- width:

  Numeric. Width in cm.

- height:

  Numeric. Height in cm.

- dpi:

  Integer. Resolution in DPI.
