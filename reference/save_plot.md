# Save Plot with Configurable Options

Wrapper around ggsave/device-based saving with support for multiple
formats and configurable resolution.

## Usage

``` r
save_plot(
  plot,
  filename,
  format = "png",
  width = 15,
  height = 10,
  dpi = 300,
  path = "."
)
```

## Arguments

- plot:

  The plot object to save, or a zero-argument drawing function.

- filename:

  Character. Output filename (without extension).

- format:

  Character. One of "png", "tiff", "jpeg", "pdf", "svg".

- width:

  Numeric. Width in cm.

- height:

  Numeric. Height in cm.

- dpi:

  Integer. Resolution in DPI.

- path:

  Character. Directory to save in. Defaults to current directory.

## Value

The file path of the saved plot (invisibly).
