# Plot Download Options Card

Creates a bslib card containing all download options for a plot.

## Usage

``` r
card_plot_download(
  prefix,
  download_id,
  height_default = 10,
  width_default = 15,
  ns = identity,
  plot_choices = NULL,
  plot_selected = NULL,
  plot_label = "Plot to download"
)
```

## Arguments

- prefix:

  Character. Prefix for all input IDs.

- download_id:

  Character. ID for the download button.

- height_default:

  Numeric. Default height in cm.

- width_default:

  Numeric. Default width in cm.

- ns:

  Namespace function for module-scoped inputs.

- plot_choices:

  Optional named character vector/list of plot choices.

- plot_selected:

  Optional default selected plot choice.

- plot_label:

  Character. Label for the optional plot-choice selector.

## Value

A bslib card UI element.
