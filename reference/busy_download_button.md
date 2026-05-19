# Download Button with Busy Feedback

Prevents repeated clicks while a long-running download is being
prepared.

## Usage

``` r
busy_download_button(
  output_id,
  label,
  class = "btn-download-soft w-100",
  busy_label = "Preparing download...",
  timeout_ms = 120000
)
```

## Arguments

- output_id:

  Download output ID.

- label:

  Button label.

- class:

  CSS class.

- busy_label:

  Label shown while the file is being prepared.

- timeout_ms:

  Fallback time before the button is re-enabled.

## Value

A download button UI element.
