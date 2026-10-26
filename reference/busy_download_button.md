# Download Button with Busy Feedback

Prevents repeated clicks while a long-running download is being
prepared.

## Usage

``` r
busy_download_button(
  output_id,
  label,
  class = "btn-download-soft w-100",
  busy_label = "Preparing file (large datasets can take several minutes)...",
  timeout_ms = NULL
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

  Optional fallback time before the button is re-enabled. By default the
  button stays busy until the server reports that the file is ready,
  because large objects can take several minutes to prepare.

## Value

A download button UI element.
