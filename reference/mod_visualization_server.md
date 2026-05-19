# Visualization Module — Server

Visualization Module — Server

## Usage

``` r
mod_visualization_server(id, seurat_obj, markers_de = NULL, de_step_number = 5)
```

## Arguments

- id:

  Module namespace ID.

- seurat_obj:

  Reactive clustered Seurat object.

- markers_de:

  Optional reactive returning the markers data frame from the DE module
  (Step 5). When supplied, users can select genes directly from DE
  results instead of uploading a file.

- de_step_number:

  Numeric DE step number to use in user messages.
