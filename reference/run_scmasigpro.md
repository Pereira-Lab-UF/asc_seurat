# Run scMaSigPro on Slingshot Output

Run scMaSigPro on Slingshot Output

## Usage

``` r
run_scmasigpro(
  sce,
  gene_subset = NULL,
  p_value = 0.05,
  rsq = 0.6,
  poly_degree = 2
)
```

## Arguments

- sce:

  A SingleCellExperiment with trajectory_path and trajectory_pseudotime
  columns.

- gene_subset:

  Optional character vector of genes to test.

- p_value:

  Numeric significance cutoff.

- rsq:

  Numeric R-squared cutoff.

- poly_degree:

  Integer polynomial degree for the model.

## Value

A list containing the fitted scMaSigPro object and a result table.
