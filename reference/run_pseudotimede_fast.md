# Run PseudotimeDE-fast on Slingshot Output

Run PseudotimeDE-fast on Slingshot Output

## Usage

``` r
run_pseudotimede_fast(sce, gene_subset = NULL, cores = 2)
```

## Arguments

- sce:

  A SingleCellExperiment with trajectory_path and trajectory_pseudotime
  columns.

- gene_subset:

  Optional character vector of genes to test.

- cores:

  Number of workers to use.

## Value

A list containing the engine name and result table.
