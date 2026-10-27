# Run PseudotimeDE-fast on Slingshot Output

Each Slingshot lineage is tested separately on its own cells, so the
result has one row per gene and lineage.

## Usage

``` r
run_pseudotimede_fast(sce, gene_subset = NULL, cores = 2)
```

## Arguments

- sce:

  A SingleCellExperiment with Slingshot results.

- gene_subset:

  Optional character vector of genes to test.

- cores:

  Number of workers to use.

## Value

A list containing the engine name and result table.
