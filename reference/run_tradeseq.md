# Run tradeSeq on Slingshot Output

Run tradeSeq on Slingshot Output

## Usage

``` r
run_tradeseq(sce, assignment, gene_subset = NULL, nknots = 6)
```

## Arguments

- sce:

  A SingleCellExperiment with Slingshot results.

- assignment:

  Slingshot assignment metadata from extract_slingshot_assignments().

- gene_subset:

  Optional character vector of genes to test.

- nknots:

  Number of knots for tradeSeq GAM fitting.

## Value

A list containing the fitted tradeSeq object and a result table.
