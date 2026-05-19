# Calculate Robust QC Limits

Uses median +/- n MADs for numeric QC metrics.

## Usage

``` r
qc_mad_limits(values, nmads = 3, lower_floor = 0, upper_ceiling = Inf)
```

## Arguments

- values:

  Numeric vector.

- nmads:

  Number of MADs from the median.

- lower_floor:

  Minimum allowed lower bound.

- upper_ceiling:

  Maximum allowed upper bound.

## Value

Named numeric vector with lower and upper values.
