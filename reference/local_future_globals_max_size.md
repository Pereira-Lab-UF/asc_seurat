# Temporarily Increase future.globals.maxSize

Computes a memory-aware limit and applies it for the current scope.

## Usage

``` r
local_future_globals_max_size(
  fraction = 0.95,
  reserve_bytes = 512 * 1024^2,
  minimum_bytes = 2 * 1024^3,
  fallback_bytes = 8 * 1024^3
)
```

## Arguments

- fraction:

  Numeric fraction of detected available memory to permit.

- reserve_bytes:

  Numeric bytes to leave unallocated.

- minimum_bytes:

  Numeric minimum limit.

- fallback_bytes:

  Numeric fallback when memory cannot be detected.

## Value

A list containing the previous option value and the applied limit.
