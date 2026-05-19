# Detect Available Memory in Bytes

Attempts to detect currently available memory, respecting container
limits when possible.

## Usage

``` r
detect_available_memory_bytes()
```

## Value

Numeric scalar with available bytes, or NA_real\_ if unavailable.
