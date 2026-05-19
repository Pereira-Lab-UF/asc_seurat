# Set Up Python Environment for PAGA

Installs and verifies the Python packages used by the supported PAGA
trajectory workflow. Docker images already include this stack.

## Usage

``` r
setup_paga(method = "auto", envname = "ascseurat-paga")
```

## Arguments

- method:

  Character. Installation method for Python. Defaults to "auto" which
  uses conda if available, otherwise virtualenv.

- envname:

  Character. Name of the Python environment. Defaults to
  "ascseurat-paga".

## Examples

``` r
if (FALSE) { # \dontrun{
ascseurat::setup_paga()
} # }
```
