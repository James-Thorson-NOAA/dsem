# Make a RAM (Reticular Action Model)

`read_model` converts SEM arrow notation to `model` describing SEM
parameters

## Usage

``` r
read_model(sem, times, variables, covs = NULL, quiet = FALSE)
```

## Arguments

- sem:

  Specification for time-series structural equation model structure
  including lagged or simultaneous effects. See Details section in
  [`make_dsem_ram`](https://james-thorson-NOAA.github.io/dsem/dev/reference/make_dsem_ram.md)
  for more description

- times:

  A character vector listing the set of times in order

- variables:

  A character vector listing the set of variables

- covs:

  A character vector listing variables for which to estimate a standard
  deviation

- quiet:

  Boolean indicating whether to print messages to terminal

## Details

See
[`make_dsem_ram`](https://james-thorson-NOAA.github.io/dsem/dev/reference/make_dsem_ram.md)
for details
