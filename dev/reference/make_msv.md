# Make text for multivariate stochastic volatility model

Make the text string for a multivariate stochastic volatility model
using arrow-lag-slope notation for DSEM.

## Usage

``` r
make_msv(
  variables,
  n_factors,
  factor_names = paste0("F", seq_len(n_factors)),
  collapse_text = TRUE
)
```

## Arguments

- variables:

  Character string of variables (i.e., column names of `tsdata`).

- n_factors:

  Number of stochastic volatility factors.

- factor_names:

  Optional character-vector of factor names, which must match NA columns
  in `tsdata`.

- collapse_text:

  whether to collapse text into long character string

## Value

A text string to be passed to
[`dsem`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem.md)
