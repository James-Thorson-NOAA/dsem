# Make text for dynamic factor analysis

Make the text string for a dynamic factor analysis expressed using
arrow-and-lag notation for DSEM.

## Usage

``` r
make_dfa(variables, n_factors, factor_names = paste0("F", seq_len(n_factors)))
```

## Arguments

- variables:

  Character string of variables (i.e., column names of `tsdata`).

- n_factors:

  Number of factors.

- factor_names:

  Optional character-vector of factor names, which must match NA columns
  in `tsdata`.

## Value

A text string to be passed to
[`dsem`](https://james-thorson-NOAA.github.io/dsem/reference/dsem.md)
