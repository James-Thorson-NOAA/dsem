# Marginal log-likelihood

Extract the (marginal) log-likelihood of a dsem model

## Usage

``` r
# S3 method for class 'dsem'
logLik(object, ...)
```

## Arguments

- object:

  Output from
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/reference/dsem.md)

- ...:

  Not used

## Value

object of class `logLik` with attributes

- val:

  log-likelihood

- df:

  number of parameters

Returns an object of class logLik. This has attributes "df" (degrees of
freedom) giving the number of (estimated) fixed effects in the model,
abd "val" (value) giving the marginal log-likelihood. This class then
allows `AIC` to work as expected.
