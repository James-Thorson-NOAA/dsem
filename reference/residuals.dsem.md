# Calculate residuals

Calculate deviance or response residuals for dsem

## Usage

``` r
# S3 method for class 'dsem'
residuals(object, type = c("deviance", "response"), ...)
```

## Arguments

- object:

  Output from
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/reference/dsem.md)

- type:

  which type of residuals to compute (only option is `"deviance"` or
  `"response"` for now)

- ...:

  Not used

## Value

A matrix of residuals, with same order and dimensions as argument
`tsdata` that was passed to `dsem`.
