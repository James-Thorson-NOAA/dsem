# Convert output from package dsem to phylopath

Convert dsem to phylopath output

## Usage

``` r
as_fitted_DAG(
  fit,
  lag = 0,
  what = c("Estimate", "Std_Error", "p_value"),
  direction = 1
)
```

## Arguments

- fit:

  Output from
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem.md)

- lag:

  which lag to output

- what:

  whether to output estimates `what="Estimate"`, standard errors
  `what="Std_Error"` or p-values `what="Std_Error"`

- direction:

  whether to include one-sided arrows `direction=1`, or both one- and
  two-sided arrows `direction=c(1,2)`

## Value

Convert output to format supplied by
[`est_DAG`](https://Ax3man.github.io/phylopath/reference/est_DAG.html)
