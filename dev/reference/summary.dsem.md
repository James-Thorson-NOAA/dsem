# summarize dsem

summarize parameters from a fitted dynamic structural equation model

## Usage

``` r
# S3 method for class 'dsem'
summary(object, ...)
```

## Arguments

- object:

  Output from
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem.md)

- ...:

  Not used

## Value

Returns a data.frame summarizing estimated path coefficients, containing
columns:

- path:

  The parsed path coefficient

- lag:

  The lag, where e.g. 1 means the predictor in time t effects the
  response in time t+1

- name:

  Parameter name

- start:

  Start value if supplied, and NA otherwise

- parameter:

  Parameter number

- first:

  Variable in path treated as predictor

- second:

  Variable in path treated as response

- direction:

  Whether the path is one-headed or two-headed

- Estimate:

  Maximum likelihood estimate

- Std_Error:

  Estimated standard error from the Hessian matrix

- z_value:

  Estimate divided by Std_Error

- p_value:

  P-value associated with z_value using a two-sided Wald test

## Details

A DSEM is specified using "arrow and lag" notation, which specifies the
set of path coefficients and exogenous variance parameters to be
estimated. Function `dsem` then estimates the maximum likelihood value
for those coefficients and parameters by maximizing the log-marginal
likelihood. Standard errors for parameters are calculated from the
matrix of second derivatives of this log-marginal likelihood (the
"Hessian matrix").

However, many users will want to associate individual parameters and
standard errors with the path coefficients that were specified using the
"arrow and lag" notation. This task is complicated in models where some
path coefficients or variance parameters are specified to share a single
value a priori, or were assigned a name of NA and hence assumed to have
a fixed value a priori (such that these coefficients or parameters have
an assigned value but no standard error). The `summary` function
therefore compiles the MLE for coefficients (including duplicating
values for any path coefficients that assigned the same value) and
standard error estimates, and outputs those in a table that associates
them with the user-supplied path and parameter names. It also outputs
the z-score and a p-value arising from a two-sided Wald test (i.e.
comparing the estimate divided by standard error against a standard
normal distribution).
