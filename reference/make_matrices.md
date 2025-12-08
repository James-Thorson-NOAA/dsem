# Make path matrices

Constructs path matrices for dynamic structural equation model (DSEM)
using a vector of parameters and specification of the DSEM

## Usage

``` r
make_matrices(beta_p, model, times, variables)
```

## Arguments

- beta_p:

  vector parameters.

- model:

  matrix or data.frame with the following columns, and one row per
  one-headed or two-headed arrow in the dynamic structural model:

  direction

  :   whether a path coefficient is one-headed (1) or two-headed (2)

  lag

  :   whether the lag associated with a given coefficient

  start

  :   starting value, used when `parameter=0`

  parameter

  :   The parameter number from `beta_p` associated with a given path

  first

  :   The variable at the tail of a given path

  second

  :   The variable at the head of a given path

- times:

  integer-vector of times to use when defining matrices

- variables:

  character-vector listing variables

## Value

A named list of matrices including:

- P_kk:

  The matrix of interactions, i.e., one-headed arrows

- G_kk:

  The matrix of exogenous covariance, i.e., two-headed arrows

## Details

When `length(times)` is \\T\\ and `length(variables)` is \\J\\,
`make_matrices` returns matrices of dimension \\TJ \times TJ\\
representing paths among \\vec(\mathbf{X})\\ where matrix \\\mathbf{X}\\
has dimension \\T \times J\\ and \\vec\\ stacks columns into a single
long vector
