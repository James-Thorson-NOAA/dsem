# Classical Runge-Kutta for system of equations

Classical Runge-Kutta of order 4.

## Usage

``` r
rk4sys(f, a, b, y0, n, Pars, ...)
```

## Arguments

- f:

  function in the differential equation \\y' = f(x, y)\\; defined as a
  function \\R \times R^m \rightarrow R^m\\, where \\m\\ is the number
  of equations.

- a:

  starting time for the interval to integrate

- b:

  ending time for the interval to integrate.

- y0:

  starting values at time `a`

- n:

  the number of steps from a to b.

- Pars:

  named list of parameters passed to f

- ...:

  additional inputs to function `f`

## Value

List with components x for grid points between a and b and y an n-by-m
matrix with solutions for variables in columns, i.e. each row contains
one time stamp.

## Details

Classical Runge-Kutta of order 4 for (systems of) ordinary differential
equations with fixed step size. Copied from pracma under GPL-3, with
small modifications to work with RTMB
