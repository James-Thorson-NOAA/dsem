# Simulate dsem

Simulate from a fitted `dsem` model

## Usage

``` r
# S3 method for class 'dsem'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  variance = c("none", "random", "both"),
  resimulate_gmrf = FALSE,
  fill_missing = FALSE,
  ...
)
```

## Arguments

- object:

  Output from
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/reference/dsem.md)

- nsim:

  number of simulated data sets

- seed:

  random seed

- variance:

  whether to ignore uncertainty in fixed and random effects, include
  estimation uncertainty in random effects, or include estimation
  uncertainty in both fixed and random effects

- resimulate_gmrf:

  whether to resimulate the GMRF based on estimated or simulated random
  effects (determined by argument `variance`)

- fill_missing:

  whether to fill in simulate all data (including values that are
  missing in the original data set)

- ...:

  Not used

## Value

Simulated data, either from `obj$simulate` where `obj` is the compiled
TMB object, first simulating a new GMRF and then calling `obj$simulate`.

## Details

This function conducts a parametric bootstrap, i.e., simulates new data
conditional upon estimated values for fixed and random effects. The user
can optionally simulate new random effects conditional upon their
estimated covariance, or simulate new fixed and random effects
conditional upon their imprecision.

Note that `simulate` will have no effect on states `x_tj` for which
there is a measurement and when those measurements are fitted using
`family="fixed"`, unless `resimulate_gmrf=TRUE`. In this latter case,
the GMRF is resimulated given estimated path coefficients
