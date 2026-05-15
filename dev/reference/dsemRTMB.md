# Fit dynamic structural equation model

Fits a dynamic structural equation model

## Usage

``` r
dsemRTMB(
  sem,
  tsdata,
  family = rep("fixed", ncol(tsdata)),
  estimate_delta0 = FALSE,
  log_prior = function(p) 0,
  control = dsem_control(),
  covs = colnames(tsdata)
)
```

## Arguments

- sem:

  Specification for time-series structural equation model structure
  including lagged or simultaneous effects. See Details section in
  [`make_dsem_ram`](https://james-thorson-NOAA.github.io/dsem/dev/reference/make_dsem_ram.md)
  for more description

- tsdata:

  time-series data, as outputted using
  [`ts`](https://rdrr.io/r/stats/ts.html), with `NA` for missing values.

- family:

  A named list of families, each returning a class `family`, including
  [`fixed()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/fixed.md),
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html),
  [`Gamma()`](https://rdrr.io/r/stats/family.html),
  [`poisson()`](https://rdrr.io/r/stats/family.html),
  [`lognormal()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/lognormal.md),
  [`tweedie()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/tweedie.md),
  or
  [`gaussian_fixed_sd()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/gaussian_fixed_sd.md)
  with names that match levels of `colnames(tsdata)` to allow different
  families by variable. Family
  [`fixed()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/fixed.md)
  specifies that states are known (i.e., measurements for that variable
  have no error). Other families allow users to supply a link function
  including `identity`, `log`, `logit`, or `cloglog`. For example
  `family = list(y = binomial("logit"), x = fixed())` would specify
  logit-linked Bernoulli distribution for variable `tsdata$y` and a
  fixed (no measurement error) distribution for `tsdata$x`. For many
  variables, it is convenient to do e.g.,
  `family = Map(function(.) gaussian(), colnames(tsdata))` rather than
  writing them all manually.

- estimate_delta0:

  Boolean indicating whether to estimate deviations from equilibrium in
  initial year as fixed effects, or alternatively to assume that
  dynamics start at some stochastic draw away from the stationary
  distribution

- log_prior:

  A user-provided function that takes as input the list of parameters
  `out$obj$env$parList()` where `out` is the output from `dsemRTMB()`,
  and returns the log-prior probability. For example
  `log_prior = function(p) dnorm( p$beta_z[1], mean=0, sd=0.1, log=TRUE)`
  specifies a normal prior probability for the first path coefficient
  with mean of zero and sd of 0.1. Note that the user must load RTMB
  using [`library(RTMB)`](https://github.com/kaskr/RTMB) prior to
  running the model.

- control:

  Output from
  [`dsem_control`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem_control.md),
  used to define user settings, and see documentation for that function
  for details.

- covs:

  optional: a character vector of one or more elements, with each
  element giving a string of variable names, separated by commas.
  Variances and covariances among all variables in each such string are
  added to the model. Warning: covs="x1, x2" and covs=c("x1", "x2") are
  not equivalent: covs="x1, x2" specifies the variance of x1, the
  variance of x2, and their covariance, while covs=c("x1", "x2")
  specifies the variance of x1 and the variance of x2 but not their
  covariance. These same covariances can be added manually via argument
  `sem`, but using argument `covs` might save time for models with many
  variables.

## Value

An object (list) of class `dsem`, fitted using RTMB

## Details

`dsemRTMB` is interchangeable with
[`dsem`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem.md),
but uses RTMB instead of TMB for estimation. Both are provided for
comparison and real-world comparison. See
[`?dsem`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem.md)
for more details
