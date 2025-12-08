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
  [`make_dsem_ram`](https://james-thorson-NOAA.github.io/dsem/reference/make_dsem_ram.md)
  for more description

- tsdata:

  time-series data, as outputted using
  [`ts`](https://rdrr.io/r/stats/ts.html), with `NA` for missing values.

- family:

  Character-vector listing the distribution used for each column of
  `tsdata`, where each element must be `fixed` (for no measurement
  error), `normal` for normal measurement error using an identity link,
  `gamma` for a gamma measurement error using a fixed CV and log-link,
  `bernoulli` for a Bernoulli measurement error using a logit-link, or
  `poisson` for a Poisson measurement error using a log-link.
  `family="fixed"` is default behavior and assumes that a given variable
  is measured exactly. Other options correspond to different
  specifications of measurement error.

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
  [`dsem_control`](https://james-thorson-NOAA.github.io/dsem/reference/dsem_control.md),
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
[`dsem`](https://james-thorson-NOAA.github.io/dsem/reference/dsem.md),
but uses RTMB instead of TMB for estimation. Both are provided for
comparison and real-world comparison. See
[`?dsem`](https://james-thorson-NOAA.github.io/dsem/reference/dsem.md)
for more details

## Examples

``` r
# Define model
sem = "
  # Link, lag, param_name
  cprofits -> consumption, 0, a1
  cprofits -> consumption, 1, a2
  pwage -> consumption, 0, a3
  gwage -> consumption, 0, a3
  cprofits -> invest, 0, b1
  cprofits -> invest, 1, b2
  capital -> invest, 0, b3
  gnp -> pwage, 0, c2
  gnp -> pwage, 1, c3
  time -> pwage, 0, c1
"

# Load data
data(KleinI, package="AER")
TS = ts(data.frame(KleinI, "time"=time(KleinI) - 1931))
tsdata = TS[,c("time","gnp","pwage","cprofits",'consumption',
               "gwage","invest","capital")]

# Fit model
fit = dsemRTMB( sem=sem,
            tsdata = tsdata,
            estimate_delta0 = TRUE,
            control = dsem_control(quiet=TRUE) )
```
