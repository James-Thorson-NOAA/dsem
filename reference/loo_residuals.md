# Calculate leave-one-out residuals

Calculates quantile residuals using the predictive distribution from a
jacknife (i.e., leave-one-out predictive distribution)

## Usage

``` r
loo_residuals(
  object,
  nsim = 100,
  what = c("quantiles", "samples", "loo"),
  track_progress = TRUE,
  ...
)
```

## Arguments

- object:

  Output from
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/reference/dsem.md)

- nsim:

  Number of simulations to use if `family!="fixed"` for some variable,
  such that simulation residuals are required.

- what:

  whether to return quantile residuals, or samples from the
  leave-one-out predictive distribution of data, or a table of
  leave-one-out predictions and standard errors for the latent state

- track_progress:

  whether to track runtimes on terminal

- ...:

  Not used

## Value

A matrix of residuals, with same order and dimensions as argument
`tsdata` that was passed to `dsem`.

## Details

Conditional quantile residuals cannot be calculated when using
`family = "fixed"`, because state-variables are fixed at available
measurements and hence the conditional distribution is a Dirac delta
function. One alternative is to use leave-one-out residuals, where we
calculate the predictive distribution for each state value when dropping
the associated observation, and then either use that as the predictive
distribution, or sample from that predictive distribution and then
calculate a standard quantile distribution for a given non-fixed family.
This appraoch is followed here. It is currently only implemented when
all variables follow `family = "fixed"`, but could be generalized to a
mix of families upon request.
