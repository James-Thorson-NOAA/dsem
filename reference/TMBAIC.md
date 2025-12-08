# Calculate marginal AIC for a fitted model

`TMBAIC` calculates AIC for a given model fit

## Usage

``` r
TMBAIC(opt, k = 2, n = Inf)
```

## Arguments

- opt:

  the output from `nlminb` or `optim`

- k:

  the penalty on additional fixed effects (default=2, for AIC)

- n:

  the sample size, for use in AICc calculation (default=Inf, for which
  AICc=AIC)

## Value

AIC, where a parsimonious model has a AIC relative to other candidate
models
