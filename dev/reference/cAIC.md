# Calculate conditional AIC

Calculates the conditional Akaike Information criterion (cAIC).

## Usage

``` r
cAIC(object, what = c("cAIC", "EDF"))
```

## Arguments

- object:

  Output from
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem.md)

- what:

  Whether to return the cAIC or the effective degrees of freedom (EDF)
  for each group of random effects.

## Value

Either the cAIC, or the effective degrees of freedom (EDF) by group of
random effects

## Details

cAIC is designed to optimize the expected out-of-sample predictive
performance for new data that share the same random effects as the
in-sample (fitted) data, e.g., spatial interpolation. In this sense, it
should be a fast approximation to optimizing the model structure based
on k-fold crossvalidation. By contrast, `AIC` calculates the marginal
Akaike Information Criterion, which is designed to optimize expected
predictive performance for new data that have new random effects, e.g.,
extrapolation, or inference about generative parameters.

cAIC also calculates as a byproduct the effective degrees of freedom,
i.e., the number of fixed effects that would have an equivalent impact
on model flexibility as a given random effect.

Both cAIC and EDF are calculated using Eq. 6 of Zheng Cadigan Thorson
2024.

Note that, for models that include profiled fixed effects, these
profiles are turned off.

## References

**Deriving the general approximation to cAIC used here**

Zheng, N., Cadigan, N., & Thorson, J. T. (2024). A note on numerical
evaluation of conditional Akaike information for nonlinear mixed-effects
models (arXiv:2411.14185). arXiv.
[doi:10.48550/arXiv.2411.14185](https://doi.org/10.48550/arXiv.2411.14185)

**The utility of EDF to diagnose hierarchical model behavior**

Thorson, J. T. (2024). Measuring complexity for hierarchical models
using effective degrees of freedom. Ecology, 105(7), e4327
[doi:10.1002/ecy.4327](https://doi.org/10.1002/ecy.4327)
