## Dynamic structural equation models

[![](https://www.r-pkg.org/badges/version/dsem)](https://cran.r-project.org/package=dsem)
[![](https://cranlogs.r-pkg.org/badges/dsem)](https://cran.r-project.org/package=dsem)
[![](https://cranlogs.r-pkg.org/badges/grand-total/dsem)](https://cran.r-project.org/package=dsem)
[![Codecov test coverage](https://codecov.io/gh/James-Thorson-NOAA/dsem/branch/test_codecov/graph/badge.svg)](https://app.codecov.io/gh/James-Thorson-NOAA/dsem?branch=test_codecov)
[![Documentation](https://img.shields.io/badge/documentation-dsem-orange.svg?colorB=E91E63)](https://james-thorson-noaa.github.io/dsem/)


Package _dsem_ fits dynamic structural equation models, which includes as nested submodels:

1. structural equation models
2. vector autoregressive models
3. dynamic factor analysis
4. state-space autoregressive integrated moving average (ARIMA) models

The model has several advantages:

* It estimates direct, indirect, and total effects among system variables, including simultaneous and lagged effects and recursive (cyclic) dependencies
* It can estimate the cumulative outcome from press or pulse experiments or initial conditions that differ from the stationary distribution of system dynamics
* It estimates structural linkages as regression slopes while jointly imputing missing values and/or measurement errors
* It is rapidly fitted as a Gaussian Markov random field (GMRF) in a Generalized Linear Mixed Model (GLMM), with speed and asymptotics associated with each
* It allows granular control over the number of parameters (and restrictions on parameters) used to structure the covariance among variables and over time,

_dsem_ is specifically intended as a minimal implementation, and uses standard packages to simplify input/output formatting:

* Input: time-series defined using class _ts_, with `NA` for missing values
* Input: structural trade-offs specified using syntax defined by package _sem_
* Output: visualizing estimated trade-offs using _igraph_
* Output: access model output using standard S3-generic functions including `summary`, `predict`, `residuals`, `simulate`, and `AIC`

Please see package vignettes for more details regarding syntax and features.

### Citation

Thorson, J. T., Andrews, A. G., Essington, T., & Large, S. (2024). Dynamic structural equation models synthesize ecosystem dynamics constrained by ecological mechanisms. Methods in Ecology and Evolution 15(4): 744-755. https://doi.org/10.1111/2041-210X.14289


