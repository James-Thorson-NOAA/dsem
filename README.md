# dsem
Dynamic structural equation models.

[![](https://cranlogs.r-pkg.org/badges/dsem)](https://cran.r-project.org/package=dsem)

[![DOI](https://zenodo.org/badge/656795688.svg)](https://zenodo.org/doi/10.5281/zenodo.10304770)

## Dynamic structural equation models
Package _dsem_ fits dynamic structural equation models, which includes as nested submodels:
* structural equation models
* vector autoregressive models
* dynamic factor analysis
* state-space autoregressive integrated moving average (ARIMA) models

The model has several advantages:
* It estimates direct, indirect, and total effects among system variables
* It can estimate the cumulative outcome from press or pulse experiments or initial conditions
* It jointly estimates structural linkages and imputes missing values
* It is rapidly fitted as a Gaussian Markov random field (GMRF) in a Generalized Linear Mixed Model (GLMM), with speed and asymptotics associated with each.

_phylosem_ is specifically intended as a minimal implementation, and uses standard packages for input/output formatting:
* Input: time-series defined using class _ts_, with `NA` for missing values
* Input: structural trade-offs specified using syntax defined by package _sem_
* Output: visualizing estimated trade-offs using _igraph_
* Output: access model output using standard S3-generic functions including `summary`, `predict`, `residuals`, and `simulate`

Please see package vignettes for more details regarding syntax and features.
