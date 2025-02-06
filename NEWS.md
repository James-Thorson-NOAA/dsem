# dsem 1.4.2

* Added `total_effect(.)` to calculate total effects
* Fixed bug where dsem crashed when `quiet = FALSE` and running without data
  e.g., as qualitative network model
* Fixed bug where dsem crashed when `quiet = FALSE` and running without data
  e.g., as qualitative network model

# dsem 1.4.1

* Edits to allow DSEM to run when all path coefficients are fixed 
  (contributed by Maurice Goodman)

# dsem 1.4.0

* Adding option for lower and upper bounds
* Adding `stepwise_selection` for automated stepwise model selection
* Adding `plot` option to use `ggraph` as alternative to previous `igraph`
  option
* Adding `convert_equations` to extend `sem::specifyEquations` and simplify
  specification for large models
* Adding argument `prior_negloglike` as interface to specify Bayesian priors
  and/or likelihood penalties in `dsem`
* Adding `dsemRTMB` using RTMB as interchangeable alternative to `dsem` using
  TMB

# dsem 1.3.0

* Adding option to specify constant marginal variance, as alternative to existing
  calculation which results in a constant conditional variance but a changing marginal 
  variance along the causal graph

# dsem 1.2.1

* removing `checkDepPackageVersion(dep_pkg="Matrix", this_pkg="TMB")` from `.onLoad()`
  as requested by K. Kristensen

# dsem 1.2.0

* Adding option to specify covariance via argument `covs`
* Adding options to specify gmrf_parameterization="projection"
* Adding vigette outlining how to fit dynamic factor analysis
* Fix bug arising when `tsdata` had two or more columns sharing a single variable name
* Adding `make_dfa` helper function
* Updating bering_sea dataset to include extra year of cold-pool, and changing vignette
  to match the published specification and results
* Updating the lag indexing for the Klein-I model in the vignette to use positive values
  for lags, and updating saved MCMC results to match that corrected specification

# dsem 1.1.0

* Adding option to specify covariance in Gamma

# dsem 1.0.2

* Eliminate `eval` usage
* Add codecov Action and badge
* Change default behavior so that all variables in `tsdata` have a standard
  deviation by default

# dsem 1.0.1

* Fix bug in `simulate.dsem` to keep up with changing interface in `dsem`
* Update CITATION to indicate accepted paper
* Remove fit_tmb to eliminate cryptic warning messages and simplify code

# dsem 1.0.0

* Initial public release
