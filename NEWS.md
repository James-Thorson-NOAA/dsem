# dsem 1.6.0

* Export and document `make_matrices`
* Remove `ggm` from Imports (because it Imports `graph` which is not on CRAN)
  and instead define a locally copy of relevant functions.  Also adding Giovanni M.
  Marchetti (the maintainer of `ggm` as contributor in DESCRIPTION)
* Update `test_dsep(.)` to allow options for imputing missing data prior to test,
  and using default `impute_data = "by_test"` based on explorations to date
* Add example for `total_effect`

# dsem 1.5.0

* Add option to simulate data with or without missingness pattern of original data
* Change `test_dsep(.)` to have option to impute missing data, and changing that
  to be the default behavior

# dsem 1.5.0

* Fix bug when adding missing variances (PR 36)
* Fix bug in detecting high ratio of variances (issue 35)
* Added `test_dsep(.)` to calculate a p-value for the probability that a
  a data set arose from the specified model (highly experimental)
* Added `total_effect(.)` to calculate total effects
* Fixed bug where dsem crashed when `quiet = FALSE` and running without data
  e.g., as qualitative network model
* Added vignette showing how to approximate diffusive movement using DSEM with
  paths among adjacent strata

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
