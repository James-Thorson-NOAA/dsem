# dsem 1.2.0

* Adding option to specify covariance via argument `covs`
* Adding options to specify gmrf_parameterization="projection"
* Adding vigette outlining how to fit dynamic factor analysis
* Fix bug arising when `tsdata` had two or more columns sharing a single variable name
* Adding `make_dfa` helper function

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
