# Detailed control for dsem structure

Define a list of control parameters. Note that the format of this input
is likely to change more rapidly than that of
[`dsem`](https://james-thorson-NOAA.github.io/dsem/reference/dsem.md)

## Usage

``` r
dsem_control(
  nlminb_loops = 1,
  newton_loops = 1,
  trace = 0,
  eval.max = 1000,
  iter.max = 1000,
  getsd = TRUE,
  quiet = FALSE,
  run_model = TRUE,
  build_model = TRUE,
  gmrf_parameterization = c("separable", "projection", "mvn_project"),
  constant_variance = c("conditional", "marginal", "diagonal"),
  use_REML = TRUE,
  profile = NULL,
  parameters = NULL,
  map = NULL,
  getJointPrecision = FALSE,
  extra_convergence_checks = TRUE,
  lower = -Inf,
  upper = Inf,
  project_k = NULL,
  suppress_nlminb_warnings = TRUE,
  stabilize_Q = FALSE
)
```

## Arguments

- nlminb_loops:

  Integer number of times to call
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- newton_loops:

  Integer number of Newton steps to do after running
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- trace:

  Parameter values are printed every `trace` iteration for the outer
  optimizer. Passed to `control` in
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- eval.max:

  Maximum number of evaluations of the objective function allowed.
  Passed to `control` in
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- iter.max:

  Maximum number of iterations allowed. Passed to `control` in
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html).

- getsd:

  Boolean indicating whether to call
  [`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html)

- quiet:

  Boolean indicating whether to run model printing messages to terminal
  or not;

- run_model:

  Boolean indicating whether to estimate parameters (the default), or
  instead to return the model inputs and compiled TMB object without
  running;

- build_model:

  Boolean indicating whether to return inputs to `MakeADFun`

- gmrf_parameterization:

  Parameterization to use for the Gaussian Markov random field, where
  the default `separable` constructs a precision matrix that must be
  full rank, `projection` constructs a full-rank and IID precision for
  variables over time and then projects this using the inverse-cholesky
  of the precision (which allows a rank-deficient GMRF, but does not
  allow `family = "fixed"`), and `mvn_project` uses the dense variance
  for the full-rank component of the GMRF and then projects values to
  the rank-deficient component (which allows `family = "fixed"` but is
  much slower).

- constant_variance:

  Whether to specify a constant conditional variance \\ \mathbf{\Gamma
  \Gamma}^t\\ using the default `constant_variance="conditional"`, which
  results in a changing marginal variance along the specified causal
  graph when lagged paths are present. Alternatively, the user can
  specify a constant marginal variance using
  `constant_variance="diagonal"` or `constant_variance="marginal"`, such
  that \\ \mathbf{\Gamma}\\ and \\\mathbf{I-P}\\ are rescaled to achieve
  this constraint. All options are equivalent when the model includes no
  lags (only simultaneous effects) and no covariances (no two-headed
  arrows). `"diagonal"` and `"marginal"` are equivalent when the model
  includes no covariances. Given some exogenous covariance,
  `constant_variance = "diagonal"` preserves the conditional correlation
  and has changing conditional variance, while
  `constant_variance = "marginal"` has changing conditional correlation
  along the causal graph.

- use_REML:

  Boolean indicating whether to treat non-variance fixed effects as
  random, either to motigate bias in estimated variance parameters or
  improve efficiency for parameter estimation given correlated fixed and
  random effects

- profile:

  Parameters to profile out of the likelihood (this subset will be
  appended to `random` with Laplace approximation disabled).

- parameters:

  list of fixed and random effects, e.g., as constructed by `dsem` and
  then modified by hand (only helpful for advanced users to change
  starting values or restart at intended values)

- map:

  list of fixed and mirrored parameters, constructed by `dsem` by
  default but available to override this default and then pass to
  [`MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html)

- getJointPrecision:

  whether to get the joint precision matrix. Passed to
  [`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html).

- extra_convergence_checks:

  Boolean indicating whether to run extra checks on model convergence.

- lower:

  vectors of lower bounds, replicated to be as long as start and passed
  to [`nlminb`](https://rdrr.io/r/stats/nlminb.html). If unspecified,
  all parameters are assumed to be unconstrained.

- upper:

  vectors of upper bounds, replicated to be as long as start and passed
  to [`nlminb`](https://rdrr.io/r/stats/nlminb.html). If unspecified,
  all parameters are assumed to be unconstrained.

- project_k:

  logical-vector only used when `gmrf_parameterization=="mvn_project"`,
  with length `dim(tsdata)` indicating which state-values should be
  projected determinstically while ignoring measurements. This is useful
  e.g., for determinstic composite variables. If `project_k=NULL` (the
  default), then it projects any variable with zero exogenous variance.
  However, this then ignores any measurements for those variables.
  Manual specification can be used to confidition upon measurements for
  determinstic variables, but identifiability conditions are then hard
  to determine automatically.

- suppress_nlminb_warnings:

  whether to suppress uniformative warnings from `nlminb` arising when a
  function evaluation is NA, which are then replaced with Inf and
  avoided during estimation

- stabilize_Q:

  add `stability_eps = 1e-10` to stabilize precision (experimental)

## Value

An S3 object of class "dsem_control" that specifies detailed model
settings, allowing user specification while also specifying default
values
