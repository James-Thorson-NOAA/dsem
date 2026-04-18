# Test d-separation

Calculate the p-value for a test of d-separation **(Experimental)**

## Usage

``` r
test_dsep(
  object,
  n_time = NULL,
  n_burnin = NULL,
  what = c("pvalue", "CIC", "all"),
  test = c("wald", "lr"),
  seed = 123456,
  order = NULL,
  impute_data = c("by_test", "single", "none")
)
```

## Arguments

- object:

  object from
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem.md)

- n_time:

  how many times to include when defining the set of conditional
  independence relationships. If missing, this value is taken from the
  maximum lag that's included in the model plus one.

- n_burnin:

  how many times to include prior to `seq_len(n_time)` when identifying
  the conditioning set that must be included when defining conditional
  independence relationships.

- what:

  whether to just get the p-value, an information criterion based on the
  conditional independence test, or a named list with these two and
  other intermediate calculations (used for diagnosing test behavior)

- test:

  whether to test each conditional-independence relationship using a
  (univariate) wald test or a (multivariate) likelihood ratio test. The
  likelihood-ratio test might be more accurate given estimation
  covariance and also faster (does not require standard errors), but
  also is not used by phylopath and therefore less supported by previous
  d-dsep testing applications.

- seed:

  random number seed used when simulating imputed data, so that results
  are reproducible.

- order:

  an optional character vector providing the order for variables to be
  tested when defining the directed acyclic graph for use in d-sep
  testing

- impute_data:

  whether to independently impute missing data for each conditional
  independence test, or to use imputed values from the original fit. The
  data are imputed separately for each conditional independence test, so
  that they are uncorrelated as expected when combining them using
  Fisher's method. Preliminary testing suggests that using imputed data
  improves test performance

## Value

A p-value representing the weight of evidence that the data arises from
the specified model, where a low p-value indicates significant evidence
for rejecting this hypothesis.

## Details

A user-specified SEM implies a set of conditional independence
relationships among variables, which can be fitted individually,
extracting the slope and associated p-value, and then combining these
p-values to define a model-wide (omnibus) p-value for the hypothesis
that a given data set arises from the specified model. This test is
modified from package:phylopath. However it is unclear exactly how to
define the set of conditional-independence assumptions in a model with
temporal autocorrelation, and the test was not developed for uses when
data are missing. At the time of writing, the function is hightly
experimental.

Note that the method is not currently designed to deal with two-headed
arrows among variables (i.e., exogenous covariance).

## References

Shipley, B. (2000). A new inferential test for path models based on
directed acyclic graphs. Structural Equation Modeling, 7(2), 206-218.
[doi:10.1207/S15328007SEM0702_4](https://doi.org/10.1207/S15328007SEM0702_4)

## Examples

``` r
# Simulate data set
set.seed(101)
a = rnorm( 100 )
b = 0.5*a + rnorm(100)
c = 1*a + rnorm(100)
d = 1*b - 0.5*c + rnorm(100)
tsdata = ts(data.frame(a=a, b=b, c=c, d=d))

# fit wrong model
wrong = dsem(
  tsdata = tsdata,
  sem = "
    a -> d, 0, a_to_d
    b -> d, 0, b_to_d
    c -> d, 0, c_to_d
  "
)
#> List of estimated fixed and random effects:
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      7  Fixed
#> 2             mu_j                      4 Random
#> Running nlminb_loop #1
#> Running newton_loop #1
test_dsep( wrong )
#> [1] 0

# fit right model
right = dsem(
  tsdata = tsdata,
  sem = "
    a -> b, 0, a_to_b
    a -> c, 0, a_to_c
    b -> d, 0, b_to_d
    c -> d, 0, c_to_d
  "
)
#> List of estimated fixed and random effects:
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      8  Fixed
#> 2             mu_j                      4 Random
#> Running nlminb_loop #1
#> Running newton_loop #1
test_dsep( right )
#> [1] 0.02636084
```
