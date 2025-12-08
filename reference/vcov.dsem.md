# Extract Variance-Covariance Matrix

extract the covariance of fixed effects, or both fixed and random
effects.

## Usage

``` r
# S3 method for class 'dsem'
vcov(object, which = c("fixed", "random", "both"), ...)
```

## Arguments

- object:

  output from `dsem`

- which:

  whether to extract the covariance among fixed effects, random effects,
  or both

- ...:

  ignored, for method compatibility

## Value

A square matrix containing the estimated covariances among the parameter
estimates in the model. The dimensions dependend upon the argument
`which`, to determine whether fixed, random effects, or both are
outputted.
