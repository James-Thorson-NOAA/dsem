# Partition variance in one variable due to another (EXPERIMENTAL)

Calculate the proportion of variance for a response variable that is
attributed to another set of predictor variables, calculated across lags
from from 0 (simultaneous effects) to a user-specified maximum lag.

## Usage

``` r
partition_variance(object, which_response, n_times = 10)
```

## Arguments

- object:

  Output from
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem.md)

- which_response:

  string matching colnames from `tsdata` identifying response variable

- n_times:

  Number of lags over which to calculate total effects

## Value

A list with two elements:

- total_variance:

  A matrix of the total variance for each variable (column) and each
  time from 1 to `n_times`

- proportion_variance_explained:

  A matrix of the proportion of variance explained for variable
  `which_response` by each model variable (column) and each time from 1
  to `n_times`

Note that in a model with lagged effects, the total_variance and
variance_explained will vary for each time (row), and the analyst might
want to either choose a time for which the value has stabilized.

## Details

This function calculates the variance for each variable and lag, and
then recalculates it when setting exogenous variance to zero for all
variables except `which_pred`. It then calculates the ratio of the
diagonal of these two. This represents the proportion of variance in the
full model that is attributable to one or more variables.

This function is under development and may still change or be removed.

## Examples

``` r
# Simulate linear model
x = rnorm(100)
y = 1 + 1 * x + rnorm(100)
data = data.frame(x=x, y=y)

# Fit as DSEM
fit = dsem( sem = "x -> y, 0, beta",
            tsdata = ts(data),
            control = dsem_control(quiet=TRUE) )

# Apply
partition_variance( fit,
                    which_response = "y",
                    n_times = 10 )
#> $total_variance
#>             x        y
#> t_1  1.104402 2.076077
#> t_2  1.104402 2.076077
#> t_3  1.104402 2.076077
#> t_4  1.104402 2.076077
#> t_5  1.104402 2.076077
#> t_6  1.104402 2.076077
#> t_7  1.104402 2.076077
#> t_8  1.104402 2.076077
#> t_9  1.104402 2.076077
#> t_10 1.104402 2.076077
#> 
#> $proportion_variance_explained
#>              x         y
#> t_1  0.5048239 0.4951761
#> t_2  0.5048239 0.4951761
#> t_3  0.5048239 0.4951761
#> t_4  0.5048239 0.4951761
#> t_5  0.5048239 0.4951761
#> t_6  0.5048239 0.4951761
#> t_7  0.5048239 0.4951761
#> t_8  0.5048239 0.4951761
#> t_9  0.5048239 0.4951761
#> t_10 0.5048239 0.4951761
#> 
```
