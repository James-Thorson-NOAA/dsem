# Simulate dsem

Plot from a fitted `dsem` model

## Usage

``` r
stepwise_selection(
  model_options,
  model_shared,
  options_initial = c(),
  quiet = FALSE,
  criterion = AIC,
  ...
)
```

## Arguments

- model_options:

  character-vector containing sem elements that could be included or
  dropped depending upon their parsimony

- model_shared:

  character-vector containing sem elements that must be included
  regardless of parsimony

- options_initial:

  character-vector containing some (possible empty) subset of
  `model_options`, where stepwise selection begins with that set of
  model options included.

- quiet:

  whether to avoid displaying progress to terminal

- criterion:

  function that computes the information criterion to be minimized,
  typically using `AIC`. However, users can instead supply a function
  that computes CIC using
  [`test_dsep`](https://james-thorson-NOAA.github.io/dsem/reference/test_dsep.md)
  and desired settings, presumably including a `set.seed` if missing
  data are being imputed

- ...:

  arguments passed to
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/reference/dsem.md),
  other than `sem` e.g., `tsdata`, `family` etc.

## Value

An object (list) that includes:

- model:

  the string with the selected SEM model

- step:

  a list, where each list element shows the models fitted during one
  step in the stepwise algorithm, from first to last step. Each step
  then lists a table, where each table row is a single fitted model, the
  first column is the AIC for that model, and the subsequent columns
  show whether each variable is included (1) or not (0)

## Details

This function conducts stepwise (i.e., forwards and backwards) model
selection using marginal AIC, while forcing some model elements to be
included and selecting among others. See `link{dsem}` for further
discussion of model selection.

## Examples

``` r
# Simulate x -> y -> z
set.seed(101)
x = rnorm(100)
y = 0.5*x + rnorm(100)
z = 1*y + rnorm(100)
tsdata = ts(data.frame(x=x, y=y, z=z))

# define candidates
model_options = c(
  "y -> z, 0, y_to_z",
  "x -> z, 0, x_to_z"
)
# define paths that are required
model_shared = "
  x -> y, 0, x_to_y
"

# Do selection
step = stepwise_selection(
  model_options = model_options,
  model_shared = model_shared,
  tsdata = tsdata,
  quiet = TRUE
)
#> List of estimated fixed and random effects:
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      4  Fixed
#> 2             mu_j                      3 Random
#> Running nlminb_loop #1
#> Running newton_loop #1
#> List of estimated fixed and random effects:
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      5  Fixed
#> 2             mu_j                      3 Random
#> Running nlminb_loop #1
#> Running newton_loop #1
#> List of estimated fixed and random effects:
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      5  Fixed
#> 2             mu_j                      3 Random
#> Running nlminb_loop #1
#> Running newton_loop #1
#> List of estimated fixed and random effects:
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      5  Fixed
#> 2             mu_j                      3 Random
#> Running nlminb_loop #1
#> Running newton_loop #1
#> List of estimated fixed and random effects:
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      4  Fixed
#> 2             mu_j                      3 Random
#> Running nlminb_loop #1
#> Running newton_loop #1
#> List of estimated fixed and random effects:
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      6  Fixed
#> 2             mu_j                      3 Random
#> Running nlminb_loop #1
#> Running newton_loop #1

# Check selected model
cat(step$model)
#> y -> z, 0, y_to_z
#> 
#>   x -> y, 0, x_to_y
```
