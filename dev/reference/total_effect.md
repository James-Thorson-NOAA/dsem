# Calculate total effects

Calculate a data frame of total effects, resulting from a pulse
experiment (i.e., an exogenous and temporary change in a single variable
in time `t=0`) or a press experiment (i.e., an exogenous and permanent
change in a single variable starting in time `t=0` and continuing for
`n_lags` times), representing the estimated effect of a change in any
variable on every other variable and any time-lag from 0 (simultaneous
effects) to a user-specified maximum lag.

## Usage

``` r
total_effect(object, n_lags = 4, type = c("pulse", "press"))
```

## Arguments

- object:

  Output from
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem.md)

- n_lags:

  Number of lags over which to calculate total effects

- type:

  Whether a pulse or press experiment is intended. A pulse experiment
  answers the question:
  `What happens if a variable is changed for only a single time-interval?" A press experiment answers the question: `What
  happens if a variable is permanently changed starting in a given
  time-interval?

## Value

A data frame listing the time-lag (lag), variable that is undergoing
some exogenous change (from), and the variable being impacted (to),
along with the total effect (total_effect) including direct and indirect
pathways, and the partial "direct" effect (direct_effect)

## Details

Total effects are taken from the Leontief matrix
\\\mathbf{(I-P)^{-1}}\\, where \\\mathbf{P}\\ is the path matrix across
variables and times. \\\mathbf{(I-P)}^{-1} \mathbf{\delta} \\ calculates
the effect of a perturbation represented by vector \\\mathbf{\delta}\\
with length \\n\_{\mathrm{lags}} \times n\_{\mathrm{J}}\\ where
\\n\_{\mathrm{J}}\\ is the number of variables. \\\mathbf{(I-P)}^{-1}
\mathbf{\delta} \\ calculates the total effect of a given variable
(from) upon any other variable (to) either in the same time (\\t=0\\),
or subsequent times (\\t \geq 1\\), where \\\mathbf{\delta} =
\mathbf{i}\_{\mathrm{T}} \otimes \mathbf{i}\_{\mathrm{J}}\\, where
\\\mathbf{i}\_{\mathrm{J}}\\ is one for the `from` variable and zero
otherwise. For a pulse experiment, \\\mathbf{i}\_{\mathrm{T}}\\ is one
at \\t=0\\ and zero for other times. For a press experiment,
\\\mathbf{i}\_{\mathrm{T}}\\ is one for all times.

We compute and list the total effect at each time from time `t=0` to
`t=n_lags-1`. For press experiments, this includes transient values as
the the total effect approaches its asymptotic value (if this exists) as
\\t\\ approaches infinity. If the analyst wants an asymptotic effect
from a press experiment, we recommend using a high lag (e.g.,
`n_lags = 100`) and then confirming that it has reached it's asymptote
(i.e., the total effect is almost identical for the last and
next-to-last lag), and then reporting the value for that last lag.

## Examples

``` r
### EXAMPLE 1
# Define linear model with slope of 0.5
sem = "
  # from, to, lag, name, starting_value
  x -> y, 0, slope, 0.5
"
# Build DSEM with specified value for path coefficients
mod = dsem(
  sem = sem,
  tsdata = ts(data.frame(x=rep(0,20),y=rep(0,20))),
  control = dsem_control( run_model = FALSE )
)
#> List of estimated fixed and random effects:
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      3  Fixed
#> 2             mu_j                      2 Random

# Show that total effect of X on Y from pulse experiment is 0.5 but does not propagate over time
pulse = total_effect(mod, n_lags = 2, type = "pulse")
subset( pulse, from=="x" & to=="y")
#>   lag to from total_effect direct_effect
#> 3   0  y    x          0.5           0.5
#> 4   1  y    x          0.0           0.0


### EXAMPLE 2
# Define linear model with slope of 0.5 and autocorrelated response
sem = "
  x -> y, 0, slope, 0.5
  y -> y, 1, ar_y, 0.8
"
mod = dsem(
  sem = sem,
  tsdata = ts(data.frame(x=rep(0,20),y=rep(0,20))),
  control = dsem_control( run_model = FALSE )
)
#> List of estimated fixed and random effects:
#>   Coefficient_name Number_of_coefficients   Type
#> 1           beta_z                      4  Fixed
#> 2             mu_j                      2 Random

# Show that total effect of X on Y from pulse experiment  is 0.5 with decay of 0.8 for each time
pulse = total_effect(mod, n_lags = 4, type = "pulse")
subset( pulse, from=="x" & to=="y")
#>   lag to from total_effect direct_effect
#> 5   0  y    x        0.500           0.5
#> 6   1  y    x        0.400           0.0
#> 7   2  y    x        0.320           0.0
#> 8   3  y    x        0.256           0.0

# Show that total effect of X on Y from press experiment  asymptotes at 2.5
press = total_effect(mod, n_lags = 50, type = "press")
subset( press, from=="x" & to=="y")
#>     lag to from total_effect direct_effect
#> 51    0  y    x     0.500000           0.5
#> 52    1  y    x     0.900000           0.5
#> 53    2  y    x     1.220000           0.5
#> 54    3  y    x     1.476000           0.5
#> 55    4  y    x     1.680800           0.5
#> 56    5  y    x     1.844640           0.5
#> 57    6  y    x     1.975712           0.5
#> 58    7  y    x     2.080570           0.5
#> 59    8  y    x     2.164456           0.5
#> 60    9  y    x     2.231565           0.5
#> 61   10  y    x     2.285252           0.5
#> 62   11  y    x     2.328201           0.5
#> 63   12  y    x     2.362561           0.5
#> 64   13  y    x     2.390049           0.5
#> 65   14  y    x     2.412039           0.5
#> 66   15  y    x     2.429631           0.5
#> 67   16  y    x     2.443705           0.5
#> 68   17  y    x     2.454964           0.5
#> 69   18  y    x     2.463971           0.5
#> 70   19  y    x     2.471177           0.5
#> 71   20  y    x     2.476942           0.5
#> 72   21  y    x     2.481553           0.5
#> 73   22  y    x     2.485243           0.5
#> 74   23  y    x     2.488194           0.5
#> 75   24  y    x     2.490555           0.5
#> 76   25  y    x     2.492444           0.5
#> 77   26  y    x     2.493955           0.5
#> 78   27  y    x     2.495164           0.5
#> 79   28  y    x     2.496131           0.5
#> 80   29  y    x     2.496905           0.5
#> 81   30  y    x     2.497524           0.5
#> 82   31  y    x     2.498019           0.5
#> 83   32  y    x     2.498415           0.5
#> 84   33  y    x     2.498732           0.5
#> 85   34  y    x     2.498986           0.5
#> 86   35  y    x     2.499189           0.5
#> 87   36  y    x     2.499351           0.5
#> 88   37  y    x     2.499481           0.5
#> 89   38  y    x     2.499585           0.5
#> 90   39  y    x     2.499668           0.5
#> 91   40  y    x     2.499734           0.5
#> 92   41  y    x     2.499787           0.5
#> 93   42  y    x     2.499830           0.5
#> 94   43  y    x     2.499864           0.5
#> 95   44  y    x     2.499891           0.5
#> 96   45  y    x     2.499913           0.5
#> 97   46  y    x     2.499930           0.5
#> 98   47  y    x     2.499944           0.5
#> 99   48  y    x     2.499955           0.5
#> 100  49  y    x     2.499964           0.5
```
