# Gaussian with known standard deviation for measurement errors

Allows using `family = gaussian_fixed_sd(link, sd)` to specify data that
have known measurement error

## Usage

``` r
gaussian_fixed_sd(link, sd)
```

## Arguments

- link:

  Link function

- sd:

  vector of known standard deviations (repeated for length of
  time-series, so be careful)
