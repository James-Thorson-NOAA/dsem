# Package index

## Specification and fitting

Core tools for model fitting and diagnostics.

- [`dsem()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem.md)
  : Fit dynamic structural equation model
- [`dsem_control()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/dsem_control.md)
  : Detailed control for dsem structure
- [`stepwise_selection()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/stepwise_selection.md)
  : Simulate dsem

## Specifying the DSEM

Tools to specify interactions among variables and over time.

- [`make_dsem_ram()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/make_dsem_ram.md)
  : Make a RAM (Reticular Action Model)
- [`make_dfa()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/make_dfa.md)
  : Make text for dynamic factor analysis
- [`convert_equations()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/convert_equations.md)
  : Convert equations notation

## Predicting and residuals

Core tools for model predictions and residuals.

- [`simulate(`*`<dsem>`*`)`](https://james-thorson-NOAA.github.io/dsem/dev/reference/simulate.dsem.md)
  : Simulate dsem
- [`loo_residuals()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/loo_residuals.md)
  : Calculate leave-one-out residuals
- [`residuals(`*`<dsem>`*`)`](https://james-thorson-NOAA.github.io/dsem/dev/reference/residuals.dsem.md)
  : Calculate residuals
- [`cAIC()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/cAIC.md)
  : Calculate conditional AIC
- [`logLik(`*`<dsem>`*`)`](https://james-thorson-NOAA.github.io/dsem/dev/reference/logLik.dsem.md)
  : Marginal log-likelihood

## Interpret output

Tools for interpreting output.

- [`summary(`*`<dsem>`*`)`](https://james-thorson-NOAA.github.io/dsem/dev/reference/summary.dsem.md)
  : summarize dsem
- [`test_dsep()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/test_dsep.md)
  : Test d-separation
- [`total_effect()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/total_effect.md)
  : Calculate total effects
- [`partition_variance()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/partition_variance.md)
  : Partition variance in one variable due to another (EXPERIMENTAL)
- [`vcov(`*`<dsem>`*`)`](https://james-thorson-NOAA.github.io/dsem/dev/reference/vcov.dsem.md)
  : Extract Variance-Covariance Matrix
- [`as_fitted_DAG()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/as_fitted_DAG.md)
  : Convert output from package dsem to phylopath
- [`as_sem()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/as_sem.md)
  : Convert dsem to sem output
- [`plot(`*`<dsem>`*`)`](https://james-thorson-NOAA.github.io/dsem/dev/reference/plot.dsem.md)
  : Simulate dsem
- [`predict(`*`<dsem>`*`)`](https://james-thorson-NOAA.github.io/dsem/dev/reference/predict.dsem.md)
  : predictions using dsem
- [`print(`*`<dsem>`*`)`](https://james-thorson-NOAA.github.io/dsem/dev/reference/print.dsem.md)
  : Print fitted dsem object
- [`make_matrices()`](https://james-thorson-NOAA.github.io/dsem/dev/reference/make_matrices.md)
  : Make path matrices

## Data sets

Data sets used for illustration and testing.

- [`bering_sea`](https://james-thorson-NOAA.github.io/dsem/dev/reference/bering_sea.md)
  : Bering Sea marine ecosystem
- [`isle_royale`](https://james-thorson-NOAA.github.io/dsem/dev/reference/isle_royale.md)
  : Isle Royale wolf and moose
- [`sea_otter`](https://james-thorson-NOAA.github.io/dsem/dev/reference/sea_otter.md)
  : Sea otter trophic cascade
- [`lake_washington`](https://james-thorson-NOAA.github.io/dsem/dev/reference/lake_washington.md)
  : Lake washington plankton
- [`paramesium_didinium`](https://james-thorson-NOAA.github.io/dsem/dev/reference/paramesium_didinium.md)
  : Paramesium-Didinium dynamics
- [`hare_lynx`](https://james-thorson-NOAA.github.io/dsem/dev/reference/hare_lynx.md)
  : Lynx-Hare dynamics
- [`pdo_departure_bay`](https://james-thorson-NOAA.github.io/dsem/dev/reference/pdo_departure_bay.md)
  : Pacific Decadal Oscillation and Departure Bay temperatures
