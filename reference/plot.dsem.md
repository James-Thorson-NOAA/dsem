# Simulate dsem

Plot from a fitted `dsem` model

## Usage

``` r
# S3 method for class 'dsem'
plot(
  x,
  y,
  edge_label = c("name", "value", "value_and_stars"),
  digits = 2,
  style = c("igraph", "ggraph"),
  ...
)
```

## Arguments

- x:

  Output from
  [`dsem`](https://james-thorson-NOAA.github.io/dsem/reference/dsem.md)

- y:

  Not used

- edge_label:

  Whether to plot parameter names, estimated values, or estimated values
  along with stars indicating significance at 0.05, 0.01, or 0.001
  levels (based on two-sided Wald tests)

- digits:

  integer indicating the number of decimal places to be used

- style:

  Whether to make a graph using `igraph` or `ggraph`

- ...:

  arguments passed to
  [`plot.igraph`](https://r.igraph.org/reference/plot.igraph.html)

## Value

Invisibly returns the output from
[`graph_from_data_frame`](https://r.igraph.org/reference/graph_from_data_frame.html)
which was passed to
[`plot.igraph`](https://r.igraph.org/reference/plot.igraph.html) for
plotting.

## Details

This function coerces output from a graph and then plots the graph.
