

set.seed(101)
library(dsem)
source( R'(C:\Users\James.Thorson\Desktop\Work files\Collaborations\2025 -- DSEM guidance paper\Check residuals\phylopath_utilities.R)')

# Confirm that graph has conditional dependencies with non-null conditioning set
library(ggm)
library(igraph)
A <- DAG(x5~ x3+x4, x3~ x2, x4~x2)
basiSet(A)
plot(graph_from_adjacency_matrix(A))

# simulate normal distribution
a = rnorm(100)
b = 1 + 0.5 * a + rnorm(100)
c = 1 + 0.5 * a + rnorm(100)
d = 2 + -0.5*b + 1*c + rnorm(100)
data = data.frame(a=a, b=b, c=c, d=d)

#
sem = "
  a -> b, 0, beta_ab
  a -> c, 0, beta_ac
  b -> d, 0, beta_bd
  c -> d, 0, beta_cd
"
fit = dsem(
  tsdata = ts(data),
  sem = sem
)
test_dsep( fit,
             max_lag = 0 )

#
model_options = c(
  "a -> b, 0, beta_ab",
  "a -> c, 0, beta_ac",
  "a -> d, 0, beta_ad",
  "b -> c, 0, beta_bc",
  "b -> d, 0, beta_bd",
  "c -> d, 0, beta_cd",

  "b -> a, 0, beta_ab",
  "a -> c, 0, beta_ac",
  "a -> d, 0, beta_ad",
  "b -> c, 0, beta_bc",
  "b -> d, 0, beta_bd",
  "c -> d, 0, beta_cd"
)
selex = stepwise_selection(
  model_options = model_options,
  model_shared = "",
  tsdata = ts(data)
)

sem_reversed = "
  d <- b, 0, beta_db
  d <- c, 0, beta_dc
  b <- a, 0, beta_ba
  c <- a, 0, beta_ca
"
fit_reversed = dsem(
  tsdata = ts(data),
  sem = sem_reversed
)
test_dsep( fit_reversed,
             max_lag = 0 )

