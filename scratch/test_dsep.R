

library(dsem)
library(ggm)
library(igraph)
library(mvtnorm)
library(dynlm)

# Confirm that graph has conditional dependencies with non-null conditioning set
A <- DAG(x5~ x3+x4, x3~ x2, x4~x2)
basiSet(A)
plot(graph_from_adjacency_matrix(A))

##############
# No time simulation
##############

# Simulation loop
pvalue_rz = array(NA, dim=c(500,3) )
for( r in seq_len(dim(pvalue_rz)[1]) ){
  set.seed(r)

  # simulate normal distribution
  a = rnorm(100)
  b = 1 + 0.5 * a + rnorm(100)
  c = 1 + 0.5 * a + rnorm(100)
  d = 2 + -0.5*b + 1*c + rnorm(100)
  data = data.frame(a=a, b=b, c=c, d=d)

  # Correct SEM
  sem1 = "
    a -> b, 0, beta_ab
    a -> c, 0, beta_ac
    b -> d, 0, beta_bd
    c -> d, 0, beta_cd
  "
  fit1 = dsem(
    tsdata = ts(data),
    sem = sem1
  )
  pvalue_rz[r,1] = test_dsep( fit1,
                              max_lag = 0 )

  # Backwards causality
  sem2 = "
    b -> a, 0, beta_ba
    c -> a, 0, beta_ca
    d -> b, 0, beta_db
    d -> c, 0, beta_dc
  "
  fit2 = dsem(
    tsdata = ts(data),
    sem = sem2
  )
  pvalue_rz[r,2] = test_dsep( fit2,
                              max_lag = 0 )

  # Flat regression structure
  sem3 = "
    a -> b, 0, beta_ab
    a -> b, 0, beta_ab
    a -> c, 0, beta_ac
  "
  fit3 = dsem(
    tsdata = ts(data),
    sem = sem2
  )
  pvalue_rz[r,3] = test_dsep( fit3,
                              max_lag = 0 )
}

# Show tests
par( mfrow=c(3,1), oma=c(2,0,0,0) )
hist(pvalue_rz[,1], breaks=seq(0,1,by=0.05), main="right model")
hist(pvalue_rz[,2], breaks=seq(0,1,by=0.05), main="backwards causality")
hist(pvalue_rz[,3], breaks=seq(0,1,by=0.05), main="regression-style structure")
mtext( side=1, text="p-value (rejecting test that model is correct)", outer=TRUE )

######### stepwise selection on both directions
model_options = c(
  "a -> b, 0, beta_ab",
  "a -> c, 0, beta_ac",
  "a -> d, 0, beta_ad",
  "b -> c, 0, beta_bc",
  "b -> d, 0, beta_bd",
  "c -> d, 0, beta_cd",

  "b -> a, 0, beta_ba",
  "c -> a, 0, beta_ca",
  "d -> a, 0, beta_da",
  "c -> b, 0, beta_cb",
  "d -> b, 0, beta_db",
  "d -> c, 0, beta_dc"
)
selex = stepwise_selection(
  model_options = model_options,
  model_shared = "",
  tsdata = ts(data)
)
cat(selex$model)

fit_selex = dsem(
  sem = selex$model,
  tsdata = ts(data)
)

#################
# VAR simulation
#################

# Simultaneous effects and covariance
P0 = rbind(
  c( 0, 1 ),
  c( 0, 0 )
)
P1 = rbind(
  c( 0, 0.5),
  c( 0, 0 )
)
V0 = solve(diag(nrow(P0))-P0) %*% t(solve(diag(nrow(P0))-P0)) # * (1 - 0.8^2)

# Simulation loop
n_replicates = 500
pvalue_rz = array(NA, dim=c(n_replicates,3) )
for( r in seq_len(n_replicates) ){
  set.seed(r)

  # simulate normal distribution
  n_obs = 100
  x_tc = array( NA, dim=c(n_obs,nrow(P0)), dimnames=list(NULL,letters[seq_len(nrow(P0))]) )
  x_tc[1,] = rmvnorm( n=1, sigma = V0 )[1,]
  for( t in 2:n_obs ){
    x_tc[t,] = (x_tc[t-1,] %*% P1) + rmvnorm( n=1, sigma = V0 )[1,]
  }
  data = data.frame(x_tc)

  # Correct SEM
  sem1 = "
    a -> b, 0, beta_ab
    a -> b, 1, rho_ab
  "
  fit1 = dsem(
    tsdata = ts(data),
    sem = sem1
  )
  test1 = test_dsep( fit1,
    max_lag = 1,
    what = "all" )
  pvalue_rz[r,1] = test1$pvalue

  # Wrong SEM
  sem2 = "
    b -> a, 0, beta_ba
    a -> b, 1, rho_ab
  "
  fit2 = dsem(
    tsdata = ts(data),
    sem = sem2
  )
  test2 = test_dsep( fit2,
    max_lag = 1,
    what = "all" )
  pvalue_rz[r,2] = test2$pvalue

  # Wrong SEM
  sem3 = "
    a -> b, 0, beta_ab
    b -> b, 1, rho_bb
  "
  fit3 = dsem(
    tsdata = ts(data),
    sem = sem3
  )
  test3 = test_dsep( fit3,
    max_lag = 1,
    what = "all" )
  pvalue_rz[r,3] = test3$pvalue

}

# Show tests
par( mfrow=c(3,1), oma=c(2,0,0,0) )
hist(pvalue_rz[,1], breaks=seq(0,1,by=0.05), main="right model")
hist(pvalue_rz[,2], breaks=seq(0,1,by=0.05), main="backwards causality")
hist(pvalue_rz[,3], breaks=seq(0,1,by=0.05), main="regression-style structure")
mtext( side=1, text="p-value (rejecting test that model is correct)", outer=TRUE )


