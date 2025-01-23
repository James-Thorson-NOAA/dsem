

setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem\scratch)' )

example = readRDS( "example.RDS" )
attach(example)

# Custom version
dGMRF = function(x, mu, Q ){
  #(2*pi)^(-n_k/2) *
  #determinant(solve(Q))^(-1/2) *
  #exp(-1/2 * (x-mu) %*% Q %*% (x-mu) )
  (-length(as.vector(x))/2) * log(2*pi) +
  (-1/2) * -1 * determinant(Q, log=TRUE)$modulus +
  (-1/2 * as.vector(x-mu) %*% Q %*% as.vector(x-mu) )[1,1]
}

#
library(RTMB)
dgmrf( y, mu=mu, Q=Q, log=TRUE )
dgmrf( as.vector(y), mu=as.vector(mu), Q=Q, log=TRUE )
mvtnorm::dmvnorm( y, mu, sigma=as.matrix(solve(Q)), log=TRUE )
mvtnorm::dmvnorm( as.vector(y), as.vector(mu), sigma=as.matrix(solve(Q)), log=TRUE )
dGMRF( y, mu=mu, Q=Q )
