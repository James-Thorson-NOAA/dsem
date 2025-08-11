
library(Matrix)
library(mvtnorm)
library(tinyVAST)   # conditional_gmrf

x_s = runif( 20 )

D_ss = outer( x_s, x_s, FUN = \(a,b) abs(a-b) )
V_ss = exp( -2 * D_ss )
Q_ss = solve(V_ss)

y_s = rmvnorm( n = 1, sigma = V_ss )[1,]


#
x_g = seq(0, 1, length=101)
x_z = c(x_g, x_s)
idx = which( x_z %in% x_s )

D_zz = outer( x_z, x_z, FUN = \(a,b) abs(a-b) )
V_zz = exp( -2 * D_zz )
Q_zz = Matrix(zapsmall(solve(V_zz)))

# Use tinyVAST simulator
y_gz = conditional_gmrf(
  Q = Q_zz,
  observed_idx = idx,
  x_obs = y_s,
  n_sims = 1000
)
yrange_gz = apply( y_gz, MARGIN = 1, FUN = quantile, prob = c(0.1,0.9) )

plot( x = x_s, y = y_s, pch = 20, cex = 2 )
#matplot( x = x_g, y = y_gz, add = TRUE, type= "p", pch = 20, col = rgb(0,0,0,0.1) )
polygon(
  x = c( x_g, rev(x_g) ),
  y = c( yrange_gz[1,], rev(yrange_gz[2,]) ),
  col = rgb(0,0,1,0.2),
  border = NA
)

#setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem\src)' )
#TMB::compile( "dsem.cpp", framework = "TMBad" )

