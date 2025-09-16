
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


#######################
# TEST OPTION
#######################

# devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)', force=TRUE )
library(dsem)

# simulate normal distribution
n_sim = 10
x = 2 + rnorm(n_sim)
y = 1 + 0.5 * x + rnorm(n_sim)
data = data.frame(x=x, y=y)

# Fit as linear model
Lm = lm( y ~ x, data=data )

# Fit as DSEM
fit = dsem( sem = "x -> y, 0, beta",
            tsdata = ts(data),
            family = c("fixed", "normal"),
            control = dsem_control( getsd = FALSE,
                                    run_model = TRUE,
                                    use_REML = FALSE ) )
# Fit as DSEM
sem = "
  x -> y, 0, beta
  x <-> x, 0, NA, 1
"
control = dsem_control( getsd = FALSE,
                                    trace = 1,
                                    run_model = TRUE,
                                    use_REML = FALSE )
control$gmrf_parameterization = "conditional_krig"
fit1 = dsem( sem = sem,
            tsdata = ts(data),
            family = c("fixed", "normal"),
            control = control )

c( AIC(fit1), AIC(Lm) )
c( fit1$internal$parhat$mu_j[2], Lm$coef[1] )
c( fit1$internal$parhat$beta_z[1], Lm$coef[2] )

#obj = fit1$obj
#opt = nlminb( obj$par, obj$fn, obj$gr, control = list(trace=1) )
#
#obj$fn( c(0,1,0) )
#rep = obj$report()
#
##
#obj = TMB::MakeADFun( data = fit1$tmb_inputs$data,
#                 parameters = fit1$tmb_inputs$par,
#                 #random = fit1$tmb_inputs$random,
#                 map = fit1$tmb_inputs$map,
#                 #profile = control$profile,
#                 DLL = "dsem",
#                 silent = TRUE )
#obj$fn(obj$par)
#obj$gr(obj$par)
#opt = nlminb( obj$par, obj$fn, obj$gr )
#H = optimHess( opt$par, obj$fn, obj$gr )
