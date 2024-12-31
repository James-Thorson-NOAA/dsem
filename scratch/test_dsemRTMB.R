
#devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\dsem)', force=TRUE, dep=FALSE )

library(dsem)
library(RTMB)
library(Matrix)

# Define model
sem = "
  # Link, lag, param_name
  cprofits -> consumption, 0, a1
  cprofits -> consumption, 1, a2
  pwage -> consumption, 0, a3
  gwage -> consumption, 0, a3
  cprofits -> invest, 0, b1
  cprofits -> invest, 1, b2
  capital -> invest, 0, b3
  gnp -> pwage, 0, c2
  gnp -> pwage, 1, c3
  time -> pwage, 0, c1
"

# Load data
data(KleinI, package="AER")
TS = ts(data.frame(KleinI, "time"=time(KleinI) - 1931))
tsdata = TS[,c("time","gnp","pwage","cprofits",'consumption',
               "gwage","invest","capital")]

# Fit model
fit0 = dsem( sem=sem,
            tsdata = tsdata,
            estimate_delta0 = TRUE,
            family = rep("normal",ncol(tsdata)),
            control = dsem_control( quiet=TRUE,
                                    run_model = FALSE,
                                    use_REML = TRUE,
                                    gmrf_parameterization = "projection" ) )
#ParHat = fit$internal$parhat
#Rep = fit$obj$report()

#
Map = fit0$tmb_inputs$map
Map$lnsigma_j = factor( rep(NA,ncol(tsdata)) )
Params = fit0$tmb_inputs$parameters
Params$lnsigma_j[] = log(0.1)

#
prior_negloglike = \(obj) -dnorm(obj$par[1],0,0.1,log=TRUE)

# Fit model
fit = dsem( sem=sem,
            tsdata = tsdata,
            estimate_delta0 = TRUE,
            family = rep("normal",ncol(tsdata)),
            prior_negloglike = prior_negloglike,
            control = dsem_control( quiet=TRUE,
                                    run_model = TRUE,
                                    use_REML = TRUE,
                                    gmrf_parameterization = "projection",
                                    map = Map,
                                    parameters = Params ) )

# RUN dsemRTMB line-by-line
if( FALSE ){
  control = dsem_control()
  covs = colnames(tsdata)
}


###################
# dsemRTMB
###################

# Files
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "make_matrices.R") )
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "compute_nll.R") )
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "read_model.R") )
source( file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\R)', "dsemRTMB.R") )

# Define prior
log_prior = function(p) dnorm( p$beta_z[1], mean=0, sd=0.01, log=TRUE)

fitRTMB = dsemRTMB( sem = sem,
            tsdata = tsdata,
            estimate_delta0 = TRUE,
            family = rep("normal",ncol(tsdata)),
            log_prior = log_prior,
            control = dsem_control( quiet = FALSE,
                                    run_model = TRUE,
                                    use_REML = TRUE,
                                    gmrf_parameterization = "projection",
                                    map = Map,
                                    parameters = Params ) )

#
obj = fitRTMB$obj
type = tapply( summary(fit)[,'direction'], INDEX=as.numeric(summary(fit)[,'parameter']), FUN=max )
lower = rep( -Inf, length(obj$par) )
lower[names(obj$par)=="beta_z"] = ifelse( type==2, 0, -Inf )
upper = rep( Inf, length(obj$par) )

# Check bounds on priors
cbind( obj$par, lower, upper )

# Run MCMC
library(tmbstan)
mcmc = tmbstan( obj,
                lower = lower,
                upper = upper,
                init = 'last.par.best',
                laplace = FALSE )
traceplot(mcmc, pars=unique(names(obj$par)), inc_warmup=TRUE)


if( FALSE ){
  obj = fitRTMB$obj
  obj$fn(obj$par)
  obj$gr(obj$par)

  #
  #obj$fn( fit$opt$par )
  rep0 = fitRTMB$obj$report( fit$obj$env$last.par.best )
  rep1 = fit$obj$report( fit$obj$env$last.par.best )

  #
  fitRTMB$obj$gr( fit$opt$par )

  range(fit$opt$par - fitRTMB$opt$par)
}

#
if( FALSE ){
  Rep = fitRTMB$obj$report()
  dgmrf( as.vector(Rep$z_tj), mu=as.vector(Rep$xhat_tj + Rep$delta_tj), Q=Rep$Q_kk, log=TRUE )
  solve(Rep$V_kk, Rep$IminusRho_kk)

  fit$tmb_inputs$parameters$x_tj
  fitRTMB$tmb_inputs$parameters$x_tj
  fit$obj$report()$jnll_gmrf
  fitRTMB$obj$report()$jnll_gmrf

  #
  sum(dnorm( fit$tmb_inputs$parameters$x_tj, log=TRUE ))
}
