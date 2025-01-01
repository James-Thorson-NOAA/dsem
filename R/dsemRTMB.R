

dsemRTMB <-
function( sem,
          tsdata,
          family = rep("fixed",ncol(tsdata)),
          estimate_delta0 = FALSE,
          log_prior = function(p) 0,
          control = dsem_control(),
          covs = colnames(tsdata) ){

  # General error checks
  if( isFALSE(is(control, "dsem_control")) ) stop("`control` must be made by `dsem_control()`")
  if( control$gmrf_parameterization=="projection" ){
    if( any(family=="fixed" & colSums(!is.na(tsdata))>0) ){
      stop("`family` cannot be `fixed` using `gmrf_parameterization=projection` for any variable with data")
    }
  }
  if( isFALSE(is(tsdata,"ts")) ) stop("`tsdata` must be a `ts` object")

  # General warnings
  if( isFALSE(control$quiet) ){
    tsdata_SD = apply( tsdata, MARGIN=2, FUN=sd, na.rm=TRUE )
    if( any((max(tsdata_SD)/min(tsdata_SD)) > 10, rm.na=TRUE) ){
      warning("Some variables in `tsdata` have much higher variance than others. Please consider rescaling variables to prevent issues with numerical convergence.")
    }
  }

  # Build RAM
  model = read_model( sem = sem,
            times = as.numeric(time(tsdata)),
            variables = colnames(tsdata),
            covs = covs,
            quiet = FALSE )
  ram = make_matrices(
            model = model,
            times = as.numeric(time(tsdata)),
            variables = colnames(tsdata) )

  #
  options = c(
    ifelse(control$gmrf_parameterization=="separable", 0, 1),
    switch(control$constant_variance, "conditional"=0, "marginal"=1, "diagonal"=2)
  )
  y_tj = tsdata

  #
  n_z = length(unique(model$parameter))
  n_t = nrow(y_tj)
  n_j = ncol(y_tj)
  n_k = prod(dim(y_tj))

  # Load data in environment for function "dBdt"
  data4 = local({
                  "c" <- ADoverload("c")
                  "[<-" <- ADoverload("[<-")
                  environment()
  })
  environment(log_prior) <- data4

  # Construct parameters
  if( is.null(control$parameters) ){
    Params = list(
      #beta_p2 = subset(model,direction==2)$start,
      #beta_p1 = subset(model,direction==1)$start,
      beta_z = 0.01 * rnorm(n_z),
      lnsigma_j = rep(0, n_j),
      mu_j = rep(0, n_j),
      delta0_j = rep(0, n_j),
      x_tj = ifelse( is.na(y_tj), 0, y_tj )
    )

    # Turn off initial conditions
    if( estimate_delta0==FALSE ){
      Params$delta0_j = numeric(0)
    }

    # Scale starting values with higher value for two-headed than one-headed arrows
    which_nonzero = which(model[,5]>0)
    beta_type = tapply( model[which_nonzero,8], INDEX=model[which_nonzero,5], max)
    Params$beta_z = ifelse(beta_type==1, 0.01, 1)

    # Override starting values if supplied
    which_nonzero = which(model[,5]>0)
    start_z = tapply( as.numeric(model[which_nonzero,4]), INDEX=model[which_nonzero,5], mean )
    Params$beta_z = ifelse( is.na(start_z), Params$beta_z, start_z)
  }else{
    Params = control$parameters
  }

  # Construct map
  if( is.null(control$map) ){
    Map = list()
    Map$x_tj = factor(ifelse( is.na(as.vector(y_tj)) | (family[col(y_tj)] %in% c("normal","binomial","poisson","gamma")), seq_len(n_k), NA ))
    Map$lnsigma_j = factor( ifelse(family=="fixed", NA, seq_along(Params$lnsigma_j)) )

    # Map off mean for latent variables
    Map$mu_j = factor( ifelse(colSums(!is.na(y_tj))==0, NA, seq_len(n_j)) )
  }else{
    Map = control$map
  }

  # Define random
  if(isTRUE(control$use_REML)){
    Random = c( "x_tj", "mu_j" )
  }else{
    Random = "x_tj"
  }

  # Build object
  cmb <- function(f, ...) function(p) f(p, ...) ## Helper to make closure
  #f(parlist, model, tsdata, family)
  obj = RTMB::MakeADFun(
          func = cmb( compute_nll,
                      model = model,
                      y_tj = y_tj,
                      family = family,
                      options = options,
                      log_prior = log_prior ),
          parameters = Params,
          random = Random,
          map = Map,
          profile = control$profile,
          silent = TRUE )

  if(control$quiet==FALSE) list_parameters(obj)
  # bundle
  internal = list(
    sem = sem,
    tsdata = tsdata,
    family = family,
    estimate_delta0 = estimate_delta0,
    control = control,
    covs = covs
  )

  # Further bundle
  out = list( "obj" = obj,
              "ram" = ram,
              "sem_full" = model,
              "tmb_inputs"=list("parameters"=Params, "random"=Random, "map"=Map),
              #"call" = match.call(),
              "internal" = internal )

  # Export stuff
  if( control$run_model==FALSE ){
    return( out )
  }

  # Fit
  #out$opt = fit_tmb( obj,
  #                   quiet = control$quiet,
  #                   control = list(eval.max=10000, iter.max=10000, trace=ifelse(control$quiet==TRUE,0,1) ),
  #                   ... )

  # Optimize
  out$opt = list( "par"=obj$par )
  for( i in seq_len(max(0,control$nlminb_loops)) ){
    if( isFALSE(control$quiet) ) message("Running nlminb_loop #", i)
    out$opt = nlminb( start = out$opt$par,
                  objective = obj$fn,
                  gradient = obj$gr,
                  upper = control$upper,
                  lower = control$lower,
                  control = list( eval.max = control$eval.max,
                                  iter.max = control$iter.max,
                                  trace = control$trace ) )
  }

  # Newtonsteps
  for( i in seq_len(max(0,control$newton_loops)) ){
    if( isFALSE(control$quiet) ) message("Running newton_loop #", i)
    g = as.numeric( obj$gr(out$opt$par) )
    h = optimHess(out$opt$par, fn=obj$fn, gr=obj$gr)
    out$opt$par = out$opt$par - solve(h, g)
    out$opt$objective = obj$fn(out$opt$par)
  }
  out$internal$parhat = obj$env$parList()

  if( isTRUE(control$extra_convergence_checks) ){
    # Gradient checks
    Grad_fixed = obj$gr( out$opt$par )
    if( isTRUE(any(Grad_fixed > 0.01)) ){
      warning("Some gradients are higher than 0.01. Some parameters might not be converged.  Consider increasing `control$newton_loops`")
    }
    # Hessian check ... condition and positive definite
    Hess_fixed = optimHess( par=out$opt$par, fn=obj$fn, gr=obj$gr, control=list(ndeps=rep(0.001,length(out$opt$par))) )
    Eigen_fixed = eigen( Hess_fixed, only.values=TRUE )
    if( (max(Eigen_fixed$values)/min(Eigen_fixed$values)) > 1e6 ){
      # See McCullough and Vinod 2003
      warning("The ratio of maximum and minimum Hessian eigenvalues is high. Some parameters might not be identifiable.")
    }
    if( isTRUE(any(Eigen_fixed$values < 0)) ){
      warning("Some Hessian eigenvalue is negative. Some parameters might not be converged.")
    }
  }else{
    Hess_fixed = NULL
  }

  # Run sdreport
  if( isTRUE(control$getsd) ){
    if( isTRUE(control$verbose) ) message("Running sdreport")
    if( is.null(Hess_fixed) ){
      Hess_fixed = optimHess( par=out$opt$par, fn=obj$fn, gr=obj$gr, control=list(ndeps=rep(0.001,length(out$opt$par)))  )
    }
    out$sdrep = TMB::sdreport( obj,
                          par.fixed = out$opt$par,
                          hessian.fixed = Hess_fixed,
                          getJointPrecision = control$getJointPrecision )
  }else{
    out$sdrep = NULL
  }

  # output
  class(out) = "dsem"
  return(out)
}
