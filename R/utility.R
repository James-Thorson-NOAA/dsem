
#' @title Calculate leave-one-out residuals
#'
#' @description Calculates quantile residuals using the predictive distribution from
#'              a jacknife (i.e., leave-one-out predictive distribution)
#'
#' @param object Output from \code{\link{dsem}}
#' @param nsim Number of simulations to use if \code{family!="fixed"} for some variable,
#'        such that simulation residuals are required.
#' @param what whether to return quantile residuals, or samples from the leave-one-out predictive
#'        distribution of data, or a table of leave-one-out predictions and standard errors for the
#'        latent state
#' @param track_progress whether to track runtimes on terminal
#' @param ... Not used
#'
#' @details
#' Conditional quantile residuals cannot be calculated when using \code{family = "fixed"}, because
#' state-variables are fixed at available measurements and hence the conditional distribution is a Dirac
#' delta function.  One alternative is to use leave-one-out residuals, where we calculate the predictive distribution
#' for each state value when dropping the associated observation, and then either use that as the
#' predictive distribution, or sample from that predictive distribution and then calculate
#' a standard quantile distribution for a given non-fixed family.  This appraoch is followed here.
#' It is currently only implemented when  all variables follow \code{family = "fixed"}, but
#' could be generalized to a mix of families upon request.
#'
#' @return
#' A matrix of residuals, with same order and dimensions as argument \code{tsdata}
#' that was passed to \code{dsem}.
#'
#' @export
loo_residuals <-
function( object,
          nsim = 100,
          what = c("quantiles","samples","loo"),
          track_progress = TRUE,
          ... ){

  # Extract and make object
  what = match.arg(what)
  tsdata = object$internal$tsdata
  parlist = object$obj$env$parList()
  df = expand.grid( time(tsdata), colnames(tsdata) )
  df$obs = as.vector(tsdata)
  df = na.omit(df)
  df = cbind( df, est=NA, se=NA )

  # Loop through observations
  for(r in 1:nrow(df) ){
    if( (r %% floor(nrow(df)/10))==1 ){
      if(isTRUE(track_progress)) message("Running leave-one-out fit ",r," of ",nrow(df)," at ",Sys.time())
    }
    ts_r = tsdata
    ts_r[match(df[r,1],time(tsdata)),match(df[r,2],colnames(tsdata))] = NA

    #
    control = object$internal$control
    control$quiet = TRUE
    control$run_model = FALSE

    # Build inputs
    #fit_r = dsem( sem = object$internal$sem,
    #                   tsdata = ts_r,
    #                   family = object$internal$family,
    #                   estimate_delta0 = object$internal$estimate_delta0,
    #                   control = control )
    fit_r = list( tmb_inputs = object$tmb_inputs )
    Data = fit_r$tmb_inputs$data
    fit_r$tmb_inputs$map$x_tj = factor(ifelse( is.na(as.vector(ts_r)) | (Data$familycode_j[col(tsdata)] %in% c(1,2,3,4)), seq_len(prod(dim(ts_r))), NA ))

    # Modify inputs
    map = fit_r$tmb_inputs$map
    parameters = fit_r$tmb_inputs$parameters
    parameters[c("beta_z","lnsigma_j","mu_j","delta0_j")] = parlist[c("beta_z","lnsigma_j","mu_j","delta0_j")]
    for( v in c("beta_z","lnsigma_j","mu_j","delta0_j") ){
      map[[v]] = factor( rep(NA,length(as.vector(parameters[[v]]))))
    }

    # Build object
    obj = TMB::MakeADFun( data = fit_r$tmb_inputs$data,
                     parameters = parameters,
                     random = NULL,
                     map = map,
                     profile = control$profile,
                     DLL="dsem",
                     silent = TRUE )

    # Rerun and record
    opt = nlminb( start = obj$par,
                  objective = obj$fn,
                  gradient = obj$gr )
    sdrep = TMB::sdreport( obj )
    df[r,'est'] = as.list(sdrep, what="Estimate")$x_tj[match(df[r,1],time(tsdata)),match(df[r,2],colnames(tsdata))]
    df[r,'se'] = as.list(sdrep, what="Std. Error")$x_tj[match(df[r,1],time(tsdata)),match(df[r,2],colnames(tsdata))]
  }

  # Compute quantile residuals
  resid_tjr = array( NA, dim=c(dim(tsdata),nsim) )
  if( all(object$internal$family == "fixed") ){
    # analytical quantile residuals
    resid_tj = array( NA, dim=dim(tsdata) )
    resid_tj[cbind(match(df[,1],time(tsdata)),match(df[,2],colnames(tsdata)))] = pnorm( q=df$obs, mean=df$est, sd=df$se )
    # Simulations from predictive distribution for use in DHARMa etc
    for(z in 1:nrow(df) ){
      resid_tjr[match(df[z,1],time(tsdata)),match(df[z,2],colnames(tsdata)),] = rnorm( n=nsim, mean=df[z,'est'], sd=df[z,'se'] )
    }
  }else{
    # Sample from leave-one-out predictive distribution for states
    resid_rz = apply( df, MARGIN=1, FUN=function(vec){rnorm(n=nsim, mean=as.numeric(vec['est']), sd=as.numeric(vec['se']))} )
    # Sample from predictive distribution of data given states
    for(r in 1:nrow(resid_rz) ){
      parameters = object$obj$env$parList()
      parameters$x_tj[which(!is.na(tsdata))] = resid_rz[r,]
      obj = TMB::MakeADFun( data = object$tmb_inputs$data,
                       parameters = parameters,
                       random = NULL,
                       map = object$tmb_inputs$map,
                       profile = NULL,
                       DLL="dsem",
                       silent = TRUE )
      resid_tjr[,,r] = obj$simulate()$y_tj
    }
    # Calculate quantile residuals
    resid_tj = apply( resid_tjr > outer(tsdata,rep(1,nsim)), MARGIN=1:2, FUN=mean )
  }
  if(what=="quantiles") return( resid_tj )
  if(what=="samples") return( resid_tjr )
  if(what=="loo") return( df )
}

#' @title Calculate total effects
#'
#' @description Calculate a data frame of total effects, representing the 
#' estimated effect of every variable on every other variable and any time-lag
#' from 0 (simultaneous effects) to a user-specified maximum lag.
#'
#' @param object Output from \code{\link{dsem}}
#' @param n_lags Number of lags over which to calculate total effects
#'
#' @details
#' Total effects are taken from the Leontief matrix \eqn{\mathbf{(I-P)^{-1}}},
#' where \eqn{\mathbf{P}} is the path matrix across variables and times.
#'
#' @return
#' A data frame listing the time-lag (lag), variable that is undergoing some 
#' exogenous change (from), and the variable being impacted (to), along with the 
#' total effect (effect) including direct and indirect pathways
#'
#' @export
total_effect <-
function( object,
          n_lags = 4 ){

  # Extract path matrix
  Z = object$internal$tsdata
  P_kk = make_matrices( 
    beta_p = object$internal$parhat$beta,
    model = object$sem_full,
    times = seq_len(n_lags),
    variables = colnames(Z)
  )$P_kk            

  # Define innovations
  delta_kj = kronecker( Diagonal(n=ncol(Z)), 
                        sparseMatrix(i=1, j=1, x=1, dims=c(n_lags,1)) )
  IminusRho_kk = Diagonal(n=nrow(P_kk)) - P_kk
  
  # Calculate total effect using sparse matrices
  Impact_kj = solve( IminusRho_kk, delta_kj )
  
  # Make into data frame
  out = expand.grid( "lag" = seq_len(n_lags)-1, 
                     "to" = colnames(Z), 
                     "from" = colnames(Z) )
  out$effect = as.vector(Impact_kj)
  return(out)
}

