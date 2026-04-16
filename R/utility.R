
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
    parameters[c("beta_z","lnsigma_z","mu_j","delta0_j")] = parlist[c("beta_z","lnsigma_z","mu_j","delta0_j")]
    for( v in c("beta_z","lnsigma_z","mu_j","delta0_j") ){
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
#' @description 
#' Calculate a data frame of total effects, resulting from a pulse experiment
#' (i.e., an exogenous and temporary change in a single variable in time \code{t=0}) or 
#' a press experiment (i.e., an exogenous and permanent change in a single variable 
#' starting in time \code{t=0} and continuing for \code{n_lags} times), representing the 
#' estimated effect of a change in any variable on every other variable and any time-lag
#' from 0 (simultaneous effects) to a user-specified maximum lag.
#'
#' @param object Output from \code{\link{dsem}}
#' @param n_lags Number of lags over which to calculate total effects
#' @param type Whether a pulse or press experiment is intended.  A pulse experiment
#' answers the question:  ``What happens if a variable is changed for only a single time-interval?"
#' A press experiment answers the question:  ``What happens if a variable is permanently changed
#' starting in a given time-interval? 
#'
#' @details
#' Total effects are taken from the Leontief matrix \eqn{\mathbf{(I-P)^{-1}}},
#' where \eqn{\mathbf{P}} is the path matrix across variables and times. 
#' \eqn{\mathbf{(I-P)}^{-1} \mathbf{\delta} }
#' calculates the effect of a perturbation represented by vector \eqn{\mathbf{\delta}}
#' with length \eqn{n_{\mathrm{lags}} \times n_{\mathrm{J}}} where \eqn{n_{\mathrm{J}}} is the number of variables.  
#' \eqn{\mathbf{(I-P)}^{-1} \mathbf{\delta} } calculates the total effect of   
#' a given variable (from)
#' upon any other variable (to) either in the same time (\eqn{t=0}), or subsequent times
#' (\eqn{t \geq 1}), where \eqn{\mathbf{\delta} = \mathbf{i}_{\mathrm{T}} \otimes \mathbf{i}_{\mathrm{J}}}, 
#' where \eqn{\mathbf{i}_{\mathrm{J}}} is one for the \code{from} variable and zero otherwise.
#' For a pulse experiment, \eqn{\mathbf{i}_{\mathrm{T}}} is one at \eqn{t=0} and zero for other times.
#' For a press experiment, \eqn{\mathbf{i}_{\mathrm{T}}} is one for all times.  
#' 
#' We compute and list the total effect at each time from time \code{t=0}
#' to \code{t=n_lags-1}.  For press experiments, this includes transient values as the the total effect 
#' approaches its asymptotic value (if this exists) as \eqn{t} approaches infinity.
#' If the analyst wants an asymptotic effect from a press experiment, we recommend
#' using a high lag (e.g., \code{n_lags = 100}) and then confirming that it has
#' reached it's asymptote (i.e., the total effect is almost identical for the last 
#' and next-to-last lag), and then reporting the value for that last lag. 
#'
#' @return
#' A data frame listing the time-lag (lag), variable that is undergoing some 
#' exogenous change (from), and the variable being impacted (to), along with the 
#' total effect (total_effect) including direct and indirect pathways, and the
#' partial "direct" effect (direct_effect)
#'
#' @examples
#' ### EXAMPLE 1
#' # Define linear model with slope of 0.5
#' sem = "
#'   # from, to, lag, name, starting_value
#'   x -> y, 0, slope, 0.5
#' "
#' # Build DSEM with specified value for path coefficients
#' mod = dsem(
#'   sem = sem,
#'   tsdata = ts(data.frame(x=rep(0,20),y=rep(0,20))),
#'   control = dsem_control( run_model = FALSE )
#' )
#'
#' # Show that total effect of X on Y from pulse experiment is 0.5 but does not propagate over time
#' pulse = total_effect(mod, n_lags = 2, type = "pulse")
#' subset( pulse, from=="x" & to=="y")
#'
#'
#' ### EXAMPLE 2
#' # Define linear model with slope of 0.5 and autocorrelated response
#' sem = "
#'   x -> y, 0, slope, 0.5
#'   y -> y, 1, ar_y, 0.8
#' "
#' mod = dsem(
#'   sem = sem,
#'   tsdata = ts(data.frame(x=rep(0,20),y=rep(0,20))),
#'   control = dsem_control( run_model = FALSE )
#' )
#'
#' # Show that total effect of X on Y from pulse experiment  is 0.5 with decay of 0.8 for each time
#' pulse = total_effect(mod, n_lags = 4, type = "pulse")
#' subset( pulse, from=="x" & to=="y")
#'
#' # Show that total effect of X on Y from press experiment  asymptotes at 2.5
#' press = total_effect(mod, n_lags = 50, type = "press")
#' subset( press, from=="x" & to=="y")
#'
#' @export
total_effect <-
function( object,
          n_lags = 4,
          type = c("pulse","press") ){

  #
  type = match.arg(type)
  
  # Unpack stuff
  Z = object$internal$tsdata
  if(is.null(object$internal$parhat)){
    object$internal$parhat = object$obj$env$parList()
  }
  # Extract path matrix
  P_kk = make_matrices(
    beta_p = object$internal$parhat$beta,
    model = object$sem_full,
    times = seq_len(n_lags),
    variables = colnames(Z)
  )$P_kk            

  # Define innovations
  if( type == "pulse" ){
    delta_kj = kronecker( Diagonal(n=ncol(Z)), 
                          sparseMatrix(i=1, j=1, x=1, dims=c(n_lags,1)) )
  }
  if( type == "press" ){
    delta_kj = kronecker( Diagonal(n=ncol(Z)), 
                          sparseMatrix(i=seq_len(n_lags), j=rep(1,n_lags), x=rep(1,n_lags), dims=c(n_lags,1)) )
  }
  IminusRho_kk = Diagonal(n=nrow(P_kk)) - P_kk
  
  # Calculate partial effect
  Partial_kj = P_kk %*% delta_kj

  # Calculate total effect using sparse matrices
  Total_kj = solve( IminusRho_kk, delta_kj )
  
  # Make into data frame
  out = expand.grid( "lag" = seq_len(n_lags)-1, 
                     "to" = colnames(Z), 
                     "from" = colnames(Z) )
  out$total_effect = as.vector(Total_kj)
  out$direct_effect = as.vector(Partial_kj)
  return(out)
}

#' @title Partition variance in one variable due to another (EXPERIMENTAL)
#'
#' @description
#' Calculate the proportion of variance for a response variable that is
#' attributed to another set of predictor variables, calculated across lags from
#' from 0 (simultaneous effects) to a user-specified maximum lag.
#'
#' @param object Output from \code{\link{dsem}}
#' @param which_response string matching colnames from \code{tsdata}
#' identifying response variable
#' @param n_times Number of lags over which to calculate total effects
#'
#' @details
#' This function calculates the variance for each variable and lag, and then
#' recalculates it when setting exogenous variance to zero for all variables except
#' \code{which_pred}.  It then calculates the ratio of the diagonal of these two.
#' This represents the proportion of variance in the full model that is attributable
#' to one or more variables.
#'
#' This function is under development and may still change or be removed.
#'
#' @return
#' A list with two elements:
#' \describe{
#'  \item{total_variance}{A matrix of the total variance for each variable (column)
#'    and each time from 1 to \code{n_times}}
#'  \item{proportion_variance_explained}{A matrix of the proportion of variance
#'    explained for variable \code{which_response} by each model variable
#'    (column) and each time from 1 to \code{n_times}}
#' }
#' Note that in a model with lagged effects, the total_variance and variance_explained
#' will vary for each time (row), and the analyst might want to either choose a time
#' for which the value has stabilized.
#'
#' @examples
#' # Simulate linear model
#' x = rnorm(100)
#' y = 1 + 1 * x + rnorm(100)
#' data = data.frame(x=x, y=y)
#'
#' # Fit as DSEM
#' fit = dsem( sem = "x -> y, 0, beta",
#'             tsdata = ts(data),
#'             control = dsem_control(quiet=TRUE) )
#'
#' # Apply
#' partition_variance( fit,
#'                     which_response = "y",
#'                     n_times = 10 )
#'
#' @export
partition_variance <-
function( object,
          which_response,
          n_times = 10 ){

  # Unpack stuff
  Z = object$internal$tsdata
  if(is.null(object$internal$parhat)){
    object$internal$parhat = object$obj$env$parList()
  }

  # Error checks
  if( !(which_response %in% colnames(Z)) ){
    stop("`which_response` not found in colnames of `tsdata`")
  }

  # Extract path matrix
  matrices = make_matrices(
    beta_p = object$internal$parhat$beta,
    model = object$sem_full,
    times = seq_len(n_times),
    variables = colnames(Z)
  )
  out = expand.grid(lag = seq_len(n_times), variable = colnames(Z) )

  #
  IminusP_kk = matrices$IminusP_kk
  invIminusP_kk = Matrix::solve(IminusP_kk)

  # Extract variance for fitted model
  G_kk = matrices$G_kk
  V_kk = Matrix::t(G_kk) %*% G_kk
  Sigma0_kk = invIminusP_kk %*% V_kk %*% Matrix::t(invIminusP_kk)

  # Zero out variances
  match_vals = which( out$variable %in% which_response )
  prop_tj = array(NA, dim=c(n_times,ncol(Z)), dimnames=list(paste0("t_",seq_len(n_times)),colnames(Z)))
  for( which_pred in colnames(Z) ){
    match_cols = which( out$variable %in% which_pred )
    G0_kk = matrices$G_kk
    G0_kk[,-match_cols] = 0
    V0_kk = Matrix::t(G0_kk) %*% G0_kk
    Sigma1_kk = invIminusP_kk %*% V0_kk %*% Matrix::t(invIminusP_kk)
    prop_tj[,which_pred] = Matrix::diag(Sigma1_kk)[match_vals] / Matrix::diag(Sigma0_kk)[match_vals]
  }

  #
  var_tj = prop_tj
  var_tj[] = Matrix::diag( Sigma0_kk )
  out = list( "total_variance" = var_tj,
              "proportion_variance_explained" = prop_tj )
  return(out)
}

#' @title
#' Classical Runge-Kutta for system of equations
#'
#' @description
#' Classical Runge-Kutta of order 4.
#'
#' @param f function in the differential equation \eqn{y' = f(x, y)};
#'        defined as a function \eqn{R \times R^m \rightarrow R^m}, where \eqn{m} is the number of equations.
#' @param a starting time for the interval to integrate
#' @param b ending time for the interval to integrate.
#' @param y0 starting values at time \code{a}
#' @param n the number of steps from a to b.
#' @param Pars named list of parameters passed to f
#' @param ... additional inputs to function \code{f}
#'
#' @details
#' Classical Runge-Kutta of order 4 for (systems of) ordinary differential
#' equations with fixed step size.
#' Copied from pracma under GPL-3, with small modifications to work with RTMB
#'
#' @return
#' List with components x for grid points between a and b and y an
#' n-by-m matrix with solutions for variables in columns, i.e.
#' each row contains one time stamp.
#'
#' @export
rk4sys <-
function( f,
          a,
          b,
          y0,
          n,
          Pars,
          ... ) {
  m <- length(y0)
  h <- (b - a)/n
  x <- seq(a + h, b, by = h)
  y <- matrix(0, nrow = n, ncol = m)
  k1 <- h * f(a, y0, Pars, ...)
  k2 <- h * f(a + h/2, y0 + k1/2, Pars, ...)
  k3 <- h * f(a + h/2, y0 + k2/2, Pars, ...)
  k4 <- h * f(a + h, y0 + k3, Pars, ...)
  y[1, ] <- y0 + k1/6 + k2/3 + k3/3 + k4/6
  for (i in seq_len(n-1)) {
    k1 <- h * f(x[i], y[i, ], Pars, ...)
    k2 <- h * f(x[i] + h/2, y[i, ] + k1/2, Pars, ...)
    k3 <- h * f(x[i] + h/2, y[i, ] + k2/2, Pars, ...)
    k4 <- h * f(x[i] + h, y[i, ] + k3, Pars, ...)
    y[i + 1, ] <- y[i, ] + k1/6 + k2/3 + k3/3 + k4/6
  }
  x <- c(a, x)
  y <- rbind(y0, y)
  return(list(x = x, y = y))
}



