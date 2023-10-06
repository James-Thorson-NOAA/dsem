#' @title Fit dynamic structural equation model
#'
#' @description Fits a dynamic structural equation model
#'
#' @param sem structural equation model structure, passed to either \code{\link[sem]{specifyModel}}
#'        or \code{\link[sem]{specifyEquations}} and then parsed to control
#'        the set of path coefficients and variance-covariance parameters.
#' @param tsdata time-series data, as outputted using \code{\link[stats]{ts}}
#' @param family Character-vector listing the distribution used for each column of \code{tsdata}, where
#'        each element must be \code{fixed} or \code{normal}.
#'        \code{family="fixed"} is default behavior and assumes that a given variable is measured exactly.
#'        Other options correspond to different specifications of measurement error.
#' @param estimate_delta0 Boolean indicating whether to estimate deviations from equilibrium in initial year
#'        as fixed effects, or alternatively to assume that dynamics start at some stochastic draw away from
#'        the stationary distribution
#' @param run_model Boolean indicating whether to estimate parameters (the default), or
#'        instead to return the model inputs and compiled TMB object without running;
#' @param quiet Boolean indicating whether to run model printing messages to terminal or not;
#' @param use_REML Boolean indicating whether to treat non-variance fixed effects as random,
#'        either to motigate bias in estimated variance parameters or improve efficiency for
#'        parameter estimation given correlated fixed and random effects
#' @param parameters list of fixed and random effects, e.g., as constructed by \code{dsem} and then modified
#'        by hand (only helpful for advanced users to change starting values or restart at intended values)
#' @param map list of fixed and mirrored parameters, constructed by \code{dsem} by default but available
#'        to override this default and then pass to \code{\link[TMB]{MakeADFun}}
#' @param ... Additional parameters passed to \code{\link{fit_tmb}}
#'
#' @importFrom TMB compile dynlib MakeADFun sdreport summary.sdreport
#' @importFrom stats .preformat.ts na.omit nlminb optimHess pnorm rnorm
#'
#' @details
#' A DSEM involves at a minimum a user-supplied specification for the path coefficients, and a
#' time-series with named variables that correspond to the user-supplied paths. Users can specify
#' a distribution for measurement errors (or assume that variables are measured without error) using
#' argument \code{family}.
#' \code{dsem} then estimates all specified coefficients via maximizing a log-marginal likelihood, while
#' also imputing the value of any variable that is specified as NA in the time-series inputs.
#' \code{summary.dsem} then assembles estimates and standard errors in an easy-to-read format.
#' Standard errors for fixed effects (path coefficients, exogenoux variance parameters, and measurement error parameters)
#' are estimated from the matrix of second derivatives of the log-marginal likelihod,
#' and standard errors for random effects (i.e., missing or state-space variables) are estimated
#' from a generalization of this method (see \code{?TMB::sdreport} for details).
#'
#' @return
#' An object (list) of class `dsem`. Elements include:
#' \describe{
#' \item{obj}{TMB object from \code{\link[TMB]{MakeADFun}}}
#' \item{ram}{RAM parsed by \code{make_ram}}
#' \item{model}{SEM model parsed from \code{sem} using \code{\link[sem]{specifyModel}} or \code{\link[sem]{specifyEquations}}}
#' \item{tmb_inputs}{The list of inputs passed to \code{\link[TMB]{MakeADFun}}}
#' \item{opt}{The output from \code{\link{fit_tmb}}}
#' }
#'
#' @references
#'
#' **Introducing the package, its features, and comparison with other software
#' (to cite when using dsem):**
#'
#' Thorson, J. T., Andrews, A., Essington, T., Large, S. (In review).
#' Dynamic structural equation models synthesize
#' ecosystem dynamics constrained by ecological mechanisms.
#'
#' @examples
#' # Define model
#' sem = "
#'   Profits -> Consumption, 0, a2
#'   Profits -> Consumption, -1, a3
#'   Priv_wage -> Consumption, 0, a4
#'   Gov_wage -> Consumption, 0, a4
#'   Consumption -> Consumption, -1, ar1
#'   Consumption -> Consumption, -2, ar2
#'   Profits -> Investment, 0, b2
#'   Profits -> Investment, -1, b3
#'   Capital_stock -> Investment, -1, b4
#'   GNP -> Priv_wage, 0, c2
#'   Taxes -> Priv_wage, 0, c2
#'   neg_Gov_wage -> Priv_wage, 0, c2
#'   GNP -> Priv_wage, -1, c3
#'   Taxes -> Priv_wage, -1, c3
#'   neg_Gov_wage -> Priv_wage, -1, c3
#'   Time -> Priv_wage, 0, c4
#' "
#'
#' # Load data
#' data(KleinI, package="AER")
#' Data = as.data.frame(KleinI)
#' Data = cbind( Data, "time" = seq(1,22)-11 )
#' colnames(Data) = sapply( colnames(Data), FUN=switch,
#'            "consumption"="Consumption", "invest"="Investment",
#'            "cprofits"="Profits", "capital"="Capital_stock", "gwage"="Gov_wage",
#'            "pwage"="Priv_wage", "gexpenditure"="Gov_expense", "taxes"="Taxes",
#'            "time"="Time", "gnp"="GNP")
#' Z = ts( cbind(Data, "neg_Gov_wage"=-1*Data[,'Gov_wage']) )
#'
#' # Fit model
#' fit = dsem( sem=sem, tsdata=Z )
#' summary( fit )
#'
#' # Plot results
#' library(ggplot2)
#' library(ggpubr)
#' library(phylopath)
#' p1 = plot(as_fitted_DAG(fit), text_size=3, type="width", show.legend=FALSE)
#' p1$layers[[1]]$mapping$edge_width = 0.5
#' p2 = plot(as_fitted_DAG(fit, lag=-1), text_size=3, type="width", show.legend=FALSE)
#' p2$layers[[1]]$mapping$edge_width = 0.25
#' ggarrange(p1 + scale_x_continuous(expand = c(0.2, 0.0)),
#'                     p2 + scale_x_continuous(expand = c(0.2, 0.0)),
#'                     labels = c("Simultaneous effects", "Lag-1 effects"),
#'                     ncol = 1, nrow = 2)
#'
#' @useDynLib dsem, .registration = TRUE
#' @export
dsem <-
function( sem,
          tsdata,
          family = rep("fixed",ncol(tsdata)),
          covs = colnames(tsdata),
          estimate_delta0 = FALSE,
          quiet = FALSE,
          run_model = TRUE,
          use_REML = TRUE,
          parameters = NULL,
          map = NULL,
          ... ){

  # (I-Rho)^-1 * Gamma * (I-Rho)^-1
  out = make_ram( sem, tsdata=tsdata, quiet=quiet, covs=covs )
  ram = out$ram

  # Error checks
  if( any((out$model[,'direction']==2) & (out$model[,2]!=0)) ){
    stop("All two-headed arrows should have lag=0")
  }
  if( !all(c(out$model[,'first'],out$model[,'second']) %in% colnames(tsdata)) ){
    stop("Some variable in `sem` is not in `tsdata`")
  }

  #
  Data = list( "RAM" = as.matrix(na.omit(ram[,1:4])),
               "RAMstart" = as.numeric(ram[,5]),
               "familycode_j" = sapply(family, FUN=switch, "fixed"=0, "normal"=1 ),
               "y_tj" = tsdata )

  # Construct parameters
  if( is.null(parameters) ){
    Params = list( "beta_z" = rep(0,max(ram[,4])),
                   "lnsigma_j" = rep(0,ncol(tsdata)),
                   "mu_j" = rep(0,ncol(tsdata)),
                   "delta0_j" = rep(0,ncol(tsdata)),
                   "x_tj" = ifelse( is.na(tsdata), 0, tsdata ))

    # Turn off initial conditions
    if( estimate_delta0==FALSE ){
      Params$delta0_j = numeric(0)
    }

    # Scale starting values with higher value for two-headed than one-headed arrows
    which_nonzero = which(ram[,4]>0)
    beta_type = tapply( ram[which_nonzero,1], INDEX=ram[which_nonzero,4], max)
    Params$beta_z = ifelse(beta_type==1, 0.01, 1)

    # Override starting values if supplied
    which_nonzero = which(ram[,4]>0)
    start_z = tapply( as.numeric(ram[which_nonzero,5]), INDEX=ram[which_nonzero,4], mean )
    Params$beta_z = ifelse( is.na(start_z), Params$beta_z, start_z)
  }else{
    Params = parameters
  }

  # Construct map
  if( is.null(map) ){
    Map = list()
    Map$x_tj = factor(ifelse( is.na(as.vector(tsdata)) | (Data$familycode_j[col(tsdata)] %in% c(1,2,3,4)), seq_len(prod(dim(tsdata))), NA ))
    Map$lnsigma_j = factor( ifelse(Data$familycode_j==0, NA, seq_along(Params$lnsigma_j)) )

    # Map off mean for latent variables
    Map$mu_j = factor( ifelse(colSums(!is.na(tsdata))==0, NA, 1:ncol(tsdata)) )
  }else{
    Map = map
  }

  # Initial run
  if(isTRUE(use_REML)){
    Random = c( "x_tj", "mu_j" )
  }else{
    Random = "x_tj"
  }
  obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL="dsem" )
  if(quiet==FALSE) list_parameters(obj)
  out = list( "obj"=obj,
              "ram"=ram,
              "sem_full"=out$model,
              "tmb_inputs"=list("data"=Data, "parameters"=Params, "random"=Random, "map"=Map),
              "call" = match.call() )

  # Export stuff
  if( run_model==FALSE ){
    return( out )
  }

  # Fit
  obj$env$beSilent()       # if(!is.null(Random))
  out$opt = fit_tmb( obj,
                     quiet = quiet,
                     control = list(eval.max=10000, iter.max=10000, trace=ifelse(quiet==TRUE,0,1) ),
                     ... )

  # output
  class(out) = "dsem"
  return(out)
}

#' @title summarize dsem
#'
#' @description summarize parameters from a fitted dynamic structural equation model
#'
#' @details
#' A DSEM is specified using "arrow and lag" notation, which specifies the set of
#' path coefficients and exogenous variance parameters to be estimated. Function \code{dsem}
#' then estimates the maximum likelihood value for those coefficients and parameters
#' by maximizing the log-marginal likelihood.  Standard errors for parameters are calculated
#' from the matrix of second derivatives of this log-marginal likelihood (the "Hessian matrix").
#'
#' However, many users will want to associate individual parameters and standard errors
#' with the path coefficients that were specified using the "arrow and lag" notation.
#' This task is complicated in
#' models where some path coefficients or variance parameters are specified to share a single value a priori,
#' or were assigned a name of NA and hence assumed to have a fixed value a priori (such that
#' these coefficients or parameters have an assigned value but no standard error).
#' The \code{summary} function therefore compiles the MLE for coefficients (including duplicating
#' values for any path coefficients that assigned the same value) and standard error
#' estimates, and outputs those in a table that associates them with the user-supplied path and parameter names.
#' It also outputs the z-score and a p-value arising from a two-sided Wald test (i.e.
#' comparing the estimate divided by standard error against a standard normal distribution).
#'
#' @param object Output from \code{\link{dsem}}
#' @param ... Not used
#'
#' @method summary dsem
#' @export
summary.dsem <-
function( object, ... ){

  # Easy of use
  model = object$sem_full
  ParHat = object$obj$env$parList()

  #
  coefs = data.frame( model, "Estimate"=c(NA,ParHat$beta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
  coefs$Estimate = ifelse( is.na(coefs$Estimate), as.numeric(model[,4]), coefs$Estimate )
  if( "SD" %in% names(object$opt) ){
    SE = as.list( object$opt$SD, report=FALSE, what="Std. Error")
    coefs = data.frame( coefs, "Std_Error"=c(NA,SE$beta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
    coefs = data.frame( coefs, "z_value"=coefs[,'Estimate']/coefs[,'Std_Error'] )
    coefs = data.frame( coefs, "p_value"=pnorm(-abs(coefs[,'z_value'])) * 2 )
  }

  return(coefs)
}

#' @title Simulate dsem
#'
#' @description Simulate from a fitted \code{dsem} model
#'
#' @param object Output from \code{\link{dsem}}
#' @param nsim number of simulated data sets
#' @param variance whether to ignore uncertainty in fixed and
#'        random effects, include estimation uncertainty in random effects,
#'        or include estimation uncertainty in both fixed and random effects
#' @param resimulate_gmrf whether to resimulate the GMRF based on estimated or
#'        simulated random effects (determined by argument \code{variance})
#' @param seed random seed
#' @param ... Not used
#'
#' @details
#' This function conducts a parametric bootstrap, i.e., simulates new data
#' conditional upon estimated values for fixed and random effects.  The user
#' can optionally simulate new random effects conditional upon their estimated
#' covariance, or simulate new fixed and random effects conditional upon their imprecision.
#'
#' Note that \code{simulate} will have no effect on states \code{x_tj} for which there
#' is a measurement and when those measurements are fitted using \code{family="fixed"}, unless
#' \code{resimulate_gmrf=TRUE}.  In this latter case, the GMRF is resimulated given
#' estimated path coefficients
#'
#' @method simulate dsem
#' @export
simulate.dsem <-
function( object,
          nsim = 1,
          seed = NULL,
          variance = c("none", "random", "both"),
          resimulate_gmrf = FALSE,
          ... ){

  # Front stuff
  set.seed(seed)
  variance = match.arg(variance)

  # Sample from GMRF using sparse precision
  rmvnorm_prec <- function(mu, prec, nsim) {
    z <- matrix(rnorm(length(mu) * nsim), ncol=nsim)
    L <- Matrix::Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    return(mu + z)
  }

  # pull out objects for easy use
  obj = object$obj
  parfull = obj$env$parList()
  tsdata = eval(object$call$tsdata)

  # Extract parameters, and add noise as desired
  par_zr = outer( obj$env$last.par.best, rep(1,nsim) )
  if( variance=="random" ){
    eps_zr = rmvnorm_prec( rep(0,length(obj$env$random)), obj$env$spHess(random=TRUE), nsim=nsim )
    par_zr[obj$env$random,] = par_zr[obj$env$random,,drop=FALSE] + eps_zr
  }
  if( variance=="both" ){
    if(is.null(object$opt$SD$jointPrecision)){
      stop("Please re-run `dsem` with `getsd=TRUE` and `getJointPrecision=TRUE`, or confirm that the model is converged")
    }
    eps_zr = rmvnorm_prec( rep(0,length(obj$env$last.par)), object$opt$SD$jointPrecision, nsim=nsim )
    par_zr = par_zr + eps_zr
  }

  # Simulate new GMRF and data conditional on simulated parameters
  out = NULL
  for( r in seq_len(nsim) ){
    if( resimulate_gmrf==TRUE ){
      # Simulate new fields
      newrep = obj$report( par=par_zr[,r] )
      newparfull = obj$env$parList()
      Q_kk = newrep$Q_kk
      tmp = rmvnorm_prec( newrep$delta_k + as.vector(newrep$xhat_tj), Q_kk, nsim=1 )
      # Modify call
      newcall = object$call
      newcall$parameters = newparfull
      newcall$parameters$x_tj[] = tmp
      # Rebuild model with new GMRF values
      newcall$run_model = FALSE
      newfit = eval(newcall)
      out[[r]] = newfit$obj$simulate()$y_tj
    }else{
      out[[r]] = obj$simulate( par_zr[,r] )
    }
    colnames(out[[r]]) = colnames(tsdata)
    tsp(out[[r]]) = attr(tsdata,"tsp")
    class(out[[r]]) = class(tsdata)
  }

  return(out)
}

#' Extract Variance-Covariance Matrix
#'
#' extract the covariance of fixed effects, or both fixed and random effects.
#'
#' @param object output from \code{dsem}
#' @param which whether to extract the covariance among fixed effects, random effects, or both
#' @param ... ignored, for method compatibility
#' @importFrom stats vcov
#' @method vcov dsem
#' @export
vcov.dsem <-
function( object,
          which = c("fixed", "random", "both"),
          ...) {

  which = match.arg(which)

  if( which=="fixed" ){
    V = object$opt$SD$cov.fixed
    if(is.null(V)){
      warning("Please re-run `dsem` with `getsd=TRUE`, or confirm that the model is converged")
    }
  }
  if( which=="random" ){
    V = solve(object$obj$env$spHess(random=TRUE))
  }
  if( which=="both" ){
    H = object$opt$SD$jointPrecision
    if(is.null(H)){
      warning("Please re-run `dsem` with `getsd=TRUE` and `getJointPrecision=TRUE`, or confirm that the model is converged")
      V = NULL
    }else{
      V = solve(H)
    }
  }

  return( V )
}

#' Calculate residuals
#'
#' @title Calculate deviance or response residuals for dsem
#'
#' @param object Output from \code{\link{dsem}}
#' @param type which type of residuals to compute (only option is \code{"deviance"} or \code{"response"} for now)
#' @param ... Note used
#'
#' @method residuals dsem
#' @export
residuals.dsem <-
function( object,
          type = c("deviance","response"),
          ... ){

  # https://stats.stackexchange.com/questions/1432/what-do-the-residuals-in-a-logistic-regression-mean
  # Normal deviance residuals
  if( FALSE ){
    x = rnorm(10)
    y = x + rnorm(10)
    Glm = glm( y ~ 1 + x, family="gaussian")
    mu = predict(Glm,type="response")
    r1 = y - mu
    r1 - resid(Glm)
  }
  # Poisson deviance residuals
  if( FALSE ){
    x = rnorm(10)
    y = rpois(10, exp(x))
    Glm = glm( y ~ 1 + x, family="poisson")
    mu = predict(Glm,type="response")
    # https://stats.stackexchange.com/questions/398098/formula-for-deviance-residuals-for-poisson-model-with-identity-link-function
    r1 = sign(y - mu) * sqrt(2*(y*log((y+1e-10)/mu) - (y-mu)))
    r1 - resid(Glm)
  }
  # Binomial deviance residuals
  if( FALSE ){
    p = 0.5
    y = rbinom(10, prob=p, size=1)
    Glm = glm( y ~ 1, family="binomial")
    mu = predict(Glm, type="response")
    r1 = sign(y - mu) * sqrt(-2*(((1-y)*log(1-mu) + y*log(mu))))
    r1 - resid(Glm)
  }
  # Gamma deviance residuals
  if( FALSE ){
    mu = 1
    cv = 0.8
    y = rgamma( n=10, shape=1/cv^2, scale=mu*cv^2 )
    Glm = glm( y ~ 1, family=Gamma(link='log'))
    mu = predict(Glm, type="response")
    r1 = sign(y - mu) * sqrt(2 * ( (y-mu)/mu - log(y/mu) ))
    r1 - resid(Glm)
  }

  # Poisson: sign(y - mu) * sqrt(2*(ifelse(y==0, 0, y*log(y/mu)) - (y-mu)))
  # Binomial:  -2 * ((1-y)*log(1-mu) + y*log(mu))
  # Gamma: 2 * ( (y-mu)/mu - log(y/mu) )

  # Easy of use
  x_tj = object$obj$env$parList()$x_tj
  y_tj = object$tmb_inputs$data$y_tj
  familycode_j = object$tmb_inputs$data$familycode_j
  report = object$obj$report()

  #
  type = match.arg(type)
  if( type == "deviance" ){
    resid_tj = report$devresid_tj
  }
  if( type == "response" ){
    resid_tj = y_tj - report$mu_tj
  }

  return(resid_tj)
}

#' predictions using dsem
#'
#' @title Predict variables given new (counterfactual) values of data, or for future or past times
#'
#' @param object Output from \code{\link{dsem}}
#' @param newdata optionally, a data frame in which to look for variables with which to predict.
#'        If omitted, the fitted data are used to create predictions. If desiring predictions after the fitted data,
#'        the user must append rows with NAs for those future times.  Similarly, if desiring predictions given counterfactual
#'        values for time-series data, then those individual observations can be edited while keeping other observations at their
#'        original fitted values.
#' @param type the type of prediction required. The default is on the scale of the linear predictors;
#'        the alternative "response" is on the scale of the response variable.
#'        Thus for a Poisson-distributed variable the default predictions are of log-intensity and type = "response" gives the predicted intensity.
#' @param ... Note used
#'
#' @method predict dsem
#' @export
predict.dsem <-
function( object,
          newdata = NULL,
          type = c("link", "response"),
          ... ){
  #
  # newdata = eval(object$call$tsdata)
  # newdata = ts( newdata[1:40,] )

  # Easy of use
  parfull = object$obj$env$parList()
  type = match.arg(type)
  report = object$obj$report()

  #
  if( is.null(newdata) ){
    if(type=="link") out = parfull$x_tj
    if(type=="response") out = report$mu_tj
  }else{
    newcall = object$call
    newcall$tsdata = newdata
    # Rebuild model with new data
    newcall$run_model = FALSE
    newfit = eval(newcall)
    # Optimize random effects given original MLE and newdata
    newfit$obj$fn( object$opt$par )
    # Return predictor
    if(type=="link") out = newfit$obj$env$parList()$x_tj
    if(type=="response") out = newfit$obj$report()$mu_tj
  }

  return(out)
}

# Extract the (marginal) log-likelihood of a dsem model
#
# @return object of class \code{logLik} with attributes
#   \item{val}{log-likelihood}
#   \item{df}{number of parameters}
#' @importFrom stats logLik
#' @export
logLik.dsem <- function(object, ...) {
  val = -1 * object$opt$objective
  df = length( object$opt$par )
  out = structure( val,
             df = df,
             class = "logLik")
  return(out)
}

#' Convert dsem to phylopath output
#'
#' @title Convert output from package dsem to phylopath
#'
#' @param fit Output from \code{\link{dsem}}
#' @param lag which lag to output
#' @param what whether to output estimates \code{what="Estimate"}, standard errors \code{what="Std_Error"}
#'        or p-values \code{what="Std_Error"}
#' @param direction whether to include one-sided arrows \code{direction=1}, or both one- and two-sided arrows \code{direction=c(1,2)}
#'
#' @return Convert output to format supplied by \code{\link[phylopath]{est_DAG}}
#'
#' @export
as_fitted_DAG <-
function( fit,
          lag = 0,
          what = "Estimate",
          direction = 1 ){

  coefs = summary( fit )
  coefs = coefs[ which(coefs[,2]==lag), ]
  coefs = coefs[ which(coefs[,'direction'] %in% direction), ]

  #
  vars = unique( c(coefs[,'first'],coefs[,'second']) )
  out = list( "coef"=array(0, dim=rep(length(vars),2), dimnames=list(vars,vars)) )
  out$coef[as.matrix(coefs[,c('first','second')])] = coefs[,what]

  class(out) = "fitted_DAG"
  return(out)
}

#' Convert dsem to sem output
#'
#' @title Convert output from package dsem to sem
#'
#' @param object Output from \code{\link{dsem}}
#'
#' @return Convert output to format supplied by \code{\link[sem]{sem}}
#'
#' @export
#as_sem <-
#function( object ){
#
#  Rho = as_fitted_DAG( object, what="Estimate", direction=1, lag=0 )$coef
#  Gamma = as_fitted_DAG( object, what="Estimate", direction=2, lag=0 )$coef
#  Gammainv = diag(1/diag(Gamma))
#  Linv = Gammainv %*% (diag(nrow(Rho))-Rho)
#  Sinv = t(Linv) %*% Linv
#  Sprime = solve(Sinv)
#  Sprime = 0.5*Sprime + 0.5*t(Sprime)
#
#  model = object$sem_full
#  model = model[model[,2]==0,c(1,3,4)]
#  out = sem::sem( model,
#             S = Sprime,
#             N = nrow(eval(object$call$tsdata)) )
#
#  # pass out
#  return(out)
#}
#
