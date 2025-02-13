#' @title Fit dynamic structural equation model
#'
#' @description Fits a dynamic structural equation model
#'
#' @param sem Specification for time-series structural equation model structure
#'        including lagged or simultaneous effects.  See Details section in
#'        \code{\link[dsem]{make_dsem_ram}} for more description
#' @param tsdata time-series data, as outputted using \code{\link[stats]{ts}}
#' @param family Character-vector listing the distribution used for each column of \code{tsdata}, where
#'        each element must be \code{fixed} (for no measurement error), 
#'        \code{normal} for normal measurement error using an identity link,
#'        \code{gamma} for a gamma measurement error using a fixed CV and log-link, 
#'        \code{bernoulli} for a Bernoulli measurement error using a logit-link, or
#'        \code{poisson} for a Poisson measurement error using a log-link.
#'        \code{family="fixed"} is default behavior and assumes that a given variable is measured exactly.
#'        Other options correspond to different specifications of measurement error.
#' @param estimate_delta0 Boolean indicating whether to estimate deviations from equilibrium in initial year
#'        as fixed effects, or alternatively to assume that dynamics start at some stochastic draw away from
#'        the stationary distribution
#' @param covs optional: a character vector of one or more elements, with each element giving a string of variable 
#'        names, separated by commas. Variances and covariances among all variables in each such string are 
#'        added to the model. Warning: covs="x1, x2" and covs=c("x1", "x2") are not equivalent: 
#'        covs="x1, x2" specifies the variance of x1, the variance of x2, and their covariance, 
#'        while covs=c("x1", "x2") specifies the variance of x1 and the variance of x2 but not their covariance.
#'        These same covariances can be added manually via argument `sem`, but using argument `covs` might
#'        save time for models with many variables.
#' @param prior_negloglike A user-provided function that takes as input the vector of fixed effects out$obj$par
#'        returns the negative log-prior probability. For example
#'        \code{prior_negloglike = function(obj) -1 * dnorm( obj$par[1], mean=0, sd=0.1, log=TRUE)}
#'        specifies a normal prior probability
#'        for the for the first fixed effect with mean of zero and logsd of 0.1.
#'        NOTE:  this implementation does not work well with \code{tmbstan} and
#'        is highly experimental.  If using priors, considering using \code{\link{dsemRTMB}}
#'        instead.  The option in \code{dsem} is mainly intended to validate its
#'        use in \code{dsemRTMB}.  Note that the user must load RTMB using
#'        \code{library(RTMB)} prior to running the model.
#' @param control Output from \code{\link{dsem_control}}, used to define user
#'        settings, and see documentation for that function for details.
#'
#' @importFrom TMB compile dynlib MakeADFun sdreport summary.sdreport
#' @importFrom stats AIC sd .preformat.ts na.omit nlminb optimHess pnorm rbinom rgamma rpois rnorm simulate time tsp<- plogis pchisq
#' @importFrom Matrix solve Cholesky sparseMatrix mat2triplet drop0 t
#' @importFrom sem sem
#' @importFrom igraph plot.igraph graph_from_data_frame with_sugiyama layout_
#' @importFrom ggraph ggraph geom_edge_arc create_layout rectangle geom_node_text theme_graph
#' @importFrom ggplot2 aes
#' @importFrom grid arrow
#' @importFrom methods is
#' @importFrom ggm basiSet findPath isAcyclic topSort
#' @importFrom utils combn
#'
#' @details
#' A DSEM involves (at a minimum):
#' \describe{
#'   \item{Time series}{a matrix \eqn{\mathbf X} where column \eqn{\mathbf x_c} for variable c is
#'         a time-series;}
#'   \item{Path diagram}{a user-supplied specification for the path coefficients, which
#'         define the precision (inverse covariance) \eqn{\mathbf Q} for a matrix of state-variables
#'         and see \code{\link{make_dsem_ram}} for more details on the math involved.}
#' }
#' The model also estimates the time-series mean \eqn{ \mathbf{\mu}_c } for each variable.
#' The mean and precision matrix therefore define a Gaussian Markov random field for \eqn{\mathbf X}:
#'
#' \deqn{ \mathrm{vec}(\mathbf X) \sim \mathrm{MVN}( \mathrm{vec}(\mathbf{I_T} \otimes \mathbf{\mu}), \mathbf{Q}^{-1}) }
#'
#' Users can the specify
#' a distribution for measurement errors (or assume that variables are measured without error) using
#' argument \code{family}.  This defines the link-function \eqn{g_c(.)} and distribution \eqn{f_c(.)}
#' for each time-series \eqn{c}:
#'
#' \deqn{ y_{t,c} \sim f_c( g_c^{-1}( x_{t,c} ), \theta_c )}
#'
#' \code{dsem} then estimates all specified coefficients, time-series means \eqn{\mu_c}, and distribution
#' measurement errors \eqn{\theta_c} via maximizing a log-marginal likelihood, while
#' also estimating state-variables \eqn{x_{t,c}}.
#' \code{summary.dsem} then assembles estimates and standard errors in an easy-to-read format.
#' Standard errors for fixed effects (path coefficients, exogenoux variance parameters, and measurement error parameters)
#' are estimated from the matrix of second derivatives of the log-marginal likelihod,
#' and standard errors for random effects (i.e., missing or state-space variables) are estimated
#' from a generalization of this method (see \code{\link[TMB]{sdreport}} for details).
#'
#' @return
#' An object (list) of class `dsem`. Elements include:
#' \describe{
#' \item{obj}{TMB object from \code{\link[TMB]{MakeADFun}}}
#' \item{ram}{RAM parsed by \code{make_dsem_ram}}
#' \item{model}{SEM structure parsed by \code{make_dsem_ram} as intermediate description of model linkages}
#' \item{tmb_inputs}{The list of inputs passed to \code{\link[TMB]{MakeADFun}}}
#' \item{opt}{The output from \code{\link[stats]{nlminb}}}
#' \item{sdrep}{The output from \code{\link[TMB]{sdreport}}}
#' \item{interal}{Objects useful for package function, i.e., all arguments
#'                passed during the call}
#' \item{run_time}{Total time to run model}
#' }
#'
#' @references
#'
#' **Introducing the package, its features, and comparison with other software
#' (to cite when using dsem):**
#'
#' Thorson, J. T., Andrews, A., Essington, T., Large, S. (2024).
#' Dynamic structural equation models synthesize
#' ecosystem dynamics constrained by ecological mechanisms.
#' Methods in Ecology and Evolution. \doi{10.1111/2041-210X.14289}
#'
#' @examples
#' # Define model
#' sem = "
#'   # Link, lag, param_name
#'   cprofits -> consumption, 0, a1
#'   cprofits -> consumption, 1, a2
#'   pwage -> consumption, 0, a3
#'   gwage -> consumption, 0, a3
#'   cprofits -> invest, 0, b1
#'   cprofits -> invest, 1, b2
#'   capital -> invest, 0, b3
#'   gnp -> pwage, 0, c2
#'   gnp -> pwage, 1, c3
#'   time -> pwage, 0, c1
#' "
#'
#' # Load data
#' data(KleinI, package="AER")
#' TS = ts(data.frame(KleinI, "time"=time(KleinI) - 1931))
#' tsdata = TS[,c("time","gnp","pwage","cprofits",'consumption',
#'                "gwage","invest","capital")]
#'
#' # Fit model
#' fit = dsem( sem=sem,
#'             tsdata = tsdata,
#'             estimate_delta0 = TRUE,
#'             control = dsem_control(quiet=TRUE) )
#' summary( fit )
#' plot( fit )
#' plot( fit, edge_label="value" )
#'
#' @useDynLib dsem, .registration = TRUE
#' @export
dsem <-
function( sem,
          tsdata,
          family = rep("fixed",ncol(tsdata)),
          estimate_delta0 = FALSE,
          prior_negloglike = NULL,
          control = dsem_control(),
          covs = colnames(tsdata) ){
  start_time = Sys.time()

  # General error checks
  if( isFALSE(is(control, "dsem_control")) ) stop("`control` must be made by `dsem_control()`")
  if( isTRUE(control$gmrf_parameterization=="projection") ){
    if( isTRUE(any(family=="fixed" & colSums(!is.na(tsdata))>0)) ){
      stop("`family` cannot be `fixed` using `gmrf_parameterization=projection` for any variable with data")
    }
  }
  if( isFALSE(is(tsdata,"ts")) ) stop("`tsdata` must be a `ts` object")

  # General warnings
  if( isFALSE(control$quiet) ){
    tsdata_SD = apply( tsdata, MARGIN=2, FUN=sd, na.rm=TRUE )
    if( isTRUE( (max(tsdata_SD,na.rm=TRUE)/min(tsdata_SD,na.rm=TRUE)) > 10) ){
       warning("Some variables in `tsdata` have much higher variance than others. Please consider rescaling variables to prevent issues with numerical convergence.")
    }
  }

  # (I-Rho)^-1 * Gamma * (I-Rho)^-1
  out = make_dsem_ram( sem,
                  times = as.numeric(time(tsdata)),
                  variables = colnames(tsdata),
                  covs = covs,
                  quiet = control$quiet )
  ram = out$ram

  # Error checks
  if( isTRUE(any((out$model[,'direction']==2) & (out$model[,2]!=0))) ){
    stop("All two-headed arrows should have lag=0")
  }
  if( isFALSE(all(c(out$model[,'first'],out$model[,'second']) %in% colnames(tsdata))) ){
    stop("Some variable in `sem` is not in `tsdata`")
  }
  if( isFALSE(ncol(tsdata) == length(unique(colnames(tsdata)))) ){
    stop("Please check `colnames(tsdata)` to confirm that all variables (columns) have a unique name")
  }
  if( isFALSE(all(control$lower == -Inf)) | isFALSE(all(control$upper == Inf)) ){
    if( isTRUE(control$newton_loops > 0) ){
      stop("If specifying `lower` or `upper`, please set `dsem_control('newton_loops'=0)`")
    }
  }

  #
  options = c(
    ifelse(control$gmrf_parameterization=="separable", 0, 1),
    switch(control$constant_variance, "conditional"=0, "marginal"=1, "diagonal"=2)
  )
  
  #
  Data = list( "options" = options,
               "RAM" = as.matrix(na.omit(ram[,1:4])),
               "RAMstart" = as.numeric(ram[,5]),
               "familycode_j" = sapply(family, FUN=switch, "fixed"=0, "normal"=1, "bernoulli"=2, "poisson"=3, "gamma"=4 ),
               "y_tj" = tsdata )

  # Construct parameters
  if( is.null(control$parameters) ){
    Params = list( "beta_z" = rep(0,max(ram[,4])),
                   "lnsigma_j" = rep(0,ncol(tsdata)),
                   "mu_j" = rep(0,ncol(tsdata)),
                   "delta0_j" = rep(0,ncol(tsdata)),
                   "x_tj" = ifelse( is.na(tsdata), 0, tsdata ) )
    #if( control$gmrf_parameterization=="separable" ){
    #  Params$x_tj = ifelse( is.na(tsdata), 0, tsdata )
    #}else{
    #  Params$eps_tj = ifelse( is.na(tsdata), 0, tsdata )
    #}

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
    Params = control$parameters
  }

  # Construct map
  if( is.null(control$map) ){
    Map = list()
    # Map off x_tj for fixed when data is available
    Map$x_tj = factor(ifelse( is.na(as.vector(tsdata)) | (Data$familycode_j[col(tsdata)] %in% c(1,2,3,4)), seq_len(prod(dim(tsdata))), NA ))
    # Map off sigma_j for fixed / bernoulli / Poisson
    Map$lnsigma_j = factor( ifelse(Data$familycode_j %in% c(0,2,3), NA, seq_along(Params$lnsigma_j)) )

    # Map off mean for latent variables
    Map$mu_j = factor( ifelse(colSums(!is.na(tsdata))==0, NA, 1:ncol(tsdata)) )
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
  obj = TMB::MakeADFun( data=Data,
                   parameters=Params,
                   random=Random,
                   map=Map,
                   profile = control$profile,
                   DLL="dsem",
                   silent = TRUE )
  if(control$quiet==FALSE) list_parameters(obj)
  # bundle
  internal = list(
    sem = sem,
    tsdata = tsdata,
    family = family,
    estimate_delta0 = estimate_delta0,
    control = control,
    covs = covs,
    prior_negloglike = prior_negloglike
  )

  # Parse priors
  if( !is.null(prior_negloglike) ){
    # prior_negloglike = function(obj) -dnorm(obj$par[1],0,1,log=TRUE)
    prior_value = tryCatch( expr = prior_negloglike(obj) )
    if( is.na(prior_value) ) stop("Check `prior_negloglike(obj$par)`")
    obj$fn_orig = obj$fn
    obj$gr_orig = obj$gr

    # BUild prior evaluator
    requireNamespace("RTMB")
    priors_obj = RTMB::MakeADFun( func = prior_negloglike, 
                                  parameters = list(par=obj$par), 
                                  silent = TRUE )
    obj$fn = function(pars) obj$fn_orig(pars) + priors_obj$fn(pars)
    obj$gr = function(pars) obj$gr_orig(pars) + priors_obj$gr(pars)
    internal$priors_obj = priors_obj
  }
  
  # Further bundle
  out = list( "obj"=obj,
              "ram"=ram,
              "sem_full"=out$model,
              "tmb_inputs"=list("data"=Data, "parameters"=Params, "random"=Random, "map"=Map),
              #"call" = match.call(),
              "internal" = internal )

  # Export stuff
  if( control$run_model==FALSE ){
    class(out) = "dsem"
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
  out$run_time = Sys.time() - start_time

  # output
  class(out) = "dsem"
  return(out)
}

#' @title Detailed control for dsem structure
#'
#' @description Define a list of control parameters.  Note that
#' the format of this input is likely to change more rapidly than that of
#' \code{\link{dsem}}
#'
#' @inheritParams TMB::MakeADFun
#'
#' @param nlminb_loops Integer number of times to call \code{\link[stats]{nlminb}}.
#' @param newton_loops Integer number of Newton steps to do after running
#'        \code{\link[stats]{nlminb}}.
#' @param trace Parameter values are printed every `trace` iteration
#'        for the outer optimizer. Passed to
#'        `control` in \code{\link[stats]{nlminb}}.
#' @param eval.max Maximum number of evaluations of the objective function
#'        allowed. Passed to `control` in \code{\link[stats]{nlminb}}.
#' @param iter.max Maximum number of iterations allowed. Passed to `control` in
#'        \code{\link[stats]{nlminb}}.
#' @param getsd Boolean indicating whether to call \code{\link[TMB]{sdreport}}
#' @param run_model Boolean indicating whether to estimate parameters (the default), or
#'        instead to return the model inputs and compiled TMB object without running;
#' @param gmrf_parameterization Parameterization to use for the Gaussian Markov 
#'        random field, where the default `separable` constructs a precision matrix
#'        that must be full rank, and the alternative `projection` constructs
#'        a full-rank and IID precision for variables over time, and then projects
#'        this using the inverse-cholesky of the precision, where this projection
#'        can be rank-deficient.
#' @param constant_variance Whether to specify a constant conditional variance 
#'        \eqn{ \mathbf{\Gamma \Gamma}^t} using the default \code{constant_variance="conditional"}, 
#'        which results in a changing marginal variance      
#'        along the specified causal graph when lagged paths are present. Alternatively, the user can
#'        specify a constant marginal variance using \code{constant_variance="diagonal"}
#'        or \code{constant_variance="marginal"},
#'        such that \eqn{ \mathbf{\Gamma}} and \eqn{\mathbf{I-P}} are rescaled to achieve this constraint.  
#'        All options
#'        are equivalent when the model includes no lags (only simultaneous effects) and
#'        no covariances (no two-headed arrows).  \code{"diagonal"} and \code{"marginal"}
#'        are equivalent when the model includes no covariances. Given some exogenous covariance, 
#'        \code{constant_variance = "diagonal"} preserves the conditional correlation and has
#'        changing conditional variance, while \code{constant_variance = "marginal"} has changing
#'        conditional correlation along the causal graph.  
#' @param quiet Boolean indicating whether to run model printing messages to terminal or not;
#' @param use_REML Boolean indicating whether to treat non-variance fixed effects as random,
#'        either to motigate bias in estimated variance parameters or improve efficiency for
#'        parameter estimation given correlated fixed and random effects
#' @param parameters list of fixed and random effects, e.g., as constructed by \code{dsem} and then modified
#'        by hand (only helpful for advanced users to change starting values or restart at intended values)
#' @param map list of fixed and mirrored parameters, constructed by \code{dsem} by default but available
#'        to override this default and then pass to \code{\link[TMB]{MakeADFun}}
#' @param getJointPrecision whether to get the joint precision matrix.  Passed
#'        to \code{\link[TMB]{sdreport}}.
#' @param extra_convergence_checks Boolean indicating whether to run extra checks on model
#'        convergence.
#' @param lower vectors of lower bounds, replicated to be as long as start and passed to \code{\link[stats]{nlminb}}.
#'        If unspecified, all parameters are assumed to be unconstrained.
#' @param upper vectors of upper bounds, replicated to be as long as start and passed to \code{\link[stats]{nlminb}}.
#'        If unspecified, all parameters are assumed to be unconstrained.
#'
#' @return
#' An S3 object of class "dsem_control" that specifies detailed model settings,
#' allowing user specification while also specifying default values
#'
#' @export
dsem_control <-
function( nlminb_loops = 1,
          newton_loops = 1,
          trace = 0,
          eval.max = 1000,
          iter.max = 1000,
          getsd = TRUE,
          quiet = FALSE,
          run_model = TRUE,
          gmrf_parameterization = c("separable", "projection"),
          constant_variance = c("conditional", "marginal", "diagonal"),
          use_REML = TRUE,
          profile = NULL,
          parameters = NULL,
          map = NULL,
          getJointPrecision = FALSE,
          extra_convergence_checks = TRUE,
          lower = -Inf,
          upper = Inf ){

  gmrf_parameterization = match.arg(gmrf_parameterization)
  constant_variance = match.arg(constant_variance)

  # Return
  structure( list(
    nlminb_loops = nlminb_loops,
    newton_loops = newton_loops,
    trace = trace,
    eval.max = eval.max,
    iter.max = iter.max,
    getsd = getsd,
    quiet = quiet,
    run_model = run_model,
    gmrf_parameterization = gmrf_parameterization,
    constant_variance = constant_variance,
    use_REML = use_REML,
    profile = profile,
    parameters = parameters,
    map = map,
    getJointPrecision = getJointPrecision,
    extra_convergence_checks = extra_convergence_checks,
    lower = lower,
    upper = upper
  ), class = "dsem_control" )
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
#' @return
#' Returns a data.frame summarizing estimated path coefficients, containing columns:
#' \describe{
#' \item{path}{The parsed path coefficient}
#' \item{lag}{The lag, where e.g. 1 means the predictor in time t effects the response in time t+1}
#' \item{name}{Parameter name}
#' \item{start}{Start value if supplied, and NA otherwise}
#' \item{parameter}{Parameter number}
#' \item{first}{Variable in path treated as predictor}
#' \item{second}{Variable in path treated as response}
#' \item{direction}{Whether the path is one-headed or two-headed}
#' \item{Estimate}{Maximum likelihood estimate}
#' \item{Std_Error}{Estimated standard error from the Hessian matrix}
#' \item{z_value}{Estimate divided by Std_Error}
#' \item{p_value}{P-value associated with z_value using a two-sided Wald test}
#' }
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
  if( "sdrep" %in% names(object) ){
    SE = as.list( object$sdrep, report=FALSE, what="Std. Error")
    coefs = data.frame( coefs, "Std_Error"=c(NA,SE$beta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
    coefs = data.frame( coefs, "z_value"=coefs[,'Estimate']/coefs[,'Std_Error'] )
    coefs = data.frame( coefs, "p_value"=pnorm(-abs(coefs[,'z_value'])) * 2 )
  }

  return(coefs)
}


#' @title Simulate dsem
#'
#' @description Plot from a fitted \code{dsem} model
#'
#' @param x Output from \code{\link{dsem}}
#' @param y Not used
#' @param edge_label Whether to plot parameter names, estimated values,
#'        or estimated values along with stars indicating significance at
#'        0.05, 0.01, or 0.001 levels (based on two-sided Wald tests)
#' @param digits integer indicating the number of decimal places to be used
#' @param style Whether to make a graph using \code{igraph} or \code{ggraph}
#' @param ... arguments passed to \code{\link[igraph]{plot.igraph}}
#'
#' @details
#' This function coerces output from a graph and then plots the graph.
#'
#' @return
#' Invisibly returns the output from \code{\link[igraph]{graph_from_data_frame}}
#' which was passed to \code{\link[igraph]{plot.igraph}} for plotting.
#'
#' @method plot dsem
#' @export
plot.dsem <-
function( x,
          y,
          edge_label = c("name","value","value_and_stars"),
          digits = 2,
          style = c("igraph","ggraph"),
          ... ){

  style = match.arg(style)
  edge_label = match.arg(edge_label)

  # Extract stuff
  out = summary(x)

  # Format inputs
  from = ifelse( out[,2]==0, out$first, paste0("lag(",out$first,",",out[,2],")"))
  vertices = union( out$second, from )
  DF = data.frame(from=from, to=out$second, label=out[,3])
  if( edge_label=="value"){
    DF$label = round(out$Estimate, digits=digits)
  }
  if( edge_label == "value_and_stars" ){
    DF$label = round(out$Estimate, digits=digits)
    add_stars = cut(out[,'p_value'], breaks=c(1,0.05,0.01,0.001,0) )
    add_stars = c("***","**","*","")[as.numeric(add_stars)]
    DF$label = paste0(DF$label, add_stars)
  }

  # Create and plotgraph
  pg <- graph_from_data_frame( d = DF,
                               directed = TRUE,
                               vertices = data.frame(vertices) )

  # Two plot styles
  if(style=="igraph"){
    coords = layout_(pg, with_sugiyama())
    plot( pg, layout = coords, ... )
  }
  if(style=="ggraph"){
    # Modified from phylopath::plot.DAG
    algorithm = 'sugiyama'
    manual_layout = NULL
    text_size = 6
    box_x = 12
    box_y = 8
    edge_width = 1
    curvature = 0
    rotation = 0
    flip_x = FALSE
    flip_y = FALSE
    label = DF$label
    l = ggraph::create_layout(pg, 'igraph', algorithm = algorithm)
    arrow = grid::arrow(type = 'closed', 18, grid::unit(15, 'points'))
    gplot = ggraph::ggraph(l) +
      ggraph::geom_edge_arc(
        aes(label = label),
        strength = curvature, arrow = arrow, edge_width = edge_width,
        end_cap = ggraph::rectangle(box_x, box_y, 'mm'),
        start_cap = ggraph::rectangle(box_x, box_y, 'mm')
      ) +
      ggraph::geom_node_text(ggplot2::aes_(label = ~name), size = text_size) +
      ggraph::theme_graph(base_family = 'sans')
    plot(gplot)
  }
  return(invisible(pg))
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
#' @param fill_missing whether to fill in simulate all data (including values
#'        that are missing in the original data set)
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
#' @return
#' Simulated data, either from \code{obj$simulate} where \code{obj} is the compiled
#' TMB object, first simulating a new GMRF and then calling \code{obj$simulate}.
#'
#' @method simulate dsem
#' @export
simulate.dsem <-
function( object,
          nsim = 1,
          seed = NULL,
          variance = c("none", "random", "both"),
          resimulate_gmrf = FALSE,
          fill_missing = FALSE,
          ... ){

  # Front stuff
  set.seed(seed)
  variance = match.arg(variance)

  # Sample from GMRF using sparse precision
  rmvnorm_prec <- function(mu, prec, nsim) {
    z <- matrix(rnorm(length(mu) * nsim), ncol=nsim)
    L <- Cholesky(prec, super=TRUE)
    z <- solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    return(mu + z)
  }

  # pull out objects for easy use
  obj = object$obj
  parfull = obj$env$parList()
  tsdata = object$internal$tsdata

  # Extract parameters, and add noise as desired
  par_zr = outer( obj$env$last.par.best, rep(1,nsim) )
  if( variance=="random" ){
    eps_zr = rmvnorm_prec( rep(0,length(obj$env$random)), obj$env$spHess(random=TRUE), nsim=nsim )
    par_zr[obj$env$random,] = par_zr[obj$env$random,,drop=FALSE] + eps_zr
  }
  if( variance=="both" ){
    if(is.null(object$sdrep$jointPrecision)){
      stop("Please re-run `dsem` with `getsd=TRUE` and `getJointPrecision=TRUE`, or confirm that the model is converged")
    }
    eps_zr = rmvnorm_prec( rep(0,length(obj$env$last.par)), object$sdrep$jointPrecision, nsim=nsim )
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
      tmp = rmvnorm_prec( as.vector(newrep$delta_tj + newrep$xhat_tj), Q_kk, nsim=1 )
      # Modify call
      #newcall = object$call
      # Get control
      #newcall$control = eval(newcall$control)
      #newcall$control$parameters = newparfull
      #newcall$control$parameters$x_tj[] = tmp
      # Rebuild model with new GMRF values
      #newcall$control$run_model = FALSE
      #newfit = eval(newcall)
      control = object$internal$control
      control$parameters = newparfull
      control$parameters$x_tj[] = tmp
      control$run_model = FALSE
      newfit = dsem( sem = object$internal$sem,
                     tsdata = object$internal$tsdata,
                     family = object$internal$family,
                     estimate_delta0 = object$internal$estimate_delta0,
                     control = control )
      out[[r]] = newfit$obj$simulate()$y_tj
    }else{
      out[[r]] = obj$simulate( par_zr[,r] )$y_tj
    }
    if(isFALSE(fill_missing)){
      # Use missingness pattern from original data
      out[[r]] = ifelse( is.na(tsdata), NA, out[[r]] )
    }
    colnames(out[[r]]) = colnames(tsdata)
    tsp(out[[r]]) = attr(tsdata,"tsp")
    class(out[[r]]) = class(tsdata)
  }

  return(out)
}

#' @title Extract Variance-Covariance Matrix
#'
#' @description extract the covariance of fixed effects, or both fixed and random effects.
#'
#' @param object output from \code{dsem}
#' @param which whether to extract the covariance among fixed effects, random effects, or both
#' @param ... ignored, for method compatibility
#' @importFrom stats vcov
#' @method vcov dsem
#'
#' @return
#' A square matrix containing the estimated covariances among the parameter estimates in the model.
#' The dimensions dependend upon the argument \code{which}, to determine whether fixed, random effects,
#' or both are outputted.
#'
#' @export
vcov.dsem <-
function( object,
          which = c("fixed", "random", "both"),
          ...) {

  which = match.arg(which)

  if( which=="fixed" ){
    V = object$sdrep$cov.fixed
    if(is.null(V)){
      warning("Please re-run `dsem` with `getsd=TRUE`, or confirm that the model is converged")
    }
  }
  if( which=="random" ){
    V = solve(object$obj$env$spHess(random=TRUE))
  }
  if( which=="both" ){
    H = object$sdrep$jointPrecision
    if(is.null(H)){
      warning("Please re-run `dsem` with `getsd=TRUE` and `getJointPrecision=TRUE`, or confirm that the model is converged")
      V = NULL
    }else{
      V = solve(H)
    }
  }

  return( V )
}

#' @title Calculate residuals
#'
#' @description Calculate deviance or response residuals for dsem
#'
#' @param object Output from \code{\link{dsem}}
#' @param type which type of residuals to compute (only option is \code{"deviance"} or \code{"response"} for now)
#' @param ... Not used
#'
#' @method residuals dsem
#'
#' @return
#' A matrix of residuals, with same order and dimensions as argument \code{tsdata}
#' that was passed to \code{dsem}.
#'
#' @export
residuals.dsem <-
function( object,
          type = c("deviance","response"),
          ... ){

  # https://stats.stackexchange.com/questions/1432/what-do-the-residuals-in-a-logistic-regression-mean
  # Normal deviance residuals
  #if( FALSE ){
  #  x = rnorm(10)
  #  y = x + rnorm(10)
  #  Glm = glm( y ~ 1 + x, family="gaussian")
  #  mu = predict(Glm,type="response")
  #  r1 = y - mu
  #  r1 - resid(Glm)
  #}
  ## Poisson deviance residuals
  #if( FALSE ){
  #  x = rnorm(10)
  #  y = rpois(10, exp(x))
  #  Glm = glm( y ~ 1 + x, family="poisson")
  #  mu = predict(Glm,type="response")
  #  # https://stats.stackexchange.com/questions/398098/formula-for-deviance-residuals-for-poisson-model-with-identity-link-function
  #  r1 = sign(y - mu) * sqrt(2*(y*log((y+1e-10)/mu) - (y-mu)))
  #  r1 - resid(Glm)
  #}
  ## Binomial deviance residuals
  #if( FALSE ){
  #  p = 0.5
  #  y = rbinom(10, prob=p, size=1)
  #  Glm = glm( y ~ 1, family="binomial")
  #  mu = predict(Glm, type="response")
  #  r1 = sign(y - mu) * sqrt(-2*(((1-y)*log(1-mu) + y*log(mu))))
  #  r1 - resid(Glm)
  #}
  ## Gamma deviance residuals
  #if( FALSE ){
  #  mu = 1
  #  cv = 0.8
  #  y = rgamma( n=10, shape=1/cv^2, scale=mu*cv^2 )
  #  Glm = glm( y ~ 1, family=Gamma(link='log'))
  #  mu = predict(Glm, type="response")
  #  r1 = sign(y - mu) * sqrt(2 * ( (y-mu)/mu - log(y/mu) ))
  #  r1 - resid(Glm)
  #}

  # Poisson: sign(y - mu) * sqrt(2*(ifelse(y==0, 0, y*log(y/mu)) - (y-mu)))
  # Binomial:  -2 * ((1-y)*log(1-mu) + y*log(mu))
  # Gamma: 2 * ( (y-mu)/mu - log(y/mu) )

  # Easy of use
  #z_tj = object$obj$report()$z_tj
  y_tj = object$tmb_inputs$data$y_tj
  #familycode_j = object$tmb_inputs$data$familycode_j
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

#' @title Print fitted dsem object
#'
#' @description Prints output from fitted dsem model
#'
#' @param x Output from \code{\link{dsem}}
#' @param ... Not used
#'
#' @return
#' No return value, called to provide clean terminal output when calling fitted
#' object in terminal.
#'
#' @method print dsem
#' @export
print.dsem <-
function( x,
          ... ){
  print(x$opt)
}

#' @title predictions using dsem
#'
#' @description Predict variables given new (counterfactual) values of data, or for future or past times
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
#' @param ... Not used
#'
#' @return
#' A matrix of predicted values with dimensions and order corresponding to
#' argument \code{newdata} is provided, or \code{tsdata} if not.
#' Predictions are provided on either link or response scale, and
#' are generated by re-optimizing random effects condition on MLE
#' for fixed effects, given those new data.
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
    #newcall = object$call
    #newcall$tsdata = newdata
    # Rebuild model with new data
    #newcall$run_model = FALSE
    #newfit = eval(newcall)
    control = object$internal$control
    control$run_model = FALSE
    if( inherits(fit,"dsemRTMB") ){
      newfit = dsemRTMB( sem = object$internal$sem,
                     tsdata = newdata,
                     family = object$internal$family,
                     estimate_delta0 = object$internal$estimate_delta0,
                     control = object$internal$control,
                     covs = object$internal$covs )
    }else{
      newfit = dsem( sem = object$internal$sem,
                     tsdata = newdata,
                     family = object$internal$family,
                     estimate_delta0 = object$internal$estimate_delta0,
                     control = object$internal$control,
                     covs = object$internal$covs )
    }
    # Optimize random effects given original MLE and newdata
    newfit$obj$fn( object$opt$par )
    # Return predictor
    if(type=="link") out = newfit$obj$env$parList()$x_tj
    if(type=="response") out = newfit$obj$report()$mu_tj
  }

  return(out)
}

#' @title Marginal log-likelihood
#'
#' @description Extract the (marginal) log-likelihood of a dsem model
#'
#' @param object Output from \code{\link{dsem}}
#' @param ... Not used
#'
#' @return object of class \code{logLik} with attributes
#'   \item{val}{log-likelihood}
#'   \item{df}{number of parameters}
#' @importFrom stats logLik
#'
#' @return
#' Returns an object of class logLik. This has attributes
#' "df" (degrees of freedom) giving the number of (estimated) fixed effects
#' in the model, abd "val" (value) giving the marginal log-likelihood.
#' This class then allows \code{AIC} to work as expected.
#'
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
          what = c("Estimate","Std_Error","p_value"),
          direction = 1 ){

  what = match.arg(what)
  coefs = summary( fit )
  coefs = coefs[ which(coefs[,2]==lag), ]
  coefs = coefs[ which(coefs[,'direction'] %in% direction), ]

  #
  #vars = unique( c(coefs[,'first'],coefs[,'second']) )
  vars = colnames(fit$tmb_inputs$data$y_tj)
  out = list( "coef"=array(0, dim=rep(length(vars),2), dimnames=list(vars,vars)) )
  out$coef[as.matrix(coefs[,c('first','second')])] = coefs[,what]

  class(out) = "fitted_DAG"
  return(out)
}

#' @title Convert dsem to sem output
#'
#' @description Convert output from package dsem to sem
#'
#' @param object Output from \code{\link{dsem}}
#' @param lag what lag to extract and visualize
#'
#' @return Convert output to format supplied by \code{\link[sem]{sem}}
#'
#' @export
as_sem <-
function( object,
          lag = 0 ){

  Rho = t(as_fitted_DAG( object, what="Estimate", direction=1, lag=lag )$coef)
  Gamma = as_fitted_DAG( object, what="Estimate", direction=2, lag=lag )$coef
  Gammainv = diag(1/diag(Gamma))
  Linv = Gammainv %*% (diag(nrow(Rho))-Rho)
  Sinv = t(Linv) %*% Linv
  Sprime = solve(Sinv)
  Sprime = 0.5*Sprime + 0.5*t(Sprime)

  model = object$sem_full
  model = model[model[,2]==0,c(1,3,4)]
  out = sem( as.matrix(model),
             S = Sprime,
             N = nrow(object$internal$tsdata) )

  # pass out
  return(out)

  #x = rnorm(10)
  #y = x + rnorm(10)
  #object = dsem( sem="x->y, 0, beta", tsdata=ts(cbind(x,y)) )
  #mysem = as_sem(object)
  #myplot = semPlot::semPlotModel( mysem )
  #semPlot::semPaths( myplot,
  #                   whatLabels = "est",
  #                   edge.label.cex = 1.5,
  #                   node.width = 4,
  #                   node.height = 2,
  #                   shapeMan = "rectangle",
  #                   edge.width = 4,
  #                   nodeLabels = myplot@Vars$name,
  #                   nDigits=4 )
}


