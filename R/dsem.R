#' Fit dynamic structural equation model
#'
#' Fits a dynamic structural equation model
#'
#' @param sem structural equation model structure, passed to either \code{\link[sem]{specifyModel}}
#'        or \code{\link[sem]{specifyEquations}} and then parsed to control
#'        the set of path coefficients and variance-covariance parameters
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
#'   Consumption <-> Consumption, 0, v1
#'   Consumption -> Consumption, -1, ar1
#'   Consumption -> Consumption, -2, ar2
#'   Profits -> Investment, 0, b2
#'   Profits -> Investment, -1, b3
#'   Capital_stock -> Investment, -1, b4
#'   Investment <-> Investment, 0, v2
#'   neg_Gov_wage <-> neg_Gov_wage, 0, v3
#'   GNP -> Priv_wage, 0, c2
#'   Taxes -> Priv_wage, 0, c2
#'   neg_Gov_wage -> Priv_wage, 0, c2
#'   GNP -> Priv_wage, -1, c3
#'   Taxes -> Priv_wage, -1, c3
#'   neg_Gov_wage -> Priv_wage, -1, c3
#'   Time -> Priv_wage, 0, c4
#'   Priv_wage <-> Priv_wage, 0, v4
#'   GNP <-> GNP, 0, v5
#'   Profits <-> Profits, 0, v6
#'   Capital_stock <-> Capital_stock, 0, v7
#'   Taxes <-> Taxes, 0, v8
#'   Time <-> Time, 0, v9
#'   Gov_wage <-> Gov_wage, 0, v10
#'   Gov_expense <-> Gov_expense, 0, v11
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
          estimate_delta0 = FALSE,
          quiet = FALSE,
          run_model = TRUE,
          use_REML = TRUE,
          parameters = NULL,
          map = NULL,
          ... ){

  # (I-Rho)^-1 * Gamma * (I-Rho)^-1
  out = make_ram( sem, tsdata=tsdata, quiet=quiet )
  ram = out$ram

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
              "model"=out$model,
              "tmb_inputs"=list("data"=Data, "parameters"=Params, "random"=Random, "map"=Map) )

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

#' summarize dsem
#'
#' @title Summarize dsem
#'
#' @param object Output from \code{\link{dsem}}
#' @param ... Note used
#'
#' @method summary dsem
#' @export
summary.dsem <-
function( object, ... ){

  # Easy of use
  model = object$model
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

#' Convert dsem to phylopath output
#'
#' @title Convert output from package dsem to phylopath
#'
#' @param fit Output from \code{\link{dsem}}
#' @param lag which lag to output
#' @param what whether to output estimates \code{what="Estimate"} or standard errors \code{what="Std_Error"}
#'
#' @return Convert output to format supplied by \code{\link[phylopath]{est_DAG}}
#'
#' @export
as_fitted_DAG <-
function( fit,
          lag = 0,
          what = "Estimate" ){

  coefs = summary( fit )
  coefs = coefs[ which(coefs[,2]==lag), ]
  coefs = coefs[ which(coefs[,'direction']==1), ]

  #
  vars = unique( c(coefs[,'first'],coefs[,'second']) )
  out = list( "coef"=array(0, dim=rep(length(vars),2), dimnames=list(vars,vars)) )
  out$coef[as.matrix(coefs[,c('first','second')])] = coefs[,what]

  class(out) = "fitted_DAG"
  return(out)
}
