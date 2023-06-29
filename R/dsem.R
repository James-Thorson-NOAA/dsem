#' Fit dynamic structural equation model
#'
#' Fits a dynamic structural equation model
#'
#' @inheritParams sem::specifyModel
#' @inheritParams fit_tmb
#'
#' @param sem structural equation model structure, passed to either \code{\link[sem]{specifyModel}}
#'        or \code{\link[sem]{specifyEquations}} and then parsed to control
#'        the set of path coefficients and variance-covariance parameters
#' @param tsdata phylogenetic structure, using class \code{\link[ape]{as.phylo}}
#' @param family Character-vector listing the distribution used for each column of \code{tsdata}, where
#'        each element must be \code{fixed} or \code{normal}.
#'        \code{family="fixed"} is default behavior and assumes that a given variable is measured exactly.
#'        Other options correspond to different specifications of measurement error.
#' @param run_model Boolean indicating whether to estimate parameters (the default), or
#'        instead to return the model inputs and compiled TMB object without running;
#' @param ... Additional parameters passed to \code{\link{fit_tmb}}
#'
#' @examples
#' \dontrun{
#' # Load data set
#' library(phylopath)
#'
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
#' data("KleinI", package = "AER")
#' Data = as.data.frame(KleinI)
#' Data = cbind( Data, "time" = seq(1,22)-11 )
#' colnames(Data) = sapply( colnames(Data), FUN=switch, "consumption"="Consumption", "invest"="Investment", "cprofits"="Profits", "capital"="Capital_stock", "gwage"="Gov_wage", "pwage"="Priv_wage", "gexpenditure"="Gov_expense", "taxes"="Taxes", "time"="Time", "gnp"="GNP")
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
#' }
#'
#' @useDynLib dsem
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
  if( FALSE ){
    cpp_dir = R'(C:\Users\James.Thorson\Desktop\Git\dsem\src)'
    setwd(cpp_dir)
    dyn.unload(dynlib("dsem"))
    compile("dsem.cpp")
    dyn.load(dynlib("dsem"))
  }

  #
  Data = list( "RAM" = as.matrix(na.omit(ram[,1:4])),
               "RAMstart" = as.numeric(ram[,5]),
               "familycode_j" = sapply(family, FUN=switch, "fixed"=0, "normal"=1 ),
               "y_tj" = tsdata )

  # Construct parameters
  if( is.null(parameters) ){
    Params = list( "beta_z" = rnorm(max(ram[,4])),
                   "lnsigma_j" = rep(0,ncol(tsdata)),
                   "mu_j" = rep(0,ncol(tsdata)),
                   "delta0_j" = rep(0,ncol(tsdata)),
                   "x_tj" = ifelse( is.na(tsdata), 0, tsdata ))

    # Turn off initial conditions
    if( estimate_delta0==FALSE ){
      Params$delta0_j = numeric(0)
    }

    # Scale starting values with higher value for two-headed than one-headed arrows
    beta_type = tapply( ram[,1], INDEX=ram[,4], max)
    Params$beta_z = ifelse(beta_type==1, 0, 1)
  }else{
    Params = parameters
  }

  # Construct map
  if( is.null(map) ){
    Map = list()
    Map$x_tj = factor(ifelse( is.na(as.vector(tsdata)) | (Data$familycode_j[col(tsdata)] %in% c(1,2,3,4)), seq_len(prod(dim(tsdata))), NA ))
    Map$lnsigma_j = factor( ifelse(Data$familycode_j==0, NA, seq_along(Params$lnsigma_j)) )
  }else{
    Map = map
  }

  # Initial run
  # delta0_j being random leads to weird behavior for wolf-moose example
  if(isTRUE(use_REML)){
    Random = c( "x_tj", "mu_j" )
  }else{
    Random = "x_tj"
  }
  obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL="dsem" )
  if(quiet==FALSE) list_parameters(obj)
  out = list( "obj"=obj, "ram"=ram, "model"=out$model,
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
#' @param x Output from \code{\link{dsem}}
#' @method summary dsem
#' @export
summary.dsem <-
function( x ){

  # Easy of use
  model = x$model
  ParHat = x$obj$env$parList()

  #
  coefs = data.frame( model, "Estimate"=c(NA,ParHat$beta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
  coefs$Estimate = ifelse( is.na(coefs$Estimate), as.numeric(model[,4]), coefs$Estimate )
  if( "SD" %in% names(x$opt) ){
    SE = as.list( x$opt$SD, report=FALSE, what="Std. Error")
    coefs = data.frame( coefs, "Std_Error"=c(NA,SE$beta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
    coefs = data.frame( coefs, "z_value"=coefs[,'Estimate']/coefs[,'Std_Error'] )
    coefs = data.frame( coefs, "p_value"=pnorm(-abs(coefs[,'z_value'])) * 2 )
  }

  return(coefs)
}
