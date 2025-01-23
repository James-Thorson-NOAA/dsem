
#' @title Calculate conditional AIC
#'
#' @description
#' Calculates the conditional Akaike Information criterion (cAIC).
#'
#' @param object Output from \code{\link{dsem}}
#' @param what Whether to return the cAIC or the effective degrees of freedom
#'        (EDF) for each group of random effects.
#'
#' @details
#' cAIC is designed to optimize the expected out-of-sample predictive
#' performance for new data that share the same random effects as the
#' in-sample (fitted) data, e.g., spatial interpolation.  In this sense,
#' it should be a fast approximation to optimizing the model structure
#' based on k-fold crossvalidation.
#' By contrast, \code{AIC} calculates the
#' marginal Akaike Information Criterion, which is designed to optimize
#' expected predictive performance for new data that have new random effects,
#' e.g., extrapolation, or inference about generative parameters.
#'
#' cAIC also calculates as a byproduct the effective degrees of freedom,
#' i.e., the number of fixed effects that would have an equivalent impact on
#' model flexibility as a given random effect.
#'
#' Both cAIC and EDF are calculated using Eq. 6 of Zheng Cadigan Thorson 2024.
#'
#' Note that, for models that include profiled fixed effects, these profiles
#' are turned off.
#'
#' @return
#' Either the cAIC, or the effective degrees of freedom (EDF) by group
#' of random effects
#'
#' @references
#'
#' **Deriving the general approximation to cAIC used here**
#'
#' Zheng, N., Cadigan, N., & Thorson, J. T. (2024).
#' A note on numerical evaluation of conditional Akaike information for
#' nonlinear mixed-effects models (arXiv:2411.14185). arXiv.
#' \doi{10.48550/arXiv.2411.14185}
#'
#' **The utility of EDF to diagnose hierarchical model behavior**
#'
#' Thorson, J. T. (2024). Measuring complexity for hierarchical
#' models using effective degrees of freedom. Ecology,
#' 105(7), e4327 \doi{10.1002/ecy.4327}
#'
#' @export
cAIC <-
function( object,
          what = c("cAIC","EDF") ){

  what = match.arg(what)
  data = object$tmb_inputs$data

  # Error checks
  if(any(is.na(object$tmb_inputs$map$x_tj))){
    stop("cAIC is not implemented when fixing states at data using family=`fixed`")
  }

  # Turn on all GMRF parameters
  map = object$tmb_inputs$map
  map$x_tj = factor(seq_len(prod(dim(data$y_tj))))

  # Make sure profile = NULL
  #if( is.null(object$internal$control$profile) ){
    obj = object$obj
  #}else{
    obj = TMB::MakeADFun( data = data,
                     parameters = object$internal$parhat,
                     random = object$tmb_inputs$random,
                     map = map,
                     profile = NULL,
                     DLL="dsem",
                     silent = TRUE )
  #}

  # Weights = 0 is equivalent to data = NA
  data$y_tj[] = NA
  # Make obj_new
  obj_new = TMB::MakeADFun( data = data,
                      parameters = object$internal$parhat,
                      map = map,
                      random = object$tmb_inputs$random,
                      DLL = "dsem",
                      profile = NULL )

  #
  par = obj$env$parList()
  parDataMode <- obj$env$last.par
  indx = obj$env$lrandom()
  q = sum(indx)
  p = length(object$opt$par)

  ## use - for Hess because model returns negative loglikelihood;
  Hess_new = -Matrix::Matrix(obj_new$env$f(parDataMode,order=1,type="ADGrad"),sparse = TRUE)
  #cov_Psi_inv = -Hess_new[indx,indx]; ## this is the marginal prec mat of REs;
  Hess_new = Hess_new[indx,indx]

  ## Joint hessian etc
  Hess = -Matrix::Matrix(obj$env$f(parDataMode,order=1,type="ADGrad"),sparse = TRUE)
  Hess = Hess[indx,indx]
  #negEDF = diag(as.matrix(solve(ddlj.r)) %*% ddlr.r)
  negEDF = Matrix::diag(Matrix::solve(Hess, Hess_new))
  #

  if(what=="cAIC"){
    jnll = obj$env$f(parDataMode)
    cnll = jnll - obj_new$env$f(parDataMode)
    cAIC = 2*cnll + 2*(p+q) - 2*sum(negEDF)
    return(cAIC)
  }
  if(what=="EDF"){
    #Sdims = object$tmb_inputs$tmb_data$Sdims
    #group = rep.int( seq_along(Sdims), times=Sdims )
    #names(negEDF) = names(obj$env$last.par)[indx]
    EDF = length(negEDF) - sum(negEDF)
    return(EDF)
  }
}
