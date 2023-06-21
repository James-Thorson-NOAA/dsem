#' Convert to fitted_DAG
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
