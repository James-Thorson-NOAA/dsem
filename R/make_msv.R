

#' @title Make text for multivariate stochastic volatility model
#'
#' @description
#' Make the text string for a multivariate stochastic volatility model using
#' arrow-lag-slope notation for DSEM.
#'
#' @param variables Character string of variables (i.e., column names of \code{tsdata}).
#' @param n_factors Number of stochastic volatility factors.
#' @param factor_names Optional character-vector of factor names,
#'        which must match NA columns in \code{tsdata}.
#' @param collapse_text whether to collapse text into long character string
#'
#' @return
#' A text string to be passed to \code{\link{dsem}}
#'
#' @export
make_msv <-
function( variables,
          n_factors,
          factor_names = paste0("F",seq_len(n_factors)),
          collapse_text = TRUE ){

  # pre-processing
  n_variables = length( variables )
  text = NULL
  collapse = function( char ){
    char = paste( char, collapse=", " )
    return(char)
  }

  # Factor SDs
  for( f in 1:n_factors ){
    SD = c( paste(factor_names[f], "<->", factor_names[f]), 0, paste0("sigma",factor_names[f]) )
    SD = collapse( SD )
    text = c( text, SD )
  }

  # Factor RWs
  for( f in 1:n_factors ){
    AR = c( paste(factor_names[f], "->", factor_names[f]), 1, paste0("rho",factor_names[f]) )
    AR = collapse( AR )
    text = c( text, AR )
  }

  # Factor loadings
  for( f in 1:n_factors ){
  for( v in f:n_variables ){
    SV = c( paste(variables[v], "<->", variables[v]), 0, factor_names[f] )
    SV = collapse( SV )
    text = c( text, SV )
  }}

  # ARs for variables
  for( v in 1:n_variables ){
    AR = c( paste(variables[v], "->", variables[v]), 1, "rho" )
    AR = collapse( AR )
    text = c( text, AR )
  }

  # Output
  if(isTRUE(collapse_text)) text = paste0( text, collapse="\n" )
  return( text )
}
