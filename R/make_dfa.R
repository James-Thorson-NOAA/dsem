
#' @title Make text for dynamic factor analysis
#'
#' @description
#' Make the text string for a dynamic factor analysis expressed using
#' arrow-and-lag notation for DSEM.
#'
#' @param variables Character string of variables (i.e., column names of \code{tsdata}).
#' @param n_factors Number of factors.
#' @param factor_names Optional character-vector of factor names,
#'        which must match NA columns in \code{tsdata}.
#'
#' @return
#' A text string to be passed to \code{\link{dsem}}
#'
#' @export
make_dfa <-
function( variables,
          n_factors,
          factor_names = paste0("F",seq_len(n_factors)) ){

  # pre-processing
  n_variables = length( variables )
  text_matrix = NULL

  # Factor SDs
  for( f in 1:n_factors ){
    SD = c( paste(factor_names[f], "<->", factor_names[f]), 0, NA, 1 )
    text_matrix = rbind( text_matrix, SD )
  }

  # Factor RWs
  for( f in 1:n_factors ){
    AR = c( paste(factor_names[f], "->", factor_names[f]), 1, NA, 1 )
    text_matrix = rbind( text_matrix, AR )
  }

  # Factor loadings
  for( f in 1:n_factors ){
  for( v in f:n_variables ){
    Load = c( paste(factor_names[f], "->", variables[v]), 0, paste0("L",f,v), 0.1 )
    text_matrix = rbind( text_matrix, Load )
  }}

  # Fix SD = 0 for additional process errors
  for( v in 1:n_variables ){
    extraSD = c( paste(variables[v], "<->", variables[v]), 0, NA, 0 )
    text_matrix = rbind( text_matrix, extraSD )
  }

  # Output
  text_vec = apply( text_matrix, paste, MARGIN=1, collapse=", " )
  text = paste0( text_vec, collapse="\n" )
  return( text )
}
