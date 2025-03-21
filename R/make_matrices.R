#' @title Make path matrices
#'
#' @description Constructs path matrices for dynamic structural equation model (DSEM)
#'        using a vector of parameters and specification of the DSEM
#'
#' @param beta_p vector parameters.
#' @param model matrix or data.frame with the following columns, and one row
#'   per one-headed or two-headed arrow in the dynamic structural model:
#' \describe{
#' \item{direction}{whether a path coefficient is one-headed (1) or two-headed (2)}
#' \item{lag}{whether the lag associated with a given coefficient}
#' \item{start}{starting value, used when \code{parameter=0}}
#' \item{parameter}{The parameter number from \code{beta_p} associated with a given path}
#' \item{first}{The variable at the tail of a given path}
#' \item{second}{The variable at the head of a given path}
#' }
#' @param times integer-vector of times to use when defining matrices
#' @param variables character-vector listing variables
#'
#' @importFrom Matrix solve Diagonal sparseMatrix drop0 kronecker
#' @importFrom RTMB ADoverload AD
#'
#' @details
#' When \code{length(times)} is \eqn{T} and \code{length(variables)} is \eqn{J},
#' \code{make_matrices} returns matrices of dimension \eqn{TJ \times TJ} representing
#' paths among \eqn{vec(\mathbf{X})} where matrix \eqn{\mathbf{X}} has dimension
#' \eqn{T \times J} and \eqn{vec} stacks columns into a single long vector
#'
#' @return
#' A named list of matrices including:
#' \describe{
#' \item{P_kk}{The matrix of interactions, i.e., one-headed arrows}
#' \item{G_kk}{The matrix of exogenous covariance, i.e., two-headed arrows}
#' }
#' @export
make_matrices <-
function( beta_p,
          model,
          times,
          variables ){

  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  model <- as.data.frame(model)
  model$parameter = as.integer(model[,'parameter'])

  # Combine fixed, estimated, and mapped parameters into vector
  beta_i = rep(0, nrow(model))
  off = which(model[,'parameter'] == 0)
  if( length(off) > 0 ){
    beta_i[off] = as.numeric(model[off,'start'])
  }
  not_off = which(model[,'parameter'] > 0)
  if( length(not_off) > 0 ){
    beta_i[not_off] = beta_p[model[not_off,'parameter']]
  }

  # Loop through paths
  P_kk = drop0(sparseMatrix( i=1, j=1, x=0, dims=rep(length(variables)*length(times),2) ))   # Make with a zero
  #P_kk = AD(P_kk)
  G_kk = (P_kk)
  for( i in seq_len(nrow(model)) ){
    lag = as.numeric(model[i,'lag'])
    # Time-lag matrix ... transpose if negative lag
    L_tt = sparseMatrix( i = seq(abs(lag)+1,length(times)),
                         j = seq(1,length(times)-abs(lag)),
                         x = 1,
                         dims = rep(length(times),2) )

    if(lag<0) L_tt = t(L_tt)
    # Interaction matrix
    P_jj = sparseMatrix( i = match(model[i,'second'],variables),
                         j = match(model[i,'first'],variables),
                         x = 1,
                         dims = rep(length(variables),2) )

    # Assemble
    tmp_kk = kronecker(P_jj, L_tt)
    if(abs(as.numeric(model[i,'direction']))==1){
      P_kk = P_kk + beta_i[i] * tmp_kk # AD(tmp_kk)
    }else{
      G_kk = G_kk + beta_i[i] * tmp_kk # AD(tmp_kk)
    }
  }

  # Diagonal component
  I_kk = Diagonal(nrow(P_kk))

  # Assemble
  IminusP_kk = AD(I_kk - P_kk)
  invV_kk = AD(G_kk)
  invV_kk@x = AD(1 / G_kk@x^2)
  #Q_kk = t(IminusP_kk) %*% invV_kk %*% IminusP_kk

  out = list(
    "P_kk" = P_kk,
    "G_kk" = G_kk,
    "invV_kk" = invV_kk,
    #"Q_kk" = Q_kk,     # NOT USED
    "IminusP_kk" = IminusP_kk
  )
  return(out)
}
