
#' @title Simulate dsem
#'
#' @description Plot from a fitted \code{dsem} model
#'
#' @param model_options character-vector containing sem elements
#'        that could be included or dropped depending upon their
#'        parsimony
#' @param model_shared character-vector containing sem elements
#'        that must be included regardless of parsimony
#' @param quiet whether to avoid displaying progress to terminal
#' @param ... arguments passed to \code{\link{dsem}},
#'        other than \code{sem} e.g., \code{tsdata}, \code{family}
#'        etc.
#'
#' @details
#' This function conducts stepwise (i.e., forwards and backwards) model
#' selection using marginal AIC, while forcing some model elements to be
#' included and selecting among others.
#'
#' @return
#' An object (list) that includes:
#' \describe{
#' \item{model}{the string with the selected SEM model}
#' \item{record}{a list showing the AIC and whether each \code{model_options} is included or not}
#' }
#'
#' @export
stepwise_selection <-
function( model_options,
          model_shared,
          quiet = FALSE,
          ... ){

  # Loop
  best = rep(0, length(model_options) )
  record = NULL
  while(TRUE){
    if(isFALSE(quiet)) message("Running with ", sum(best), " vars included: ", Sys.time() )
    df_options = outer( rep(1,length(best)), best )
    which_diag = cbind( seq_len(nrow(df_options)), seq_len(nrow(df_options)) )
    df_options[which_diag] = 1 - df_options[which_diag]
    df_options = rbind(best, df_options)
    AIC_i = rep(NA, nrow(df_options))

    for(i in 1:nrow(df_options)){
      model = paste( paste(model_options[which(df_options[i,]==1)],collapse="\n"),
                     paste(model_shared,collapse="\n"),
                     sep="\n" )
      myfit = dsem( sem = model, ... )
      AIC_i[i] = AIC(myfit)
    }

    # Save record and decide whether to continue
    record[[length(record)+1]] = cbind( AIC_i, df_options )
    if(which.min(AIC_i)==1){
      break()
    }else{
      best = df_options[which.min(AIC_i),]
    }
  }

  # Return best
  model = paste( paste(model_options[which(best==1)],collapse="\n"),
                 paste(model_shared,collapse="\n"),
                 sep="\n" )
  return(list(model=model, record=record))
}
