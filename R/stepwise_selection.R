
#' @title Simulate dsem
#'
#' @description Plot from a fitted \code{dsem} model
#'
#' @param model_options character-vector containing sem elements
#'        that could be included or dropped depending upon their
#'        parsimony
#' @param model_shared character-vector containing sem elements
#'        that must be included regardless of parsimony
#' @param options_initial character-vector containing some (possible empty)
#'        subset of \code{model_options}, where stepwise selection begins
#'        with that set of model options included.
#' @param quiet whether to avoid displaying progress to terminal
#' @param criterion function that computes the information criterion to be
#'        minimized, typically using \code{AIC}.  However, users can instead supply
#'        a function that computes CIC using
#'        \code{\link{test_dsep}} and desired settings, presumably including a \code{set.seed}
#'        if missing data are being imputed
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
#' @examples
#' # Simulate x -> y -> z
#' set.seed(101)
#' x = rnorm(100)
#' y = 0.5*x + rnorm(100)
#' z = 1*y + rnorm(100)
#' tsdata = ts(data.frame(x=x, y=y, z=z))
#'
#' # define candidates
#' model_options = c(
#'   "y -> z, 0, y_to_z",
#'   "x -> z, 0, x_to_z"
#' )
#' # define paths that are required
#' model_shared = "
#'   x -> y, 0, x_to_y
#' "
#'
#' # Do selection
#' step = stepwise_selection(
#'   model_options = model_options,
#'   model_shared = model_shared,
#'   tsdata = tsdata,
#'   quiet = TRUE
#' )
#'
#' # Check selected model
#' cat(step$model)
#'
#' @export
stepwise_selection <-
function( model_options,
          model_shared,
          options_initial = c(),
          quiet = FALSE,
          criterion = AIC,
          ... ){

  # Loop
  best = (model_options %in% options_initial)
  step = NULL
  while(TRUE){
    if(isFALSE(quiet)) message("Running with ", sum(best), " vars included: ", Sys.time() )
    df_options = outer( rep(1,length(best)), best )
    which_diag = cbind( seq_len(nrow(df_options)), seq_len(nrow(df_options)) )
    df_options[which_diag] = 1 - df_options[which_diag]
    df_options = rbind(best, df_options)
    IC_i = rep(NA, nrow(df_options))

    for(i in 1:nrow(df_options)){
      model = paste( paste(model_options[which(df_options[i,]==1)],collapse="\n"),
                     paste(model_shared,collapse="\n"),
                     sep="\n" )
      myfit = dsem( sem = model, ... )
      IC_i[i] = criterion(myfit)
    }

    # Save step and decide whether to continue
    step[[length(step)+1]] = cbind( IC_i, df_options )
    if(which.min(IC_i)==1){
      break()
    }else{
      best = df_options[which.min(IC_i),]
    }
  }

  # Return best
  model = paste( paste(model_options[which(best==1)],collapse="\n"),
                 paste(model_shared,collapse="\n"),
                 sep="\n" )
  return(list(model=model, step=step))
}
