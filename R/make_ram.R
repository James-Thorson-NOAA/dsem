#' @title Make a RAM (Reticular Action Model)
#'
#' @description \code{make_ram} converts SEM arrow notation to \code{ram} describing SEM parameters
#'
#' @inheritParams dsem
#' @param sem See Details section in \code{\link[dsem]{make_ram}}
#' @param times A character vector listing the set of times in order
#' @param variables A character vector listing the set of variables
#' @param quiet Boolean indicating whether to print messages to terminal
#' @param remove_na Boolean indicating whether to remove NA values from RAM (default) or not.
#'            \code{remove_NA=FALSE} might be useful for exploration and diagnostics for
#'            advanced users
#'
#' @details
#' Each line of the RAM specification for \code{\link[dsem]{make_ram}} consists of four (unquoted) entries,
#' separated by commas:
#'
#' \describe{
#'   \item{1. Arrow specification:}{This is a simple formula, of the form
#'     \code{A -> B} or, equivalently, \code{B <- A} for a regression
#'     coefficient (i.e., a single-headed or directional arrow);
#'     \code{A <-> A} for a variance or \code{A <-> B} for a covariance
#'     (i.e., a double-headed or bidirectional arrow). Here, \code{A} and
#'     \code{B} are variable names in the model. If a name does not correspond
#'     to an observed variable, then it is assumed to be a latent variable.
#'     Spaces can appear freely in an arrow specification, and
#'     there can be any number of hyphens in the arrows, including zero: Thus,
#'     e.g., \code{A->B}, \code{A --> B}, and \code{A>B} are all legitimate
#'     and equivalent.}
#'   \item{2. Lag (using negative values):}{An integer specifying whether the linkage
#'     is simultaneous (\code{lag=0}) or lagged (e.g., \code{lag=-1}), where
#'     only one-headed arrows can be lagged.}
#'   \item{3. Parameter name:}{The name of the regression coefficient, variance,
#'     or covariance specified by the arrow. Assigning the same name to two or
#'     more arrows results in an equality constraint. Specifying the parameter name
#'     as \code{NA} produces a fixed parameter.}
#'   \item{4. Value:}{start value for a free parameter or value of a fixed parameter.
#'     If given as \code{NA} (or simply omitted), \code{sem} will compute the start value.}
#' }
#'
#'   Lines may end in a comment following \code{#}.#'
#' Lines may end in a comment following #.
#'
#' Copied from package `sem` under licence GPL (>= 2) with permission from John Fox
#'
#' @return A reticular action module (RAM) describing dependencies
#'
#' @export
make_ram <-
function( sem,
          times,
          variables,
          covs = NULL,
          quiet = FALSE,
          remove_na = TRUE ){

  ####### Define location functions
  # helper function
  match_row = function( df, x ) which( df[1]==x[1] & df[2]==x[2] )
  #
  add.variances <- function() {
      variables <- need.variance()
      nvars <- length(variables)
      if (nvars == 0)
          return(model)
      message("NOTE: adding ", nvars, " variances to the model")
      paths <- character(nvars)
      par.names <- character(nvars)
      for (i in 1:nvars) {
          paths[i] <- paste(variables[i], "<->", variables[i])
          par.names[i] <- paste("V[", variables[i], "]", sep = "")
      }
      model.2 <- cbind(c(model[, 1], paths), c(model[,2], rep(0,nvars)), c(model[, 3],
          par.names), c(model[, 4], rep(NA, length(paths))))
      model.2
  }
  need.variance <- function() {
      all.vars <- classify_variables(model)
      exo.vars <- all.vars$exogenous
      end.vars <- all.vars$endogenous
      variables <- logical(0)
      for (i in seq_len(nrow(model))) {
          paths = model[i,1]
          lag = model[i,2]
          vars <- gsub(pattern=" ", replacement="", x=paths)
          vars <- sub("-*>", "->", sub("<-*", "<-", vars))
          vars <- sub("<->|<-", "->", vars)
          vars <- strsplit(vars, "->")[[1]]
          if ((vars[1] != vars[2]) | (lag != 0)) {
              for (a.variable in vars) {
                if (is.na(variables[a.variable]))
                  variables[a.variable] <- TRUE
              }
          }
          else {
              variables[vars[1]] <- FALSE
          }
      }
      if (!exog.variances && length(exo.vars) > 0)
          variables[exo.vars] <- FALSE
      if (!endog.variances && length(end.vars) > 0)
          variables[end.vars] <- FALSE
      names(variables)[variables]
  }

  ####### Step 2 -- Make RAM
  # convert to data frame
  model = scan( text = sem,
                what = list(path = "", lag = 1, par = "", start = 1, dump = ""),
                sep = ",",
                strip.white = TRUE,
                comment.char = "#",
                fill = TRUE,
                quiet = quiet)
  model$path <- gsub("\\t", " ", model$path)
  model$par[model$par == ""] <- NA
  model <- cbind( model$path, model$lag, model$par, model$start)

  if( !is.null(covs) ){
    for (cov in covs) {
      vars <- strsplit(cov, "[ ,]+")[[1]]
      nvar <- length(vars)
      for (i in 1:nvar) {
      for (j in i:nvar) {
        p1 = paste(vars[i], "<->", vars[j])
        p2 = if (i==j) paste("V[", vars[i], "]", sep = "") else paste("C[",vars[i], ",", vars[j], "]", sep = "")
        p3 = NA
        row <- c(p1, 0, p2, p3)
        if( any((row[1]==model[,1]) & (row[2]==model[,2])) ){
          next
        }else{
          model <- rbind(model, row, deparse.level = 0)
        }
      }}
    }
  }

  exog.variances = endog.variances = TRUE
  model = add.variances()

  ####### Step 2 -- Make RAM

  # Global stuff
  #Q_dimnames = dimnames(.preformat.ts(tsdata))
  Q_dimnames = list( times, variables )
  if(any(sapply(Q_dimnames,is.null))) stop("Check dimnames")
  Q_names = expand.grid(Q_dimnames)
  ram = NULL  # heads, to, from, parameter
  vars = Q_dimnames[[2]]

  # Deal with fixed values
  par.names = model[, 3]
  pars = na.omit(unique(par.names))
  par.nos = apply(outer(pars, par.names, "=="), 2, which)
  #par.nos = ifelse( sapply(par.nos,length)==0, 0, unlist(par.nos) )
  par.nos = unlist(sapply( par.nos, FUN=\(x) ifelse(length(x)==0, 0, x) ))
  model = cbind( model, "parameter"=par.nos )
  startvalues = model[,4]

  # Add incidence to model
  model = cbind( model, first=NA, second=NA, direction=NA )
  for( i in seq_len(nrow(model)) ){
    path = parse_path(model[i,1])
    model[i,c('first','second','direction')] = unlist( path[c('first','second','direction')] )
  }

  # Loop through paths
  for( i in seq_len(nrow(model)) ){
  for( t in seq_along(times) ){
    lag = as.numeric(model[i,2])
    par.no = par.nos[i]
    # Get index for "from"
    from = c( Q_dimnames[[1]][t], model[i,'first'] )
    from_index = match_row( Q_names, from )
    from_index = ifelse( length(from_index)==0, NA, from_index )
    # Get index for "to"
    to = c( Q_dimnames[[1]][t-lag], model[i,'second'] )
    to_index = match_row( Q_names, to )
    to_index = ifelse( length(to_index)==0, NA, to_index )
    ram_new = data.frame( "heads"=abs(as.numeric(model[i,'direction'])), "to"=to_index, "from"=from_index, "parameter"=par.no, "start"=startvalues[i] )
    ram = rbind( ram, ram_new )
  }}
  rownames(ram) = NULL

  #
  if( isTRUE(remove_na) ){
    which_keep = which(apply( ram[,1:4], MARGIN=1, FUN=\(x)!any(is.na(x)) ))
    ram = ram[ which_keep, ]
  }

  #
  out = list( "model"=model, "ram"=ram)
  return(out)
}
