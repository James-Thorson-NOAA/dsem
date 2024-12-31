
#' @title Make a RAM (Reticular Action Model)
#'
#' @description \code{read_model} converts SEM arrow notation to \code{model} describing SEM parameters
#'
#' @inheritParams dsem
#' @param times A character vector listing the set of times in order
#' @param variables A character vector listing the set of variables
#' @param covs A character vector listing variables for which to estimate a standard deviation
#' @param quiet Boolean indicating whether to print messages to terminal
#'
#' @details
#' \strong{RAM specification using arrow-and-lag notation}
#'
#'
#' @export
read_model <-
function( sem,
          times,
          variables,
          covs = NULL,
          quiet = FALSE ){

  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  ####### Error checks
  if( !is.numeric(times) ) stop("`times` must be numeric in `make_dsem_ram`")

  ####### Define local functions
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
      model.2 <- cbind(
        'path' = c(model[, 1], paths),
        'lag' = c(model[,2], rep(0,nvars)),
        'name' = c(model[, 3], par.names),
        'start' = c(model[, 4], rep(NA, length(paths))) )
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
  model <- data.frame( "path"=model$path, "lag"=model$lag,
                       "name"=model$par, "start"=model$start)

  # Adding a SD automatically
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

  # Add parameter column
  par.names = model[, 3]
  pars = na.omit(unique(par.names))
  par.nos = apply(outer(pars, par.names, "=="), 2, which)
  par.nos = unlist(sapply( par.nos, FUN=\(x) ifelse(length(x)==0, 0, x) ))
  model = cbind( model, "parameter"=par.nos )

  # Add incidence to model
  model = cbind( model, first=NA, second=NA, direction=NA )
  for( i in seq_len(nrow(model)) ){
    path = parse_path(model[i,1])
    model[i,c('first','second','direction')] = unlist( path[c('first','second','direction')] )
  }

  ####### Step 2 -- Make RAM

  # Deal with fixed values
  #if( max(model$parameter) != length(beta_p) ) stop("Check beta_p")

  return(model)
}

