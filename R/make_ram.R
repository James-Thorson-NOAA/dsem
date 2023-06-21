#' Make RAM
#' @export
make_ram <-
function( sem,
          tsdata,
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
      all.vars <- sem:::classifyVariables(model)
      exo.vars <- all.vars$exogenous
      end.vars <- all.vars$endogenous
      variables <- logical(0)
      for (i in seq_len(nrow(model))) {
          paths = model[i,1]
          lag = model[i,2]
          vars <- sem:::strip.white(paths)
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
  exog.variances = endog.variances = TRUE
  model = add.variances()

  ####### Step 2 -- Make RAM

  # Global stuff
  Q_dimnames = dimnames(.preformat.ts(tsdata))
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
    path = sem:::parse.path(model[i,1])
    model[i,c('first','second','direction')] = unlist( path[c('first','second','direction')] )
  }

  # Loop through paths
  for( i in seq_len(nrow(model)) ){
  for( t in 1:nrow(Z) ){
    lag = as.numeric(model[i,2])
    par.no = par.nos[i]
    # Get index for "from"
    from = c( Q_dimnames[[1]][t], model[i,'first'] )
    from_index = match_row( Q_names, from )
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
