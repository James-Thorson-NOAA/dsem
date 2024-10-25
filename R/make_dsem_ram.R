#' @title Make a RAM (Reticular Action Model)
#'
#' @description \code{make_dsem_ram} converts SEM arrow notation to \code{ram} describing SEM parameters
#'
#' @inheritParams dsem
#' @param times A character vector listing the set of times in order
#' @param variables A character vector listing the set of variables
#' @param covs A character vector listing variables for which to estimate a standard deviation
#' @param quiet Boolean indicating whether to print messages to terminal
#' @param remove_na Boolean indicating whether to remove NA values from RAM (default) or not.
#'            \code{remove_NA=FALSE} might be useful for exploration and diagnostics for
#'            advanced users
#'
#' @details
#' \strong{RAM specification using arrow-and-lag notation}
#'
#' Each line of the RAM specification for \code{\link[dsem]{make_dsem_ram}} consists of four (unquoted) entries,
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
#'   \item{2. Lag (using positive values):}{An integer specifying whether the linkage
#'     is simultaneous (\code{lag=0}) or lagged (e.g., \code{X -> Y, 1, XtoY}
#'     indicates that X in time T affects Y in time T+1), where
#'     only one-headed arrows can be lagged. Using positive values to indicate lags
#'      then matches the notational convention used in package \pkg{dynlm}.}
#'   \item{3. Parameter name:}{The name of the regression coefficient, variance,
#'     or covariance specified by the arrow. Assigning the same name to two or
#'     more arrows results in an equality constraint. Specifying the parameter name
#'     as \code{NA} produces a fixed parameter.}
#'   \item{4. Value:}{start value for a free parameter or value of a fixed parameter.
#'     If given as \code{NA} (or simply omitted), the model is provide a default
#'     starting value.}
#' }
#'
#' Lines may end in a comment following #. The function extends code copied from package
#' `sem` under licence GPL (>= 2) with permission from John Fox.
#'
#' \strong{Simultaneous autoregressive process for simultaneous and lagged effects}
#'
#' This text then specifies linkages in a multivariate time-series model for variables \eqn{\mathbf X}
#' with dimensions \eqn{T \times C} for \eqn{T} times and \eqn{C} variables.
#' \code{make_dsem_ram} then parses this text to build a path matrix \eqn{\mathbf{P}} with
#' dimensions \eqn{TC \times TC}, where element \eqn{\rho_{k_2,k_1}}
#' represents the impact of \eqn{x_{t_1,c_1}} on \eqn{x_{t_2,c_2}}, where \eqn{k_1=T c_1+t_1}
#' and \eqn{k_2=T c_2+t_2}.  This path matrix defines a simultaneous equation
#'
#' \deqn{ \mathrm{vec}(\mathbf X) = \mathbf P \mathrm{vec}(\mathbf X) + \mathrm{vec}(\mathbf \Delta)}
#'
#' where \eqn{\mathbf \Delta} is a matrix of exogenous errors with covariance \eqn{\mathbf{V = \Gamma \Gamma}^t},
#' where \eqn{\mathbf \Gamma} is the Cholesky of exogenous covariance.  This
#' simultaneous autoregressive (SAR) process then results in \eqn{\mathbf X} having covariance:
#'
#' \deqn{ \mathrm{Cov}(\mathbf X) = \mathbf{(I - P)}^{-1} \mathbf{\Gamma \Gamma}^t \mathbf{((I - P)}^{-1})^t }
#'
#' Usefully, computing the inverse-covariance (precision) matrix \eqn{\mathbf{Q = V}^{-1}} does not require inverting \eqn{\mathbf{(I - P)}}:
#'
#' \deqn{ \mathbf{Q} = (\mathbf{\Gamma}^{-1} \mathbf{(I - P)})^t \mathbf{\Gamma}^{-1} \mathbf{(I - P)} }
#'
#' \strong{Example: univariate first-order autoregressive model}
#'
#' This simultaneous autoregressive (SAR) process across variables and times
#' allows the user to specify both simutanous effects (effects among variables within
#' year \eqn{T}) and lagged effects (effects among variables among years \eqn{T}).
#' As one example, consider a univariate and first-order autoregressive process where \eqn{T=4}.
#' with independent errors.  This is specified by passing \code{ sem = "X -> X, 1, rho \n X <-> X, 0, sigma" } to \code{make_dsem_ram}.
#' This is then parsed to a RAM:
#'
#' \tabular{rrrrr}{
#'   \strong{heads} \tab \strong{to} \tab \strong{from} \tab \strong{paarameter} \tab \strong{start} \cr
#'   1 \tab 2 \tab 1 \tab 1 \tab <NA>\cr
#'   1 \tab 3 \tab 2 \tab 1 \tab <NA>\cr
#'   1 \tab 4 \tab 3 \tab  1 \tab <NA>\cr
#'   2 \tab 1 \tab 1 \tab 2 \tab <NA>\cr
#'   2 \tab 2 \tab 2 \tab  2 \tab <NA>\cr
#'   2 \tab 3 \tab 3 \tab 2 \tab <NA>\cr
#'   2 \tab 4 \tab 4 \tab 2 \tab <NA>
#' }
#'
#' Rows of this RAM where \code{heads=1} are then interpreted to construct the path matrix \eqn{\mathbf P}, where column "from"
#' in the RAM indicates column number in the matrix, column "to" in the RAM indicates row number in the matrix:
#'
#'     \deqn{ \mathbf P = \begin{bmatrix}
#'     0 & 0 & 0 & 0 \\
#'     \rho & 0 & 0 & 0 \\
#'     0 & \rho & 0 & 0 \\
#'     0 & 0 & \rho & 0\\
#'     \end{bmatrix} }
#'
#' While rows where \code{heads=2} are interpreted to construct the Cholesky of exogenous covariance \eqn{\mathbf \Gamma}
#' and column "parameter" in the RAM associates each nonzero element of those
#' two matrices with an element of a vector of estimated parameters:
#'
#'     \deqn{ \mathbf \Gamma = \begin{bmatrix}
#'     \sigma & 0 & 0 & 0 \\
#'     0 & \sigma & 0 & 0 \\
#'     0 & 0 & \sigma & 0 \\
#'     0 & 0 & 0 & \sigma\\
#'     \end{bmatrix} }
#'
#' with two estimated parameters \eqn{\mathbf \beta = (\rho, \sigma) }. This then results in covariance:
#'
#'     \deqn{ \mathrm{Cov}(\mathbf X) = \sigma^2 \begin{bmatrix}
#'     1      & \rho^1              & \rho^2                        & \rho^3                  \\
#'     \rho^1 & 1 + \rho^2          & \rho^1 (1 + \rho^2)           & \rho^2 (1 + \rho^2)     \\
#'     \rho^2 & \rho^1 (1 + \rho^2) & 1 + \rho^2 + \rho^4           & \rho^1 (1 + \rho^2 + \rho^4)                 \\
#'     \rho^3 & \rho^2 (1 + \rho^2) & \rho^1 (1 + \rho^2 + \rho^4)  & 1 + \rho^2 + \rho^4 + \rho^6 \\
#'     \end{bmatrix} }
#'
#' Which converges on the stationary covariance for an AR1 process for times \eqn{t>>1}:
#'
#'     \deqn{ \mathrm{Cov}(\mathbf X) = \frac{\sigma^2}{1+\rho^2} \begin{bmatrix}
#'     1 & \rho^1 & \rho^2 & \rho^3 \\
#'     \rho^1 & 1 & \rho^1 & \rho^2 \\
#'     \rho^2 & \rho^1 & 1 & \rho^1 \\
#'     \rho^3 & \rho^2 & \rho^1 & 1\\
#'     \end{bmatrix} }
#'
#' except having a lower pointwise variance for the initial times, which arises as a "boundary effect".
#'
#' Similarly, the arrow-and-lag notation can be used to specify a SAR representing
#' a conventional structural equation model (SEM), cross-lagged (a.k.a. vector autoregressive)
#' models (VAR), dynamic factor analysis (DFA), or many other time-series models.
#'
#' @return A reticular action module (RAM) describing dependencies
#'
#' @examples
#' # Univariate AR1
#' sem = "
#'   X -> X, 1, rho
#'   X <-> X, 0, sigma
#' "
#' make_dsem_ram( sem=sem, variables="X", times=1:4 )
#'
#' # Univariate AR2
#' sem = "
#'   X -> X, 1, rho1
#'   X -> X, 2, rho2
#'   X <-> X, 0, sigma
#' "
#' make_dsem_ram( sem=sem, variables="X", times=1:4 )
#'
#' # Bivariate VAR
#' sem = "
#'   X -> X, 1, XtoX
#'   X -> Y, 1, XtoY
#'   Y -> X, 1, YtoX
#'   Y -> Y, 1, YtoY
#'   X <-> X, 0, sdX
#'   Y <-> Y, 0, sdY
#' "
#' make_dsem_ram( sem=sem, variables=c("X","Y"), times=1:4 )
#'
#' # Dynamic factor analysis with one factor and two manifest variables
#' # (specifies a random-walk for the factor, and miniscule residual SD)
#' sem = "
#'   factor -> X, 0, loadings1
#'   factor -> Y, 0, loadings2
#'   factor -> factor, 1, NA, 1
#'   X <-> X, 0, NA, 0.01       # Fix at negligible value
#'   Y <-> Y, 0, NA, 0.01       # Fix at negligible value
#' "
#' make_dsem_ram( sem=sem, variables=c("X","Y","factor"), times=1:4 )
#'
#' # ARIMA(1,1,0)
#' sem = "
#'   factor -> factor, 1, rho1 # AR1 component
#'   X -> X, 1, NA, 1          # Integrated component
#'   factor -> X, 0, NA, 1
#'   X <-> X, 0, NA, 0.01      # Fix at negligible value
#' "
#' make_dsem_ram( sem=sem, variables=c("X","factor"), times=1:4 )
#'
#' # ARIMA(0,0,1)
#' sem = "
#'   factor -> X, 0, NA, 1
#'   factor -> X, 1, rho1     # MA1 component
#'   X <-> X, 0, NA, 0.01     # Fix at negligible value
#' "
#' make_dsem_ram( sem=sem, variables=c("X","factor"), times=1:4 )
#'
#' @export
make_dsem_ram <-
function( sem,
          times,
          variables,
          covs = NULL,
          quiet = FALSE,
          remove_na = TRUE ){
  # Docs : https://roxygen2.r-lib.org/articles/formatting.html

  # MATH CHECK IN ROXYGEN DOCS ABOVE
  if( FALSE ){
    #rho = 0.8
    #sigma = 0.5
    #Rho = Gamma = matrix(0, nrow=4, ncol=4)
    #Rho[cbind(2:4,1:3)] = rho
    #Gamma = I = diag(4)
    #diag(Gamma)[] = sigma
    ## DSEM covariance
    #solve(I-Rho) %*% Gamma %*% t(Gamma) %*% t(solve(I-Rho))
    ## Stated covariance
    #sigma^2 * rbind(
    #  c(1, rho, rho^2, rho^3),
    #  c(rho, 1+rho^2, rho*(1+rho^2), rho^2*(1+rho^2) ),
    #  c(rho^2, rho*(1+rho^2), 1+rho^2+rho^4, rho*(1+rho^2+rho^4) ),
    #  c(rho^3, rho^2*(1+rho^2), rho*(1+rho^2+rho^4), 1+rho^2+rho^4+rho^6 )
    #)
  }

  ####### Error checks
  if( !is.numeric(times) ) stop("`times` must be numeric in `make_dsem_ram`")

  ####### Define local functions
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
  model <- cbind( "path"=model$path, "lag"=model$lag, "name"=model$par, "start"=model$start)

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

  ####### Step 2 -- Make RAM

  # Global stuff
  Q_names = expand.grid( times, variables )
  ram = NULL  # heads, to, from, parameter

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
  P_kk = drop0(sparseMatrix( i=1, j=1, x=0, dims=rep(length(variables)*length(times),2) ))   # Make with a zero
  #P_kk = new("dgCMatrix")
  #P_kk = Matrix()
  #P_kk@Dim <- as.integer(rep(length(variables)*length(times),2))
  #P_kk@p = integer(length(variables)*length(times)+1L)
  #G_kk = P_kk
  for( i in seq_len(nrow(model)) ){
    lag = as.numeric(model[i,2])
    L_tt = sparseMatrix( i = seq(lag+1,length(times)),
                         j = seq(1,length(times)-lag),
                         x = 1,
                         dims = rep(length(times),2) )
    P_jj = sparseMatrix( i = match(model[i,'second'],variables),
                         j = match(model[i,'first'],variables),
                         x = 1,
                         dims = rep(length(variables),2) )
    tmp_kk = kronecker(P_jj, L_tt)
    if(abs(as.numeric(model[i,'direction']))==1){
      P_kk = P_kk + tmp_kk * par.nos[i]
    }else{
      G_kk = G_kk + tmp_kk * par.nos[i]
    }
    #for( t in seq_along(times) ){
    #  # Get index for "from"
    #  from = c( times[t], model[i,'first'] )
    #  from_index = match_row( Q_names, from )
    #  from_index = ifelse( length(from_index)==0, NA, from_index )
    #  # Get index for "to"
    #  to = c( times[t+lag], model[i,'second'] )
    #  to_index = match_row( Q_names, to )
    #  to_index = ifelse( length(to_index)==0, NA, to_index )
    #  ram_new = data.frame( "heads"=abs(as.numeric(model[i,'direction'])), "to"=to_index, "from"=from_index, "parameter"=par.no, "start"=startvalues[i] )
    #  ram = rbind( ram, ram_new )
    #}
  }
  #rownames(ram) = NULL
  #f = \(x) sapply(mat2triplet(drop0(x)),cbind)
  f = \(x) matrix(unlist(mat2triplet(x)),ncol=3)
  ram = rbind( cbind(1, f(P_kk)),
               cbind(2, f(G_kk)) )
  ram = data.frame( ram, startvalues[ram[,4]] )
  colnames(ram) = c("heads", "to", "from", "parameter", "start")

  #
  if( isTRUE(remove_na) ){
    which_keep = which(apply( ram[,1:4], MARGIN=1, FUN=\(x)!any(is.na(x)) ))
    ram = ram[ which_keep, ]
  }

  #
  out = list( "model"=model,
              "ram"=ram,
              "variables" = variables,
              "times" = times )
  class(out) = "dsem_ram"
  return(out)
}
