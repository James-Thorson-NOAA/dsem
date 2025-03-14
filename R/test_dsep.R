
get_P <-
function( object,
          times ){

  Summary = summary(object)
  variables = colnames(object$internal$tsdata)
  P_kk = make_matrices(
    beta_p = object$internal$parhat$beta,
    model = object$sem_full,
    times = times,
    variables = variables
  )$P_kk
  varnames = expand.grid( times - 1, variables )
  varnames = paste0( varnames[,2], "_lag", varnames[,1] )
  dimnames(P_kk) = list(to=varnames, from=varnames)
  return(P_kk)
}

# Modified from phylopath
find_consensus_order <-
function(model_set) {

  # If the fully combined model is acyclic, then we use that.
  model_set_same_order <- lapply(model_set, function(x) {
    x[rownames(model_set[[1]]), colnames(model_set[[1]])]
  } )
  full_model <- sign(Reduce('+', model_set_same_order))
  if (isAcyclic(full_model)) {
    return(rownames(topSort(full_model)))
  }
  # Otherwise we find the most common orderings and use those.
  # Make sure all models are ordered:
  model_set <- lapply(model_set, topSort)
  vars <- lapply(model_set, colnames)
  combs <- as.data.frame(t(combn(vars[[1]], 2)), stringsAsFactors = FALSE)
  names(combs) <- c('node1', 'node2')
  combs$count <- 0
  for (i in seq_along(vars)) {
    v <- apply(combs, 1, function(x) {
      which(vars[[i]] == x[1]) < which(vars[[i]] == x[2])
    } )
    combs$count <- combs$count + v
  }

  # If node1 is commonly ordered above node2, leave as is, otherwise swap them around
  tmp <- combs$node1
  combs$node1 <- ifelse(combs$count > 0.5 * length(model_set), combs$node1, combs$node2)
  combs$node2 <- ifelse(combs$count > 0.5 * length(model_set), combs$node2, tmp)

  # Now we order the nodes by how many nodes they are above, this should go from n:1
  combs$n <- table(combs$node1)[combs$node1]
  combs <- combs[order(-combs$n), ]
  res <- unlist(c(unique(combs$node1), utils::tail(combs, 1)[, 'node2']))
  names(res) <- NULL
  res
}

#basiSet <-
#function( amat ){
#  amat <- topSort(amat)
#  nod <- rownames(amat)
#  dv <- length(nod)
#  ind <- NULL
#  for (r in 1:dv) {
#  for (s in r:dv) {
#    if ((amat[r, s] != 0) | (s == r)) {
#      next
#    }
#    else{
#      ed <- nod[c(r, s)]
#      pa.r <- nod[amat[, r] == 1]
#      pa.s <- nod[amat[, s] == 1]
#      dsep <- union(pa.r, pa.s)
#      dsep <- setdiff(dsep, ed)
#      b <- list(c(ed, dsep))
#      ind <- c(ind, b)
#    }
#  }}
#  ind
#}

# Modified from phylopath
find_paths <-
function( A,
          order ) {

  s <- basiSet(A) #[order,order])
  #s <- basiSet(A)
  if (is.null(s)) {
    stop('One or some of your models are fully connected, and cannot be tested.')
  }
  s <- lapply(s, function(x) {
    # define whether there are existing paths between the two nodes in both directions.
    path1 <- !is.null(findPath(A, which(rownames(A) == x[1]), which(rownames(A) == x[2])))
    path2 <- !is.null(findPath(A, which(rownames(A) == x[2]), which(rownames(A) == x[1])))
    if (path1 & !path2) {
      # the first vertex is upstream, so we do not re-order
      return(x)
    }
    if ((path2 & !path1) | (path1 & path2)) {
      # these conditions should not occur, the first means basiSet is returning the wrong order,
      # the second should only occur if there are cycles.
      stop('If you get this error, please contact the maintainer.')
    }
    if (!path1 & !path2) {
      # check whether the order is according to `order`
      if (which(order == x[1]) < which(order == x[2])) {
        return(x)
      } else {
        return(c(x[2], x[1], x[-(1:2)]))
      }
    }
  } )
  return(s)
}

remove_paths <-
function( paths,
          n_burnin ){

  out = NULL
  for( i in seq_along(paths) ){
    lags = sapply( paths[[i]], function(char) as.numeric(strsplit(char,"_lag")[[1]][2]) )
    if( all(lags[1:2] >= n_burnin) ){
      out = c( out, paths[i] )
    }
  }
  return(out)
}

# Modified from phylopath
path_to_arrow <-
function(x) {

  dep <- x[2]
  ind <- x[1]
  cond <- x[c(-1, -2)]

  out = paste0( ind, " -> ", dep, ", 0, target" )
  for( i in seq_along(cond)){
    row = paste0( cond[i], " -> ", dep, ", 0, nuissance_", i )
    #out = paste0( out, " \n ", row )
    out = c( out, row )
  }
  return(out)
}

convert_path <-
function( path ){

  arrow_and_lag = path
  terms = strsplit(arrow_and_lag, ", " )[[1]]
  term_one = strsplit(terms[1], " -> " )[[1]]
  lag_one = sapply( term_one, function(x){strsplit(x, "_lag")[[1]][2]} )
  term_one = sapply( term_one, function(x){strsplit(x, "_lag")[[1]][1]} )
  terms[1] = paste0( term_one, collapse = " -> " )
  terms[2] = diff(as.numeric(lag_one))
  arrow_and_lag = paste0( terms, collapse = ", " )
  return(arrow_and_lag)
}

fit_dsem <-
function( object,
          sem,
          impute_data = TRUE,
          seed,
          tsdata,
          getsd = TRUE ){

  # Modify controls
  control = object$internal$control
    control$quiet = TRUE
    control$extra_convergence_checks = TRUE
    control$newton_loops = 0
    control$getsd = getsd

  if( impute_data == "none" ){
    tsdata = object$internal$tsdata
  }
  if( impute_data == "by_test" ){
    # Simulate random effects from joint precision, and measurement errors from states
    tsdata = simulate( object,
                       variance = ifelse(length(object$obj$env$random)==0,"none","random"),
                       seed = seed,
                       fill_missing = TRUE )[[1]]
  }
  if( impute_data == "single" ){
    # use tsdata passed from test_dsep
  }

  # Refit
  if( inherits(object,"dsemRTMB") ){
    fit = try(dsemRTMB( sem = paste0(sem, collapse=" \n "),
                  tsdata = tsdata,
                  family = object$internal$family,
                  estimate_delta0 = object$internal$estimate_delta0,
                  log_prior = object$internal$log_prior,
                  control = control ))
  }else{
    fit = try(dsem( sem = paste0(sem, collapse=" \n "),
                  tsdata = tsdata,
                  family = object$internal$family,
                  estimate_delta0 = object$internal$estimate_delta0,
                  prior_negloglike = object$internal$prior_negloglike,
                  control = control ))
  }
  return(fit)
}

#' @title Test d-separation
#'
#' @description
#' Calculate the p-value for a test of d-separation \strong{(Experimental)}
#'
#' @param object object from \code{\link{dsem}}
#' @param n_time how many times to include when defining the set of conditional
#'        independence relationships. If missing, this value is taken from
#'        the maximum lag that's included in the model plus one.
#' @param n_burnin how many times to include prior to \code{seq_len(n_time)} when
#'        identifying the conditioning set that must be included when defining
#'        conditional independence relationships.
#' @param what whether to just get the p-value, an information criterion
#'        based on the conditional independence test, or a named list with these two
#'        and other intermediate calculations (used for diagnosing test behavior)
#' @param test whether to test each conditional-independence relationship
#'        using a (univariate) wald test or a (multivariate) likelihood ratio test.
#'        The likelihood-ratio test might be more accurate given estimation covariance
#'        and also faster (does not require standard errors), but also is not
#'        used by phylopath and therefore less supported by previous d-dsep
#'        testing applications.
#' @param impute_data whether to independently impute missing data for each
#'        conditional independence test, or to use imputed values from the original
#'        fit.  The data are imputed separately for each conditional independence
#'        test, so that they are uncorrelated as expected when combining them
#'        using Fisher's method.  Preliminary testing suggests
#'        that using imputed data improves test performance
#' @param order an optional character vector providing the order for variables to be
#'        tested when defining the directed acyclic graph for use in d-sep testing
#' @param seed random number seed used when simulating imputed data, so that
#'        results are reproducible.
#'
#' @details
#' A user-specified SEM implies a set of conditional independence relationships
#' among variables, which can be fitted individually, extracting the
#' slope and associated p-value, and then combining these p-values to define
#' a model-wide (omnibus) p-value for the hypothesis that a given data set arises
#' from the specified model.  This test is modified from package:phylopath.  However
#' it is unclear exactly how to define the set of conditional-independence assumptions
#' in a model with temporal autocorrelation, and the test was not developed for
#' uses when data are missing.  At the time of writing, the function is hightly
#' experimental.
#'
#' Note that the method is not currently designed to deal with two-headed arrows
#' among variables (i.e., exogenous covariance).
#'
#' @return
#' A p-value representing the weight of evidence that the data arises
#' from the specified model, where a low p-value indicates
#' significant evidence for rejecting this hypothesis.
#'
#' @references
#' Shipley, B. (2000). A new inferential test
#' for path models based on directed acyclic graphs. Structural
#'   Equation Modeling, 7(2), 206-218. \doi{10.1207/S15328007SEM0702_4}
#'
#' @examples
#' # Simulate data set
#' set.seed(101)
#' a = rnorm( 100 )
#' b = 0.5*a + rnorm(100)
#' c = 1*a + rnorm(100)
#' d = 1*b - 0.5*c + rnorm(100)
#' tsdata = ts(data.frame(a=a, b=b, c=c, d=d))
#'
#' # fit wrong model
#' wrong = dsem(
#'   tsdata = tsdata,
#'   sem = "
#'     a -> d, 0, a_to_d
#'     b -> d, 0, b_to_d
#'     c -> d, 0, c_to_d
#'   "
#' )
#' test_dsep( wrong )
#'
#' # fit right model
#' right = dsem(
#'   tsdata = tsdata,
#'   sem = "
#'     a -> b, 0, a_to_b
#'     a -> c, 0, a_to_c
#'     b -> d, 0, b_to_d
#'     c -> d, 0, c_to_d
#'   "
#' )
#' test_dsep( right )
#' @export
test_dsep <-
function( object,
          n_time = NULL,
          n_burnin = NULL,
          what = c("pvalue","CIC","all"),
          test = c("wald","lr"),
          seed = 123456,
          order = NULL,
          impute_data = c("by_test","single","none") ){

  # Check inputs
  what = match.arg(what)
  test = match.arg(test)
  impute_data = match.arg(impute_data)
  if( test=="lr" & isTRUE(impute_data) ) stop("LR test is not designed to work when imputing data")
  out = list( n_time = n_time,
              n_burnin = n_burnin )

  # TO CHECK:
  # Using a single imputed data set for all conditional independence tests is a
  # bad idea because it induces a correlation among tests, which the Fisher method ignores,
  # such that p-values are skewed towards zero even for the right model
  #out$tsdata = object$internal$tsdata
  if( impute_data == "single" ){
    # Simulate random effects from joint precision, and measurement errors from states
    out$tsdata = simulate( object,
                       variance = ifelse(length(object$obj$env$random)==0,"none","random"),
                       fill_missing = TRUE,
                       seed = seed )[[1]]
  }

  # Detect n_time and n_burnin
  if( is.null(out$n_burnin) ){
    if( is.null(out$n_time) ){
      out$n_burnin = max( summary(object)$lag )
    }else{
      out$n_burnin = 0
    }
  }
  if(is.null(out$n_time)){
    out$n_time = max( summary(object)$lag ) + 1
  }else{
    warning( "Please check `n_time` carefully")
  }
  times = seq_len(out$n_time + out$n_burnin)

  # Get adjacency
  P_kk = get_P( object = object, times = times )
  d = t(as.matrix(P_kk))
  # replace with adjacency matrix
  A = ifelse( d==0, 0, 1)

  # Re-order
  if( is.null(order) ){
    out$order = find_consensus_order(list(A))
  }else{
    out$order = order
  }
  # Find paths
  out$paths = find_paths( A,
                      order = out$order )
  # Remove paths from initial "burn-in" buffer in time-stepping
  out$paths = remove_paths( out$paths,
                            n_burnin = out$n_burnin )
  if(length(out$paths)==0){
    stop("No paths remain. The model appears to have no conditional independence relationships to test.")
  }
  # convert to arrow-and-lag notation in variable-by-lag
  out$arrows = lapply( out$paths,
                       FUN = path_to_arrow)
  # Convert to SEM notation
  out$sems = lapply( out$arrows,
                 FUN = function(x){ sapply( x,
                                            FUN = function(y){convert_path(y)},
                                            USE.NAMES = FALSE) } )
  # Eliminate duplicates
  out$sems = unique(out$sems)
  # Fit models
  #out$fits = lapply( out$sems,
  #                   FUN = fit_dsem,
  #                   object = object,
  #                   getsd = ifelse( test=="lr", FALSE, TRUE),
  #                   #seed = seq_along(out$sems_null),
  #                   impute_data = impute_data )
  for(i in seq_along(out$sems)){
    out$fits[[i]] = fit_dsem(
              object = object,
              sem=out$sems[[i]],
              impute_data = impute_data,
              tsdata = out$tsdata,
              seed = seed + i,
              getsd = ifelse( test=="lr", FALSE, TRUE) )
  }
  # out$fits[[2]]$obj$env$data$y_tj

  if( test == "lr" ){
    # eliminate target variable and refit
    # Requires a fixed seed for each paired comparison of out$fits and out$fits_null
    out$sems_null = lapply( out$sems,
                        FUN = function(vec){vec[-1]} )
    #out$fits_null = lapply( out$sems_null,
    #                        FUN = fit_dsem,
    #                        object = object,
    #                        getsd = FALSE,
    #                        #seed = seq_along(out$sems_null),
    #                        impute_data = impute_data )
    for(i in seq_along(out$sems_null)){
      out$fits[[i]] = fit_dsem(
                object = object,
                sem=out$sems_null[[i]],
                impute_data = impute_data,
                tsdata = out$tsdata,
                seed = seed + i,
                getsd = FALSE )
    }
    # Compare objectives as likelihood ratio test (chi-squared with 1 degree of freedom)
    objectives = sapply( out$fits,
                         FUN = function(list) list$opt$obj )
    objectives_null = sapply( out$fits_null,
                              FUN = function(list) list$opt$obj )
    out$pvalues = 1 - pchisq( 2*objectives_null - 2*objectives, df = 1 )
  }else{
    summaries = lapply( out$fits, summary )
    get_pvalue = function(l) if(is.null(nrow(l))){NA}else{l[1,'p_value']}
    out$pvalues = sapply( summaries,
                          FUN = get_pvalue )
  }
  #for(i in seq_along(summaries)) get_pvalue( summaries[[i]] )

  # Compute test statistics
  C_p = -2 * sum( log(out$pvalues) )
  # Fisher's method: sum(log(pvalues)) follows chi-squared distribution with 2*length(pvalues) degrees of freedom
  out$pvalue = 1 - pchisq( C_p, df = 2 * length(out$pvalues) )
  out$CIC = C_p + 2 * length(object$opt$par)

  # Returns
  if( what == "pvalue" ){
    out = out$pvalue
  }else if( what == "CIC" ){
    out = out$CIC
  }
  return(out)
}
