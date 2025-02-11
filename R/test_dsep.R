
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

# Modified from phylopath
set_to_paths <-
function(x) {

  dep <- x[2]
  ind <- x[1]
  cond <- x[c(-1, -2)]

  out = paste0( ind, " -> ", dep, ", 0, target_", dep, "_", ind )
  for( i in seq_along(cond)){
    row = paste0( cond[i], " -> ", dep, ", 0, nuissance_", dep, "_", cond[i] )
    #out = paste0( out, " \n ", row )
    out = c( out, row )
  }
  return(out)
}

# Modified from phylopath
find_paths <-
function( A,
          order ) {

  s <- basiSet(A[order,order])
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

get_P <-
function( object,
          max_lag ){

  Summary = summary(object)
  times = seq_len( max_lag + 1 )
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

convert_path <-
function( path ){

  arrow_and_lag = path
  terms = strsplit(arrow_and_lag, ", " )[[1]]
  term_one = strsplit(terms[1], " -> " )[[1]]
  lag_one = sapply( term_one, \(x){strsplit(x, "_lag")[[1]][2]} )
  term_one = sapply( term_one, \(x){strsplit(x, "_lag")[[1]][1]} )
  terms[1] = paste0( term_one, collapse = " -> " )
  terms[2] = diff(as.numeric(lag_one))
  arrow_and_lag = paste0( terms, collapse = ", " )
  return(arrow_and_lag)
}

fit_dsem <-
function( object,
          sem ){

  control = object$internal$control
  fit = dsem( sem = paste0(sem," \n "),
                tsdata = object$internal$tsdata,
                family = object$internal$family,
                estimate_delta0 = object$internal$estimate_delta0,
                control = dsem_control( use_REML = control$use_REML,
                                        quiet = TRUE ) )
  return(fit)
}

#' @title Test d-separation
#'
#' @description
#' Calculate the p-value for a test of d-separation \strong{(Experimental)}
#'
#' @param object object from \code{\link{dsem}}
#' @param max_lag how many lags to include when defining the set of conditional
#'        independence relationships. If missing, this value is taken from
#'        the maximum lag that's included in the model.
#' @param what whether to just get the p-value, or a named list with the p-value,
#' @param conditional_test whether to test each conditional-independence relationship
#'        using a (univariate) Wald test or a (multivariate) likelihood ratio test
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
          max_lag = NULL,
          what = c("pvalue","CIC","all"),
          conditional_test = c("wald","lr") ){

  # Get amat
  what = match.arg(what)
  conditional_test = match.arg(conditional_test)
  if(is.null(max_lag)){
    max_lag = max( summary(object)$lag )
  }
  P_kk = get_P( object = object, max_lag = max_lag )
  d = t(as.matrix(P_kk))
  # replace with adjacency matrix
  A = ifelse( d==0, 0, 1)

  # find_formulas
  order = find_consensus_order(list(A))
  paths = find_paths( A, order=order )
  arrow_and_lags = lapply(paths, set_to_paths)
  sems = lapply( arrow_and_lags, function(x){ sapply(x, function(y){convert_path(y)}, USE.NAMES=FALSE) } )
  fits = lapply( sems, fit_dsem, object=object )

  if( conditional_test == "lr" ){
    # eliminate target variable and refit
    sems_null = lapply( sems, function(vec){vec[-1]} )
    fits_null = lapply( sems_null, fit_dsem, object=object )
    # Compare objectives as likelihood ratio test
    objectives = sapply( fits, function(l) l$opt$obj )
    objectives_null = sapply( fits_null, function(l)l$opt$obj )
    pvalues = 1 - pchisq( 2*objectives_null - 2*objectives, df = 1 )
  }else{
    summaries = lapply( fits, summary )
    pvalues = sapply( summaries, function(l) l[1,'p_value'] )
  }

  #
  C_p = -2 * sum( log(pvalues) )
  pvalue = 1 - pchisq( C_p, df = 2 * length(pvalues) )
  CIC = C_p + 2 * length(object$opt$par)
  if( what == "pvalue" ){
    return(pvalue)
  }else if( what == "CIC" ){
    return(CIC)
  }else{
    return(list(
      pvalue = pvalue,
      CIC = CIC,
      #paths = paths,
      arrow_and_lags = arrow_and_lags,
      sems = sems,
      fits = fits,
      pvalues = pvalues
    ))
  }
}
