
partition_variance <-
function( object,
          which_pred,
          which_response,
          n_lags = 10 ){

  # Unpack stuff
  Z = object$internal$tsdata
  if(is.null(object$internal$parhat)){
    object$internal$parhat = object$obj$env$parList()
  }

  # Extract path matrix
  matrices = make_matrices(
    beta_p = object$internal$parhat$beta,
    model = object$sem_full,
    times = seq_len(n_lags),
    variables = colnames(Z)
  )
  out = expand.grid(lag = seq_len(n_lags) - 1, variable = colnames(Z) )

  #
  IminusP_kk = matrices$IminusP_kk
  invIminusP_kk = Matrix::solve(IminusP_kk)
  G0_kk = G_kk = matrices$G_kk

  # Extract variance for fitted model
  V_kk = Matrix::t(G_kk) %*% G_kk
  Sigma1_kk = invIminusP_kk %*% V_kk %*% Matrix::t(invIminusP_kk)

  # Zero out variances
  match_cols = which( out$variable %in% which_pred )
  G0_kk[,-match_cols] = 0
  V0_kk = Matrix::t(G0_kk) %*% G0_kk
  Sigma0_kk = invIminusP_kk %*% V0_kk %*% Matrix::t(invIminusP_kk)

  #
  match_vals = which( out$variable %in% which_response )
  var_ratio = Matrix::diag(Sigma0_kk)[match_vals] / Matrix::diag(Sigma1_kk)[match_vals]
  return(var_ratio)
}

library(dsem)

##########
# Linear model example
##########

x = rnorm(100)
y = 1 + 1 * x + rnorm(100)
data = data.frame(x=x, y=y)

# Fit as linear model
Lm = lm( y ~ x, data=data )

# Fit as DSEM
fit = dsem( sem = "x -> y, 0, beta",
            tsdata = ts(data),
            control = dsem_control(quiet=TRUE) )

#
partition_variance( fit,
                    which_pred = "x",
                    which_response = "y",
                    n_lags = 10 )

#################
# Klein example
#################
data(KleinI, package="AER")
TS = ts(data.frame(KleinI, "time"=time(KleinI) - 1931))

# Specify by declaring each arrow and lag
sem = "
  # Link, lag, param_name
  cprofits -> consumption, 0, a1
  cprofits -> consumption, 1, a2
  pwage -> consumption, 0, a3
  gwage -> consumption, 0, a3

  cprofits -> invest, 0, b1
  cprofits -> invest, 1, b2
  capital -> invest, 0, b3

  gnp -> pwage, 0, c2
  gnp -> pwage, 1, c3
  time -> pwage, 0, c1
"
tsdata = TS[,c("time","gnp","pwage","cprofits",'consumption',
               "gwage","invest","capital")]
fit = dsem( sem=sem,
            tsdata = tsdata,
            estimate_delta0 = TRUE,
            control = dsem_control(
              quiet = TRUE,
              newton_loops = 0) )

#
partition_variance( fit,
                    which_pred = "gnp",
                    which_response = "pwage",
                    n_lags = 10 )
