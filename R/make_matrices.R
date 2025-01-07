make_matrices <-
function( beta_p,
          model,
          times,
          variables ){

  if(missing(beta_p)){
    model_unique = model[match(unique(model$parameter),model$parameter),]
    beta_p =  as.numeric(model_unique$start)
  }

  # Loop through paths
  P_kk = drop0(sparseMatrix( i=1, j=1, x=0, dims=rep(length(variables)*length(times),2) ))   # Make with a zero
  #P_kk = AD(P_kk)
  G_kk = (P_kk)
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

    # Assemble
    if(abs(as.numeric(model[i,'direction']))==1){
      tmp_kk = (kronecker(P_jj, L_tt))
      P_kk = P_kk + beta_p[model$parameter[i]] * tmp_kk # AD(tmp_kk)
    }else{
      tmp_kk = (kronecker(P_jj, L_tt))
      G_kk = G_kk + beta_p[model$parameter[i]] * tmp_kk # AD(tmp_kk)
    }
  }

  # Diagonal component
  I_kk = Diagonal(nrow(P_kk))

  # Assemble
  IminusP_kk = AD(I_kk - P_kk)
  invV_kk = AD(G_kk)
  invV_kk@x = 1 / G_kk@x^2
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
