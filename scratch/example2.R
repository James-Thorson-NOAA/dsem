  library(Matrix)
  library(RTMB)

  # Load minimal example
  example = readRDS( file=file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\scratch)',"example.RDS"))
  attach(example)
  set.seed(101)
  z = rnorm(nrow(Rho_kk))

  # Define function with both versions
  f = function(p){
    #
    sparse_solve = function(x){
      invx = solve(x)
      if( RTMB:::ad_context() ){
        out = sparseMatrix(
                      i = row(invx),
                      j = col(invx),
                      x = 1,
                 )
        out = AD(out)
        out@x = invx
        #out = drop0(out)    # drop0 doesn't work
        return(out)
      }else{
        return(invx)
      }
    }

    Gamma_kk = AD( p[1] * Gamma_kk )
    Rho_kk = AD( p[2] * Rho_kk )
    V_kk = t(Gamma_kk) %*% Gamma_kk
    IminusRho_kk = Diagonal(nrow(Rho_kk)) - Rho_kk

    # Option-1 ... works
    invV_kk = Gamma_kk
    invV_kk@x = 1 / Gamma_kk@x^2
    Q1_kk = t(IminusRho_kk) %*% invV_kk %*% IminusRho_kk
    REPORT( Q1_kk )

    # Option-2
    invV_kk = solve(V_kk)
    Q2_kk = t(IminusRho_kk) %*% invV_kk %*% IminusRho_kk
    REPORT( Q2_kk )

    # Option-3
    invV_kk = solve(V_kk)
    #if( RTMB:::ad_context() ){
    #  temp_kk = sparseMatrix(
    #                i = row(invV_kk),
    #                j = col(invV_kk),
    #                x = 1,
    #           )
    #  temp_kk = AD(temp_kk)
    #  temp_kk@x = invV_kk
    #  invV_kk = temp_kk
    #}
    invV_kk = sparse_solve(V_kk)
    Q3_kk = t(IminusRho_kk) %*% invV_kk %*% IminusRho_kk
    REPORT( Q3_kk )

    # Option-1 WORKS for diagonal matrix only
    #jnll_gmrf = -1 * dgmrf( z, mu=rep(0,nrow(Rho_kk)), Q=Q1_kk, log=TRUE )

    # Option-2 DOES NOT WORK
    #jnll_gmrf = -1 * dgmrf( z, mu=rep(0,nrow(Rho_kk)), Q=Q2_kk, log=TRUE )

    # Option-3 WORKS generally, but requires a manual coersion from dense to sparse matrix
    jnll_gmrf = -1 * dgmrf( z, mu=rep(0,nrow(Rho_kk)), Q=Q3_kk, log=TRUE )
    return(jnll_gmrf)
  }

  # Works for Option-1
  f(c(1,1))
  obj = MakeADFun( f, c(1,1) )
  obj$fn(obj$par)

  # Show they're identical if compiled using Option-1
  Rep = obj$report()
  identical( Rep$Q1_kk, Rep$Q2_kk )
  identical( Rep$Q1_kk, Rep$Q3_kk )
