  library(Matrix)
  library(RTMB)

  # Load minimal example
  example = readRDS( file=file.path(R'(C:\Users\James.Thorson\Desktop\Git\dsem\scratch)',"example.RDS"))
  attach(example)
  set.seed(101)
  z = rnorm(nrow(Rho_kk))

  # Define function with both versions
  f = function(p){
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
    #trip = mat2triplet(invV_kk)
    #print(trip)
    #print(invV_kk@i)
    #print(invV_kk@p)
    #print(invV_kk@x)
    #p = invV_kk@p
    #n.p = length(p)
    #dp <- p[-1L] - p[-n.p]
    #n.dp = length(dp)
    #j = rep.int(seq.int(from = 0L, length.out = n.dp), dp) + 1
    #print(j)
    #invV_kk = sparseMatrix(
    #            i = invV_kk@i + 1,
    #            p = invV_kk@p,
    #            #j = j,
    #            x = invV_kk@x
    #          )
    if( RTMB:::ad_context() ){
      #print(invV_kk[1:10,1:10])
      #class(invV_kk)
      #print(row(invV_kk))

      #vec = as.matrix(invV_kk)[1,1]
      #print(vec)

      #Dim = rep(nrow(V_kk),2)
      #invV2_kk = drop0(sparseMatrix( i=1, j=1, x=0, dims=Dim ))   # Make with a zero
      #for(i in seq_len(nrow(invV_kk)) ){
      #for(j in seq_len(ncol(invV_kk)) ){
      #  Ind_kk = sparseMatrix(
      #              i = i,
      #              j = j,
      #              x = 1,
      #              dims = Dim
      #            )
      #  invV2_kk = invV2_kk + invV_kk[i,j] * Ind_kk
      #}}
      invV2_kk = sparseMatrix(
                    i = row(invV_kk),
                    j = col(invV_kk),
                    x = 1,
               )
      invV2_kk = AD(invV2_kk)
      invV2_kk@x = invV_kk
      REPORT(invV2_kk)
    }else{
      invV2_kk = invV_kk
    }
    #invV2_kk@x = invV_kk@x
    #invV_kk = AD(invV_kk)
    #print(length(invV2_kk@x))
    #print(length(invV_kk@x))
    #print(invV2_kk)
    Q2_kk = t(IminusRho_kk) %*% invV2_kk %*% IminusRho_kk
    REPORT( Q2_kk )
    #print(Q2_kk)

    # Option-3 ... try rebuilding
    #Q3_kk = AD(sparseMatrix(i=Q2_kk@i, p=Q2_kk@p, x=Q2_kk@x))
    #REPORT( Q3_kk )

    # Option-1 WORKS
    #jnll_gmrf = -1 * dgmrf( z, mu=rep(0,nrow(Rho_kk)), Q=Q1_kk, log=TRUE )

    # Option-2 DOES NOT WORK
    jnll_gmrf = -1 * dgmrf( z, mu=rep(0,nrow(Rho_kk)), Q=Q2_kk, log=TRUE )
    #return(jnll_gmrf)
  }

  # Works for Option-1
  f(c(1,1))
  obj = MakeADFun( f, c(1,1) )
#  obj$fn(obj$par)

  # Show they're identical if compiled using Option-1
#  Rep = obj$report()
#  identical( Rep$Q1_kk, Rep$Q2_kk )
#  attr(Rep$Q1_kk,"Dim")
