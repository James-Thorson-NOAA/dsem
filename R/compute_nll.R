# model function
compute_nll <-
function( parlist,
          model,
          y_tj,
          family,
          options,
          log_prior,
          simulate_data = FALSE,
          simulate_gmrf = FALSE ) {
  # options[1] -> 0: full rank;  1: rank-reduced GMRF
  # options[2] -> 0: constant conditional variance;  1: constant marginal variance

  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  # Temporary fix for solve(adsparse) returning dense-matrix
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

  # Unpack parameters explicitly
  beta_z = parlist$beta_z
  delta0_j = parlist$delta0_j
  mu_j = parlist$mu_j
  sigma_j = exp( parlist$lnsigma_j )
  x_tj = parlist$x_tj

  #
  n_z = length(unique(model$parameter))
  #n_p2 = length(unique(subset(model,direction==2)$parameter))
  #n_p1 = length(unique(subset(model,direction==1)$parameter))
  n_t = nrow(y_tj)
  n_j = ncol(y_tj)
  n_k = prod(dim(y_tj))

  # Unpack
  #model_unique = model[match(unique(model$parameter),model$parameter),]
  #beta_p[which(model_unique$direction==1)] = beta_p1
  #beta_p[which(model_unique$direction==2)] = beta_p2

  # Build matrices
  ram = make_matrices(
            beta_p = beta_z,
            model = model,
            times = as.numeric(time(y_tj)),
            variables = colnames(y_tj) )

  # Assemble
  Rho_kk = AD(ram$P_kk)
  IminusRho_kk = Diagonal(n_k) - Rho_kk
  # Assemble variance
  Gamma_kk = AD(ram$G_kk)
  V_kk = t(Gamma_kk) %*% Gamma_kk

  # Rescale I-Rho and Gamma if using constant marginal variance options
  if( (options[2]==1) || (options[2]==2) ){
    invIminusRho_kk = sparse_solve(IminusRho_kk)     # solve(adsparse) returns dense-matrix
    #print(class(invIminusRho_kk))

    # Hadamard squared LU-decomposition
    # See: https://eigen.tuxfamily.org/dox/group__QuickRefPage.html
    squared_invIminusRho_kk = invIminusRho_kk
    squared_invIminusRho_kk@x = squared_invIminusRho_kk@x^2

    if( options[2] == 1 ){
      # 1-matrix
      ones_k1 = matrix(1, nrow=n_k, ncol=1)

      # Calculate diag( t(Gamma) * Gamma )
      squared_Gamma_kk = Gamma_kk
      squared_Gamma_kk@x = squared_Gamma_kk@x^2
      sigma2_k1 = t(squared_Gamma_kk) %*% ones_k1;

      # Rowsums
      margvar_k = solve(squared_invIminusRho_kk, sigma2_k1)

      # Rescale IminusRho_kk and Gamma
      invmargsd_kk = invsigma_kk = AD(Diagonal(n_k))
      invmargsd_kk@x = 1 / sqrt(margvar_k)
      invsigma_kk@x = 1 / sqrt(sigma2_k1)
      IminusRho_kk = invmargsd_kk %*% IminusRho_kk;
      Gamma_kk = invsigma_kk %*% Gamma_kk;
    }else{
      # calculate diag(Gamma)^2
      targetvar_k = diag(Gamma_kk)^2

      # Rescale Gamma
      margvar_k = solve(squared_invIminusRho_kk, targetvar_k)
      diag(Gamma_kk) = sqrt(margvar_k)
    }
  }

  # Calculate effect of initial condition -- SPARSE version
  delta_tj = matrix( 0, nrow=n_t, ncol=n_j )
  if( length(delta0_j)>0 ){
    delta_tj[1,] = delta0_j
    delta_k = solve( IminusRho_kk, as.vector(delta_tj) )
    delta_tj[] = delta_k
  }

  #
  xhat_tj = outer( rep(1,n_t), mu_j )

  # Probability of GMRF
  if( options[1] == 0 ){
    # Full rank GMRF
    z_tj = x_tj

    # Works for diagonal
    #invV_kk = AD(ram$G_kk)
    #invV_kk@x = 1 / Gamma_kk@x^2
    #Q_kk = t(IminusRho_kk) %*% invV_kk %*% IminusRho_kk

    # Works in general
    invV_kk = sparse_solve( V_kk )
    Q_kk = t(IminusRho_kk) %*% invV_kk %*% IminusRho_kk

    # Fine from here
    jnll_gmrf = -1 * dgmrf( as.vector(z_tj), mu=as.vector(xhat_tj + delta_tj), Q=Q_kk, log=TRUE )
    REPORT( Q_kk )

    # Only simulate GMRF if also simulating new data
    if( isTRUE(simulate_data) & isTRUE(simulate_gmrf) ){
      x_tj[] = z_tj[] = rgmrf( mu=as.vector(xhat_tj + delta_tj), Q=Q_kk )
    }
  }else{
    # Reduced rank projection .. dgmrf is lower precision than GMRF in CPP
    jnll_gmrf = -1 * sum( dnorm(x_tj, mean=0, sd=1, log=TRUE) )

    # Only simulate GMRF if also simulating new data
    if( isTRUE(simulate_data) & isTRUE(simulate_gmrf) ){
      x_tj[] = rnorm(n=prod(dim(x_tj)), mean=0, sd=1)
    }

    #
    z_k1 = as.vector(x_tj)
    z_k2 = Gamma_kk %*% z_k1
    z_k3 = solve(IminusRho_kk, z_k2)
    z_tj = matrix(as.vector(z_k3), nrow=n_t, ncol=n_j) + xhat_tj + delta_tj
    REPORT( z_k1 )
    REPORT( z_k2 )
    REPORT( z_k3 )
  }

  # Likelihood of data | random effects
  loglik_tj = mu_tj = devresid_tj = matrix( 0, nrow=n_t, ncol=n_j )
  pow = function(a,b) a^b
  for( t in 1:n_t ){
  for( j in 1:n_j ){
    # familycode = 0 :  don't include likelihood
    if( family[j]=="fixed" ){
      mu_tj[t,j] = z_tj[t,j];
      if(!is.na(y_tj[t,j])){
        if( isTRUE(simulate_data) ){
          y_tj[t,j] = mu_tj[t,j];
        }
      }
      devresid_tj[t,j] = 0;
    }
    # familycode = 1 :  normal
    if( family[j]=="normal" ){
      mu_tj[t,j] = z_tj[t,j];
      if(!is.na(y_tj[t,j])){
        loglik_tj[t,j] = dnorm( y_tj[t,j], mu_tj[t,j], sigma_j[j], TRUE );
        if( isTRUE(simulate_data) ){
          y_tj[t,j] = rnorm( n=1, mean=mu_tj[t,j], sd=sigma_j[j] )
        }
      }
      devresid_tj[t,j] = y_tj[t,j] - mu_tj[t,j];
    }
    # familycode = 2 :  binomial
    if( family[j]=="binomial" ){
      mu_tj[t,j] = invlogit(z_tj[t,j]);
      if(!is.na(y_tj[t,j])){
        loglik_tj[t,j] = dbinom( y_tj[t,j], Type(1.0), mu_tj[t,j], TRUE );
        if( isTRUE(simulate_data) ){
          y_tj[t,j] = rbinom( n=1, size=1, prob=mu_tj[t,j] );
        }
      }
      devresid_tj[t,j] = sign(y_tj[t,j] - mu_tj[t,j]) * pow(-2*(((1-y_tj[t,j])*log(1-mu_tj[t,j]) + y_tj[t,j]*log(mu_tj[t,j]))), 0.5);
    }
    # familycode = 3 :  Poisson
    if( family[j]=="poisson" ){
      mu_tj[t,j] = exp(z_tj[t,j]);
      if(!is.na(y_tj[t,j])){
        loglik_tj[t,j] = dpois( y_tj[t,j], mu_tj[t,j], TRUE );
        if( isTRUE(simulate_data) ){
          y_tj[t,j] = rpois( n=1, lambda=mu_tj[t,j] );
        }
      }
      devresid_tj[t,j] = sign(y_tj[t,j] - mu_tj[t,j]) * pow(2*(y_tj[t,j]*log((Type(1e-10) + y_tj[t,j])/mu_tj[t,j]) - (y_tj[t,j]-mu_tj[t,j])), 0.5);
    }
    # familycode = 4 :  Gamma:   shape = 1/CV^2; scale = mean*CV^2
    if( family[j]=="gamma" ){
      mu_tj[t,j] = exp(z_tj[t,j]);
      if(!is.na(y_tj[t,j])){
        loglik_tj[t,j] = dgamma( y_tj[t,j], pow(sigma_j[j],-2), mu_tj[t,j]*pow(sigma_j[j],2), TRUE );
        if( isTRUE(simulate_data) ){
          y_tj[t,j] = rgamma( n=1, shape=pow(sigma_j[j],-2), scale=mu_tj[t,j]*pow(sigma_j[j],2) );
        }
      }
      devresid_tj[t,j] = sign(y_tj[t,j] - mu_tj[t,j]) * pow(2 * ( (y_tj[t,j]-mu_tj[t,j])/mu_tj[t,j] - log(y_tj[t,j]/mu_tj[t,j]) ), 0.5);
    }
  }}

  # Calculate priors
  log_prior_value = log_prior( parlist )

  jnll = -1 * sum(loglik_tj)
  jnll = jnll + jnll_gmrf - sum(log_prior_value)

  #
  REPORT( loglik_tj )
  REPORT( jnll_gmrf )
  REPORT( xhat_tj ) # needed to simulate new GMRF in R
  #REPORT( delta_k ) # FIXME>  Eliminate in simulate.dsem
  REPORT( delta_tj ) # needed to simulate new GMRF in R
  REPORT( Rho_kk )
  REPORT( Gamma_kk )
  REPORT( mu_tj )
  REPORT( devresid_tj )
  REPORT( IminusRho_kk )
  REPORT( V_kk )
  REPORT( jnll )
  REPORT( loglik_tj )
  REPORT( jnll_gmrf )
  REPORT( log_prior_value )
  #SIMULATE{
  #  REPORT( y_tj )
  #}
  REPORT( z_tj )
  ADREPORT( z_tj )

  if( isTRUE(simulate_data) ){
    list( y_tj=y_tj, z_tj=z_tj)
  }else{
    jnll
  }
}

if( FALSE ){
  family = rep("fixed",ncol(tsdata))
  y_tj = tsdata
  estimate_delta0 = FALSE
}

