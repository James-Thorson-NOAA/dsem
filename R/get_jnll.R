# model function
get_jnll <-
function( parlist,
          model,
          tsdata,
          family ) {

  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

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
  n_t = nrow(tsdata)
  n_j = ncol(tsdata)
  n_k = prod(dim(tsdata))

  # Unpack
  #model_unique = model[match(unique(model$parameter),model$parameter),]
  #beta_p[which(model_unique$direction==1)] = beta_p1
  #beta_p[which(model_unique$direction==2)] = beta_p2

  # Build matrices
  ram = make_matrices(
            beta_p = beta_z,
            model = model,
            times = as.numeric(time(tsdata)),
            variables = colnames(tsdata) )

  # Assemble
  Rho_kk = AD(ram$P_kk)
  IminusRho_kk = Diagonal(n_k) - Rho_kk
  # Assemble variance
  Gamma_kk = invV_kk = AD(ram$G_kk)
  invV_kk@x = 1 / Gamma_kk@x^2

  # Calculate effect of initial condition -- SPARSE version
  delta_tj = matrix( 0, nrow=n_t, ncol=n_j )
  if( length(delta0_j)>0 ){
    delta_tj[1,] = delta0_j
    delta_k = solve( IminusRho_kk, as.vector(delta_tj) )
    delta_tj[] = delta_k
  }

  #
  xhat_tj = outer( rep(1,n_t), mu_j )
  z_tj = x_tj

  # Probability of GMRF
  Q_kk = t(IminusRho_kk) %*% invV_kk %*% IminusRho_kk
  jnll_gmrf = -1 * dgmrf( as.vector(z_tj), mu=as.vector(xhat_tj + delta_tj), Q=Q_kk, log=TRUE )
  REPORT( Q_kk )

  # Likelihood of data | random effects
  loglik_tj = mu_tj = devresid_tj = matrix( 0, nrow=n_t, ncol=n_j )
  pow = function(a,b) a^b
  for( t in 1:n_t ){
  for( j in 1:n_j ){
    # familycode = 0 :  don't include likelihood
    if( family[j]=="fixed" ){
      mu_tj[t,j] = z_tj[t,j];
      if(!is.na(y_tj[t,j])){
        #SIMULATE{
        #  y_tj[t,j] = mu_tj[t,j];
        #}
      }
      devresid_tj[t,j] = 0;
    }
    # familycode = 1 :  normal
    if( family[j]=="normal" ){
      mu_tj[t,j] = z_tj[t,j];
      if(!is.na(y_tj[t,j])){
        loglik_tj[t,j] = dnorm( y_tj[t,j], mu_tj[t,j], sigma_j[j], TRUE );
        #SIMULATE{
        #  y_tj[t,j] = rnorm( mu_tj[t,j], sigma_j[j] );
        #}
      }
      devresid_tj[t,j] = y_tj[t,j] - mu_tj[t,j];
    }
    # familycode = 2 :  binomial
    if( family[j]=="binomial" ){
      mu_tj[t,j] = invlogit(z_tj[t,j]);
      if(!is.na(y_tj[t,j])){
        loglik_tj[t,j] = dbinom( y_tj[t,j], Type(1.0), mu_tj[t,j], TRUE );
        #SIMULATE{
        #  y_tj[t,j] = rbinom( Type(1), mu_tj[t,j] );
        #}
      }
      devresid_tj[t,j] = sign(y_tj[t,j] - mu_tj[t,j]) * pow(-2*(((1-y_tj[t,j])*log(1-mu_tj[t,j]) + y_tj[t,j]*log(mu_tj[t,j]))), 0.5);
    }
    # familycode = 3 :  Poisson
    if( family[j]=="poisson" ){
      mu_tj[t,j] = exp(z_tj[t,j]);
      if(!is.na(y_tj[t,j])){
        loglik_tj[t,j] = dpois( y_tj[t,j], mu_tj[t,j], TRUE );
        #SIMULATE{
        #  y_tj[t,j] = rpois( mu_tj[t,j] );
        #}
      }
      devresid_tj[t,j] = sign(y_tj[t,j] - mu_tj[t,j]) * pow(2*(y_tj[t,j]*log((Type(1e-10) + y_tj[t,j])/mu_tj[t,j]) - (y_tj[t,j]-mu_tj[t,j])), 0.5);
    }
    # familycode = 4 :  Gamma:   shape = 1/CV^2; scale = mean*CV^2
    if( family[j]=="gamma" ){
      mu_tj[t,j] = exp(z_tj[t,j]);
      if(!is.na(y_tj[t,j])){
        loglik_tj[t,j] = dgamma( y_tj[t,j], pow(sigma_j[j],-2), mu_tj[t,j]*pow(sigma_j[j],2), TRUE );
        #SIMULATE{
        #  y_tj[t,j] = rgamma( pow(sigma_j[j],-2), mu_tj[t,j]*pow(sigma_j[j],2) );
        #}
      }
      devresid_tj[t,j] = sign(y_tj[t,j] - mu_tj[t,j]) * pow(2 * ( (y_tj[t,j]-mu_tj[t,j])/mu_tj[t,j] - log(y_tj[t,j]/mu_tj[t,j]) ), 0.5);
    }
  }}
  jnll = -1 * sum(loglik_tj);
  jnll = jnll + jnll_gmrf;

  #
  REPORT( loglik_tj )
  REPORT( jnll_gmrf )
  REPORT( xhat_tj ); # needed to simulate new GMRF in R
  #REPORT( delta_k ); # FIXME>  Eliminate in simulate.dsem
  REPORT( delta_tj ); # needed to simulate new GMRF in R
  REPORT( Rho_kk );
  REPORT( Gamma_kk );
  REPORT( mu_tj );
  REPORT( devresid_tj );
  REPORT( IminusRho_kk );
  REPORT( jnll );
  REPORT( loglik_tj );
  REPORT( jnll_gmrf );
  #SIMULATE{
  #  REPORT( y_tj );
  #}
  REPORT( z_tj );
  ADREPORT( z_tj );

  jnll
}

if( FALSE ){
  family = rep("fixed",ncol(tsdata))
  y_tj = tsdata
  estimate_delta0 = FALSE
}

