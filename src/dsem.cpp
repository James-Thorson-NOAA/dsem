
#define TMB_LIB_INIT R_init_dsem
#include <TMB.hpp>

// get sign of double, only for REPORT use
template<class Type>
Type sign(Type x){
  return x / pow(pow(x,2),0.5);
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  //using namespace Eigen;
  using namespace density;

  // Data
  //DATA_INTEGER( resimulate_gmrf );
  DATA_IMATRIX( RAM );
  DATA_VECTOR( RAMstart );
  DATA_IVECTOR( familycode_j );
  DATA_ARRAY( y_tj );

  // Parameters
  PARAMETER_VECTOR( beta_z );
  PARAMETER_VECTOR( lnsigma_j );
  PARAMETER_VECTOR( mu_j );
  PARAMETER_VECTOR( delta0_j );
  PARAMETER_ARRAY( x_tj );
  // vector of real observations (NA's removed) for OSA residuals, k2<=k
  DATA_VECTOR(y_k2); 		
  DATA_VECTOR_INDICATOR(keep, y_k2);

  // Indices
  int n_t = x_tj.rows();
  int n_j = x_tj.cols();
  int n_k = n_t * n_j;      // data

  // globals
  Type jnll = 0;
  Type jnll_gmrf = 0;
  matrix<Type> loglik_tj( n_t, n_j );
  loglik_tj.setZero();
  vector<Type> sigma_j( n_j );
  sigma_j = exp( lnsigma_j );

  // Assemble precision
  // Using Gamma_kk seems to crash when fixing values
  Eigen::SparseMatrix<Type> Q_kk( n_k, n_k );
  // SEM
  Eigen::SparseMatrix<Type> Linv_kk(n_k, n_k);
  Eigen::SparseMatrix<Type> Rho_kk(n_k, n_k);
  //Eigen::SparseMatrix<Type> Gamma_kk(n_k, n_k);
  Eigen::SparseMatrix<Type> Gammainv_kk(n_k, n_k);
  //matrix<Type> Gammainv2_kk(n_k, n_k);
  Eigen::SparseMatrix<Type> I_kk( n_k, n_k );
  Rho_kk.setZero();
  //Gamma_kk.setZero();
  Gammainv_kk.setZero();
  //Gammainv2_kk.setZero();
  I_kk.setIdentity();
  Type tmp;
  for(int r=0; r<RAM.rows(); r++){
    // Extract estimated or fixed value
    if(RAM(r,3)>=1){
      tmp = beta_z(RAM(r,3)-1);
    }else{
      tmp = RAMstart(r);
    }
    if(RAM(r,0)==1) Rho_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
    //if(RAM(r,0)==2) Gamma_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = beta_z(RAM(r,3)-1); // Cholesky of covariance, so -Inf to Inf;
    if(RAM(r,0)==2) Gammainv_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = 1 / tmp;
  }
  //Gammainv2_kk = invertSparseMatrix( Gamma_kk );
  //Linv_kk = asSparseMatrix(Gammainv2_kk) * ( I_kk - Rho_kk );
  Linv_kk = Gammainv_kk * ( I_kk - Rho_kk );
  Q_kk = Linv_kk.transpose() * Linv_kk;


  // Calculate effect of initial condition -- SPARSE version
  vector<Type> delta_k( n_k );
  delta_k.setZero();
  if( delta0_j.size() > 0 ){
    // Compute delta_k
    matrix<Type> delta0_k1( n_k, 1 );
    delta0_k1.setZero();
    int k;
    for(int j=0; j<n_j; j++){
      k = j * n_t;
      delta0_k1(k,0) = delta0_j(j);
    }

    // DENSE version
    //matrix<Type> Dense_kk = I_kk - Rho_kk;
    //matrix<Type> invIminusRho_kk = atomic::matinv( Dense_kk );
    //delta_k = (invIminusRho_kk * delta0_k1).array();

    // SPARSE version
    // See C:\Users\James.Thorson\Desktop\Work files\AFSC\2023-06 -- Sparse inverse-product\Kasper example\lu.cpp
    Eigen::SparseMatrix<Type> IminusRho_kk = I_kk - Rho_kk;
    Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > lu;
    lu.compute(IminusRho_kk);
    matrix<Type> x = lu.solve(delta0_k1);

    REPORT( delta0_k1 );
    delta_k = x.array();
  }
  REPORT( delta_k );

  // Centered GMRF
  array<Type> xhat_tj( n_t, n_j );
  for(int t=0; t<n_t; t++){
  for(int j=0; j<n_j; j++){
    xhat_tj(t,j) = mu_j(j);
  }}
  jnll_gmrf = GMRF(Q_kk)( x_tj - xhat_tj - delta_k );
  //SIMULATE{
  //  if( resimulate_gmrf >= 1 ){
  //    //x_tj = GMRF(Q_kk).simulate(x_tj);
  //    //x_tj += xhat_tj + delta_k;
  //  }
  //  REPORT( x_tj );
  //}

  // Distribution for data
  array<Type> devresid_tj( n_t, n_j );
  array<Type> mu_tj( n_t, n_j );
  bool is_not_na;

  int k2=0;			// counter for non NA observations for OSA
  for(int t=0; t<n_t; t++){
  for(int j=0; j<n_j; j++){
    is_not_na=!R_IsNA(asDouble(y_tj(t,j)));
    
    // familycode = 0 :  don't include likelihood
    if( familycode_j(j)==0 ){
      mu_tj(t,j) = x_tj(t,j);
      if(is_not_na){
        SIMULATE{
          y_tj(t,j) = mu_tj(t,j);
        }
      }
      devresid_tj(t,j) = 0;
    }
    // familycode = 1 :  normal
    if( familycode_j(j)==1 ){
      mu_tj(t,j) = x_tj(t,j);
      if(is_not_na){
        loglik_tj(t,j) = keep(k2)*dnorm( y_k2(k2), mu_tj(t,j), sigma_j(j), true );
	loglik_tj(t,j) += keep.cdf_lower(k2)*log(pnorm(y_k2(k2), mu_tj(t,j), sigma_j(j)));
	loglik_tj(t,j) += keep.cdf_upper(k2)*log(1.0-pnorm(y_k2(k2), mu_tj(t,j), sigma_j(j)));
	k2++;
        SIMULATE{
          y_tj(t,j) = rnorm( mu_tj(t,j), sigma_j(j) );
        }
      }
      devresid_tj(t,j) = y_tj(t,j) - mu_tj(t,j);
    }
    // familycode = 2 :  binomial
    if( familycode_j(j)==2 ){
      mu_tj(t,j) = invlogit(x_tj(t,j));
      if(is_not_na){
	loglik_tj(t,j) = keep(k2)*dbinom( y_k2(k2), Type(1.0), mu_tj(t,j), true );
	//loglik_tj(t,j) += keep.cdf_lower(k2)*log(pbinom(y_k2(k2), Type(1.0), mu_tj(t,j)));
	//loglik_tj(t,j) += keep.cdf_upper(k2)*log(1.0-pbinom(y_k2(k2), Type(1.0), mu_tj(t,j)));
	//   loglik_tj(t,j) = dbinom( y_tj(t,j), Type(1.0), mu_tj(t,j), true );
	k2++;
        SIMULATE{
          y_tj(t,j) = rbinom( Type(1), mu_tj(t,j) );
        }
      }
      devresid_tj(t,j) = sign(y_tj(t,j) - mu_tj(t,j)) * pow(-2*(((1-y_tj(t,j))*log(1-mu_tj(t,j)) + y_tj(t,j)*log(mu_tj(t,j)))), 0.5);
    }
    // familycode = 3 :  Poisson
    if( familycode_j(j)==3 ){
      mu_tj(t,j) = exp(x_tj(t,j));
      if(is_not_na){
	loglik_tj(t,j) = keep(k2)*dpois( y_k2(k2), mu_tj(t,j), true );
	loglik_tj(t,j) += keep.cdf_lower(k2)*log(ppois(y_k2(k2), mu_tj(t,j)));
	loglik_tj(t,j) += keep.cdf_upper(k2)*log(1.0-ppois(y_k2(k2), mu_tj(t,j)));
	k2++;
        loglik_tj(t,j) = dpois( y_tj(t,j), mu_tj(t,j), true );
        SIMULATE{
          y_tj(t,j) = rpois( mu_tj(t,j) );
        }
      }
      devresid_tj(t,j) = sign(y_tj(t,j) - mu_tj(t,j)) * pow(2*(y_tj(t,j)*log((Type(1e-10) + y_tj(t,j))/mu_tj(t,j)) - (y_tj(t,j)-mu_tj(t,j))), 0.5);
    }
    // familycode = 4 :  Gamma:   shape = 1/CV^2; scale = mean*CV^2
    if( familycode_j(j)==4 ){
      mu_tj(t,j) = exp(x_tj(t,j));
      if(is_not_na){
	Type tmp1=pow(sigma_j(j),-2);
	Type tmp2=mu_tj(t,j)*pow(sigma_j(j),2);
	loglik_tj(t,j) = keep(k2)*dgamma(y_k2(k2), tmp1, tmp2, true );
	loglik_tj(t,j) += keep.cdf_lower(k2)*log(pgamma(y_k2(k2), tmp1, tmp2));
	loglik_tj(t,j) += keep.cdf_upper(k2)*log(1.0-pgamma(y_k2(k2), tmp1, tmp2));
	k2++;
        //loglik_tj(t,j) = dgamma( y_tj(t,j), tmp1 , tmp2, true );
        SIMULATE{
          y_tj(t,j) = rgamma(tmp1, tmp2);
        }
      }
      devresid_tj(t,j) = sign(y_tj(t,j) - mu_tj(t,j)) * pow(2 * ( (y_tj(t,j)-mu_tj(t,j))/mu_tj(t,j) - log(y_tj(t,j)/mu_tj(t,j)) ), 0.5);
    }
  }}
  jnll -= loglik_tj.sum();
  jnll += jnll_gmrf;

  // Reporting
  //REPORT( V_kk );
  REPORT( Q_kk );
  REPORT( xhat_tj ); // needed to simulate new GMRF in R
  REPORT( delta_k ); // needed to simulate new GMRF in R
  REPORT( Rho_kk );
  REPORT( mu_tj );
  REPORT( devresid_tj );
  //REPORT( Gammainv_kk );
  REPORT( jnll );
  REPORT( loglik_tj );
  REPORT( jnll_gmrf );
  SIMULATE{
    REPORT( y_tj );
  }
  return jnll;
}
