
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
  Eigen::SparseMatrix<Type> Q_kk( n_k, n_k );
  // SEM
  Eigen::SparseMatrix<Type> Linv_kk(n_k, n_k);
  Eigen::SparseMatrix<Type> Rho_kk(n_k, n_k);
  Eigen::SparseMatrix<Type> Gamma_jj(n_j, n_j);
  Eigen::SparseMatrix<Type> Gamma_kk(n_k, n_k);
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
    //if((RAM(r,0)==2) && (RAM(r,1)==(n_t)) && (RAM(r,2)<=n_j)){
    //  Gamma_jj.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
    //}
    if(RAM(r,0)==2){
      Gamma_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = tmp; // Cholesky of covariance, so -Inf to Inf;
    }
    if(RAM(r,0)==2) Gammainv_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = 1 / tmp;
  }
  // Option-1
  //Eigen::SparseMatrix<Type> Q1_kk( n_k, n_k );
  //Linv_kk = Gammainv_kk * ( I_kk - Rho_kk );
  //Q1_kk = Linv_kk.transpose() * Linv_kk;

  //Eigen::SparseMatrix<Type> Q2_kk( n_k, n_k );
  //REPORT( Gamma_jj );
  //matrix<Type> Gammainv_jj( n_j, n_j );
  //Gammainv_jj = invertSparseMatrix( Gamma_jj );  // Returns dense matrix (trying to return sparse throws compiler error)
  //REPORT( Gammainv_jj );
  //Eigen::SparseMatrix<Type> Gammainv2_jj( n_j, n_j );
  //Gammainv2_jj = asSparseMatrix( Gammainv_jj );
  //REPORT( Gammainv2_jj );
  //Eigen::SparseMatrix<Type> I_tt( n_t, n_t );
  //I_tt.setIdentity();
  //REPORT( I_tt );
  //Eigen::SparseMatrix<Type> Gammainv2_kk( n_k, n_k );
  //Gammainv2_kk = kronecker( Gammainv2_jj, I_tt );
  //REPORT( Gammainv2_kk );
  //REPORT( Gammainv_kk );
  //Linv_kk = asSparseMatrix(Gammainv2_kk) * ( I_kk - Rho_kk );

  // Option-3
  Eigen::SparseMatrix<Type> V_kk( n_k, n_k );
  V_kk = Gamma_kk.transpose() * Gamma_kk;
  matrix<Type> Vinv_kk( n_k, n_k );
  Vinv_kk = invertSparseMatrix( V_kk );
  Eigen::SparseMatrix<Type> Vinv2_kk( n_k, n_k );
  Vinv2_kk = asSparseMatrix( Vinv_kk );
  REPORT( Gamma_kk );
  REPORT( Vinv2_kk );
  Eigen::SparseMatrix<Type> Linv2_kk(n_k, n_k);
  Linv2_kk = I_kk - Rho_kk;
  Q_kk = Linv2_kk.transpose() * Vinv2_kk * Linv2_kk;

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
  for(int t=0; t<n_t; t++){
  for(int j=0; j<n_j; j++){
    // familycode = 0 :  don't include likelihood
    if( familycode_j(j)==0 ){
      mu_tj(t,j) = x_tj(t,j);
      if(!R_IsNA(asDouble(y_tj(t,j)))){
        SIMULATE{
          y_tj(t,j) = mu_tj(t,j);
        }
      }
      devresid_tj(t,j) = 0;
    }
    // familycode = 1 :  normal
    if( familycode_j(j)==1 ){
      mu_tj(t,j) = x_tj(t,j);
      if(!R_IsNA(asDouble(y_tj(t,j)))){
        loglik_tj(t,j) = dnorm( y_tj(t,j), mu_tj(t,j), sigma_j(j), true );
        SIMULATE{
          y_tj(t,j) = rnorm( mu_tj(t,j), sigma_j(j) );
        }
      }
      devresid_tj(t,j) = y_tj(t,j) - mu_tj(t,j);
    }
    // familycode = 2 :  binomial
    if( familycode_j(j)==2 ){
      mu_tj(t,j) = invlogit(x_tj(t,j));
      if(!R_IsNA(asDouble(y_tj(t,j)))){
        loglik_tj(t,j) = dbinom( y_tj(t,j), Type(1.0), mu_tj(t,j), true );
        SIMULATE{
          y_tj(t,j) = rbinom( Type(1), mu_tj(t,j) );
        }
      }
      devresid_tj(t,j) = sign(y_tj(t,j) - mu_tj(t,j)) * pow(-2*(((1-y_tj(t,j))*log(1-mu_tj(t,j)) + y_tj(t,j)*log(mu_tj(t,j)))), 0.5);
    }
    // familycode = 3 :  Poisson
    if( familycode_j(j)==3 ){
      mu_tj(t,j) = exp(x_tj(t,j));
      if(!R_IsNA(asDouble(y_tj(t,j)))){
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
      if(!R_IsNA(asDouble(y_tj(t,j)))){
        loglik_tj(t,j) = dgamma( y_tj(t,j), pow(sigma_j(j),-2), mu_tj(t,j)*pow(sigma_j(j),2), true );
        SIMULATE{
          y_tj(t,j) = rgamma( pow(sigma_j(j),-2), mu_tj(t,j)*pow(sigma_j(j),2) );
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
