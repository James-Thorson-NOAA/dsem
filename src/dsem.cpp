
#define TMB_LIB_INIT R_init_dsem
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //using namespace Eigen;
  using namespace density;

  // Data
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

  // Distribution for data
  for(int t=0; t<n_t; t++){
  for(int j=0; j<n_j; j++){
    if( !R_IsNA(asDouble(y_tj(t,j))) ){
      // familycode = 0 :  don't include likelihood
      // familycode = 1 :  normal
      if( familycode_j(j)==1 ){
        loglik_tj(t,j) = dnorm( y_tj(t,j), x_tj(t,j), sigma_j(j), true );
      }
      // familycode = 2 :  binomial
      if( familycode_j(j)==2 ){
        loglik_tj(t,j) = dbinom( y_tj(t,j), Type(1.0), invlogit(x_tj(t,j)), true );
      }
      // familycode = 3 :  Poisson
      if( familycode_j(j)==3 ){
        loglik_tj(t,j) = dpois( y_tj(t,j), exp(x_tj(t,j)), true );
      }
      // familycode = 4 :  Gamma:   shape = 1/CV^2; scale = mean*CV^2
      if( familycode_j(j)==4 ){
        loglik_tj(t,j) = dgamma( y_tj(t,j), pow(sigma_j(j),-2), exp(x_tj(t,j))*pow(sigma_j(j),2), true );
      }
    }
  }}
  jnll -= loglik_tj.sum();
  jnll += jnll_gmrf;

  // Reporting
  //REPORT( V_kk );
  REPORT( Q_kk );
  REPORT( Rho_kk );
  //REPORT( Gammainv_kk );
  REPORT( jnll );
  REPORT( loglik_tj );
  REPORT( jnll_gmrf );
  return jnll;
}
