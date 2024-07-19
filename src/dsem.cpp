
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
  DATA_IVECTOR( options ); 
  // options(0) -> 0: full rank;  1: rank-reduced GMRF
  // options(1) -> 0: constant conditional variance;  1: constant marginal variance
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
  int n_t = y_tj.rows();
  int n_j = y_tj.cols();
  int n_k = n_t * n_j;      // data
  int k = 0;

  // globals
  Type jnll = 0;
  Type jnll_gmrf = 0;
  matrix<Type> loglik_tj( n_t, n_j );
  loglik_tj.setZero();
  vector<Type> sigma_j( n_j );
  sigma_j = exp( lnsigma_j );

  // Assemble precision
  // SEM
  Eigen::SparseMatrix<Type> Rho_kk(n_k, n_k);
  //Eigen::SparseMatrix<Type> Gamma_jj(n_j, n_j);
  Eigen::SparseMatrix<Type> Gamma_kk(n_k, n_k);
  //Eigen::SparseMatrix<Type> Gammainv_kk(n_k, n_k);
  Eigen::SparseMatrix<Type> I_kk( n_k, n_k );
  Rho_kk.setZero();
  //Gammainv_kk.setZero();
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
    if(RAM(r,0)==2){
      Gamma_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = tmp; // Cholesky of covariance, so -Inf to Inf;
    }
    //if(RAM(r,0)==2) Gammainv_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = 1 / tmp;
  }
  Eigen::SparseMatrix<Type> IminusRho_kk = I_kk - Rho_kk;
  
  // Compute inverse LU-decomposition
  Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > inverseIminusRho_kk;
  inverseIminusRho_kk.compute(IminusRho_kk);

  // Rescale I-Rho and Gamma if using constant marginal variance options
  if( (options(1)==1) || (options(1)==2) ){
    Eigen::SparseMatrix<Type> invIminusRho_kk;
    
    // WORKS:  Based on: https://github.com/kaskr/adcomp/issues/74
    invIminusRho_kk = inverseIminusRho_kk.solve(I_kk);
    
    // Hadamard squared LU-decomposition
    // See: https://eigen.tuxfamily.org/dox/group__QuickRefPage.html
    Eigen::SparseMatrix<Type> squared_invIminusRho_kk(n_k, n_k);
    squared_invIminusRho_kk = invIminusRho_kk.cwiseProduct(invIminusRho_kk);
    Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > invsquared_invIminusRho_kk;
    invsquared_invIminusRho_kk.compute(squared_invIminusRho_kk);
    
    if( options(1) == 1 ){
      // 1-matrix
      matrix<Type> ones_k1( n_k, 1 );
      ones_k1.setOnes();

      // Calculate diag( t(Gamma) * Gamma )
      Eigen::SparseMatrix<Type> squared_Gamma_kk = Gamma_kk.cwiseProduct(Gamma_kk);
      matrix<Type> sigma2_k1 = squared_Gamma_kk.transpose() * ones_k1;

      // Rowsums
      matrix<Type> margvar_k1 = invsquared_invIminusRho_kk.solve(sigma2_k1);
      
      // Rescale IminusRho_kk and Gamma
      Eigen::SparseMatrix<Type> invmargsd_kk(n_k, n_k);
      Eigen::SparseMatrix<Type> invsigma_kk(n_k, n_k);
      for( int k=0; k<n_k; k++ ){
        invmargsd_kk.coeffRef(k,k) = pow( margvar_k1(k,0), -0.5 );
        invsigma_kk.coeffRef(k,k) = pow( sigma2_k1(k,0), -0.5 );
      }
      IminusRho_kk = invmargsd_kk * IminusRho_kk;
      Gamma_kk = invsigma_kk * Gamma_kk;
      
      // Recompute inverse LU-decomposition
      inverseIminusRho_kk.compute(IminusRho_kk);
    }else{
      // calculate diag(Gamma)^2
      matrix<Type> targetvar_k1( n_k, 1 );
      for( int k=0; k<n_k; k++ ){
        targetvar_k1(k,0) = Gamma_kk.coeffRef(k,k) * Gamma_kk.coeffRef(k,k);
      }
      
      // Rescale Gamma
      matrix<Type> margvar_k1 = invsquared_invIminusRho_kk.solve(targetvar_k1);
      for( int k=0; k<n_k; k++ ){
        Gamma_kk.coeffRef(k,k) = pow( margvar_k1(k,0), 0.5 );
      }
    }
  }
  
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
    matrix<Type> x = inverseIminusRho_kk.solve(delta0_k1);

    REPORT( delta0_k1 );
    delta_k = x.array();
  }
  REPORT( delta_k );

  // Format mu_j
  array<Type> xhat_tj( n_t, n_j );
  array<Type> delta_tj( n_t, n_j );
  for(int j=0; j<n_j; j++){
  for(int t=0; t<n_t; t++){
    k = j*n_t + t;
    xhat_tj(t,j) = mu_j(j);
    delta_tj(t,j) = delta_k(k);
  }}

  // Apply GMRF
  array<Type> z_tj( n_t, n_j );
  if( options(0)==0 ){
    // Only compute Vinv_kk if Gamma_kk is full rank
    Eigen::SparseMatrix<Type> V_kk = Gamma_kk.transpose() * Gamma_kk;
    matrix<Type> Vinv_kk = invertSparseMatrix( V_kk );
    Eigen::SparseMatrix<Type> Vinv2_kk = asSparseMatrix( Vinv_kk );
    Eigen::SparseMatrix<Type> Q_kk = IminusRho_kk.transpose() * Vinv2_kk * IminusRho_kk;
    
    // Centered GMRF
    jnll_gmrf = GMRF(Q_kk)( x_tj - xhat_tj - delta_tj );
    z_tj = x_tj;
    REPORT( Q_kk );
  }else{
    // Rank-deficient (projection) method
    jnll_gmrf += GMRF(I_kk)( x_tj );

    // Forward-format matrix
    matrix<Type> z_k1( n_t*n_j, int(1) );
    for(int j=0; j<n_j; j++){
    for(int t=0; t<n_t; t++){
      k = j*n_t + t;
      z_k1(k,0) = x_tj(t,j);
    }}

    // (I-Rho)^{-1} * Gamma * Epsilon
    matrix<Type> z2_k1 = Gamma_kk * z_k1;
    matrix<Type> z3_k1 = inverseIminusRho_kk.solve(z2_k1);

    // Back-format vector
    for(int j=0; j<n_j; j++){
    for(int t=0; t<n_t; t++){
      k = j*n_t + t;
      z_tj(t,j) = z3_k1(k,0);
    }}
    
    // Add back mean and deviation
    z_tj += xhat_tj + delta_tj;
  }
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
      mu_tj(t,j) = z_tj(t,j);
      if(!R_IsNA(asDouble(y_tj(t,j)))){
        SIMULATE{
          y_tj(t,j) = mu_tj(t,j);
        }
      }
      devresid_tj(t,j) = 0;
    }
    // familycode = 1 :  normal
    if( familycode_j(j)==1 ){
      mu_tj(t,j) = z_tj(t,j);
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
      mu_tj(t,j) = invlogit(z_tj(t,j));
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
      mu_tj(t,j) = exp(z_tj(t,j));
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
      mu_tj(t,j) = exp(z_tj(t,j));
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
  REPORT( xhat_tj ); // needed to simulate new GMRF in R
  REPORT( delta_k ); // FIXME>  Eliminate in simulate.dsem
  REPORT( delta_tj ); // needed to simulate new GMRF in R
  REPORT( Rho_kk );
  REPORT( Gamma_kk );
  REPORT( mu_tj );
  REPORT( devresid_tj );
  REPORT( IminusRho_kk );
  //REPORT( Gammainv_kk );
  REPORT( jnll );
  REPORT( loglik_tj );
  REPORT( jnll_gmrf );
  SIMULATE{
    REPORT( y_tj );
  }
  REPORT( z_tj );
  ADREPORT( z_tj );
  return jnll;
}
