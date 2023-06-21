
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// SparseMatrix for Ornstein-Uhlenbeck network correlations
//template<class Type>
//Eigen::SparseMatrix<Type> Q_network( Type log_alpha,
//                                     int n_s,
//                                     vector<int> parent_s,
//                                     vector<int> child_s,
//                                     vector<Type> dist_s ){
//
//  Eigen::SparseMatrix<Type> Q( n_s, n_s );
//  Type alpha = exp( log_alpha );
//  for(int s=0; s<n_s; s++){
//    Q.coeffRef( s, s ) = Type(1.0);
//  }
//  for(int s=1; s<parent_s.size(); s++){
//    if( exp(-dist_s(s))!=0 ){
//      Q.coeffRef( parent_s(s), child_s(s) ) = -exp(-alpha*dist_s(s)) / (1-exp(-2*alpha*dist_s(s)));
//      Q.coeffRef( child_s(s), parent_s(s) ) = Q.coeffRef( parent_s(s), child_s(s) );
//      Q.coeffRef( parent_s(s), parent_s(s) ) += exp(-2*alpha*dist_s(s)) / (1-exp(-2*alpha*dist_s(s)));
//      Q.coeffRef( child_s(s), child_s(s) ) += exp(-2*alpha*dist_s(s)) / (1-exp(-2*alpha*dist_s(s)));
//    }
//  }
//  return Q;
//}

// ICAR for BM with sum-to-zero constraint
//  See: C:\Users\James.Thorson\Desktop\Work files\AFSC\2022-09 -- ICAR specification

// Precision of evolutionary covariance
//template<class Type>
//Eigen::SparseMatrix<Type> Q_sem( vector<Type> beta_z,
//                                 matrix<int> RAM,
//                                 int n_vars ){
//
//  // Define temporary objects
//  Eigen::SparseMatrix<Type> Q_vv( n_vars, n_vars );
//  // SEM
//  Eigen::SparseMatrix<Type> Linv_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> Rho_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> Gamma_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> Gammainv_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> I_vv( n_vars, n_vars );
//  Rho_vv.setZero();
//  Gamma_vv.setZero();
//  I_vv.setIdentity();
//  for(int zI=0; zI<RAM.rows(); zI++){
//    if(RAM(zI,0)==1) Rho_vv.coeffRef( RAM(zI,1)-1, RAM(zI,2)-1 ) = beta_z(RAM(zI,3)-1);
//    if(RAM(zI,0)==2) Gamma_vv.coeffRef( RAM(zI,1)-1, RAM(zI,2)-1 ) = beta_z(RAM(zI,3)-1); // Cholesky of covariance, so -Inf to Inf;
//  }
//  Gammainv_vv = atomic::matinv( Gamma_vv );
//  Linv_vv = Gammainv_vv * ( I_vv - Rho_vv );
//  Q_vv = Linv_vv.transpose() * Linv_vv;
//  return Q_vv;
//}

// Evolutionary covariance
//template<class Type>
//matrix<Type> V_sem( vector<Type> beta_z,
//                                 matrix<int> RAM,
//                                 vector<Type> RAMstart,
//                                 int n_vars ){
//
//  // Define temporary objects
//  matrix<Type> V_vv(n_vars, n_vars);
//  // SEM
//  matrix<Type> L_vv(n_vars, n_vars);
//  matrix<Type> Rho_vv(n_vars, n_vars);
//  matrix<Type> Gamma_vv(n_vars, n_vars);
//  matrix<Type> I_vv( n_vars, n_vars );
//  Rho_vv.setZero();
//  Gamma_vv.setZero();
//  I_vv.setIdentity();
//  Type tmp;
//  for(int r=0; r<RAM.rows(); r++){
//    // Extract estimated or fixed value
//    if(RAM(r,3)>=1){
//      tmp = beta_z(RAM(r,3)-1);
//    }else{
//      tmp = RAMstart(r);
//    }
//    // Assign to proper matrix
//    if(RAM(r,0)==1){
//      Rho_vv( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
//    }else{
//      Gamma_vv( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
//    }
//  }
//  L_vv = I_vv - Rho_vv;
//  L_vv = atomic::matinv( L_vv );
//  L_vv = L_vv * Gamma_vv;
//  V_vv = L_vv * L_vv.transpose();
//  return V_vv;
//}

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

//  // Assemble covariance
//  // SEM
//  matrix<Type> V_kk( n_k, n_k );
//  matrix<Type> L_kk(n_k, n_k);
//  matrix<Type> Rho_kk(n_k, n_k);
//  matrix<Type> Gamma_kk(n_k, n_k);
//  matrix<Type> I_kk( n_k, n_k );
//  Rho_kk.setZero();
//  Gamma_kk.setZero();
//  I_kk.setIdentity();
//  Type tmp;
//  for(int r=0; r<RAM.rows(); r++){
//    // Extract estimated or fixed value
//    //if(RAM(r,3)>=1){
//      tmp = beta_z(RAM(r,3)-1);
//    //}//else{
//      // tmp = RAMstart(r);
//    //}
//    // Assign to proper matrix
//    if(RAM(r,0)==1){
//      Rho_kk( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
//    }else{
//      Gamma_kk( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
//    }
//  }
//  L_kk = I_kk - Rho_kk;
//  L_kk = atomic::matinv( L_kk );
//  L_kk = L_kk * Gamma_kk;
//  V_kk = L_kk * L_kk.transpose();
//  jnll += MVNORM(V_kk)( x_j );

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

  // Centered GMRF
  array<Type> xhat_tj( n_t, n_j );
  for(int t=0; t<n_t; t++){
  for(int j=0; j<n_j; j++){
    xhat_tj(t,j) = mu_j(j);
  }}
  jnll_gmrf = GMRF(Q_kk)( x_tj - xhat_tj );

  // Distribution for data
  for(int t=0; t<n_t; t++){
  for(int j=0; j<n_j; j++){
    if( !isNA(y_tj(t,j)) ){
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
