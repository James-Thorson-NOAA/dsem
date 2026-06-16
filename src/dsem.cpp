
#define TMB_LIB_INIT R_init_dsem
#include <TMB.hpp>

// get sign of double, only for REPORT use
template<class Type>
Type sign(Type x){
  return x / pow(pow(x,2),0.5);
}

// Deviance for the Tweedie
// https://en.wikipedia.org/wiki/Tweedie_distribution#Properties
template<class Type>
Type devresid_tweedie( Type y,
                       Type mu,
                       Type p ){

  Type c1 = pow( y, 2.0-p ) / (1.0-p) / (2.0-p);
  Type c2 = y * pow( mu, 1.0-p ) / (1.0-p);
  Type c3 = pow( mu, 2.0-p ) / (2.0-p);
  Type deviance = 2 * (c1 - c2 + c3 );
  Type devresid = sign( y - mu ) * pow( deviance, 0.5 );
  return devresid;
}

// dlnorm
template<class Type>
Type dlnorm( Type x,
             Type meanlog,
             Type sdlog,
             int give_log=0){

  //return 1/(sqrt(2*M_PI)*sd) * exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// Get sparse submatrix, for use in dgmrf_conditional
// Modified from chatGPT-5
template<class Type>
Eigen::SparseMatrix<Type> get_submatrix( Eigen::SparseMatrix<Type> A,
                                            vector<int> row_idx,
                                            vector<int> col_idx ){

  // Build submatrix manually
  Eigen::SparseMatrix<Type> sub(row_idx.size(), col_idx.size());

  for (int k = 0; k < A.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<Type>::InnerIterator it(A, k); it; ++it) {
      // find if row and col are in the selection
      auto row_pos = std::find(row_idx.data(), row_idx.data() + row_idx.size(), it.row());
      auto col_pos = std::find(col_idx.data(), col_idx.data() + col_idx.size(), it.col());
      if (row_pos != row_idx.data() + row_idx.size() && col_pos != col_idx.data() + col_idx.size()) {
        int new_row = row_pos - row_idx.data();
        int new_col = col_pos - col_idx.data();
        sub.coeffRef(new_row, new_col) = it.value();
      }
    }
  }
  return sub;
}

// Evaluate negative log-density from conditional-GMRF
// see scratch/simulate_conditional_gmrf.R
// modified from tinyVAST::conditional_gmrf
//template<class Type>
//Type GMRF_conditional( vector<Type> x,  // x[obs_idx] is the observed values
//                       Eigen::SparseMatrix<Type> Q,
//                       vector<int> obs_idx,
//                       vector<int> unobs_idx ){
//  using namespace density;
//  Type out;
//
//  // Only compute if some are conditional
//  if( obs_idx.size() > 0 ){
//    // Partition Q
//    Eigen::SparseMatrix<Type> Q_uo = get_submatrix( Q, unobs_idx, obs_idx );
//    Eigen::SparseMatrix<Type> Q_uu = get_submatrix( Q, unobs_idx, unobs_idx );
//
//    // Extract observed values
//    vector<Type> obs_x( obs_idx.size() );
//    for (int i = 0; i < obs_idx.size(); i++) {
//      obs_x(i) = x(obs_idx(i));
//    }
//    vector<Type> unobs_x( unobs_idx.size() );
//    for (int i = 0; i < unobs_idx.size(); i++) {
//      unobs_x(i) = x(unobs_idx(i));
//    }
//
//    // sparseLU
//    Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > inverseQ_uu;
//    inverseQ_uu.compute(Q_uu);
//
//    // Compute conditional mean and covariance
//    // mu_cond <- -Q_uu_inv %*% Q_uo %*% x_obs
//    matrix<Type> projx = Q_uo * obs_x;
//    matrix<Type> mu_cond = -1 * inverseQ_uu.solve(projx);
//
//    //
//    vector<Type> diff_x = unobs_x - mu_cond.array();
//    out = GMRF( Q_uu )( diff_x );
//  }else{
//    out = GMRF( Q )( x );
//  }
//  return out;
//}

template<class Type>
Type objective_function<Type>::operator() ()
{
  //using namespace Eigen;
  using namespace density;

  // Data
  DATA_IVECTOR( options ); 
  // options(0) -> 0: full rank;  1: rank-reduced GMRF;  2: conditional krigging
  // options(1) -> 0: constant conditional variance;  1: constant marginal variance
  // options(2) -> 0: use GMRF(Q);  1: use GMRF(Q + 1e-14 \times I)
  // options(3) -> 0: natural-scale moderating variance;  1: log-scale moderating variance
  DATA_IMATRIX( RAM );
  DATA_VECTOR( RAMstart );
  DATA_IVECTOR( familycode_j );
  DATA_IVECTOR( linkcode_j );
  DATA_IVECTOR( sigmastart_j );
  DATA_ARRAY( eps_tj );
  DATA_ARRAY( y_tj );

  // Parameters
  PARAMETER_VECTOR( beta_z );
  PARAMETER_VECTOR( lnsigma_z );
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
  vector<Type> sigma_z = exp( lnsigma_z );

  // Assemble precision
  // SEM
  Eigen::SparseMatrix<Type> Rho_kk(n_k, n_k);
  Eigen::SparseMatrix<Type> Gamma_kk(n_k, n_k);
  Eigen::SparseMatrix<Type> I_kk( n_k, n_k );
  Rho_kk.setZero();
  I_kk.setIdentity();
  Type tmp;
  // RAM uses R-style indices (i.e., starting at 1)
  for(int r=0; r<RAM.rows(); r++){
    // Extract estimated or fixed value
    if(RAM(r,3)>=1){
      tmp = beta_z(RAM(r,3)-1);
    }else{
      tmp = RAMstart(r);
    }
    if(RAM(r,0) == 1){
      Rho_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
    }
    if(RAM(r,0) == 2){
      Gamma_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = tmp; // Cholesky of covariance, so -Inf to Inf;
    }
    if(RAM(r,0) == 3){
      Rho_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = x_tj( RAM(r,4)-1, RAM(r,5)-1 );
    }
    if(RAM(r,0) == 4){
      // Trying to decide whether to use log or natural-space
      if( options(3) == 0) Gamma_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = x_tj( RAM(r,4)-1, RAM(r,5)-1 );
      if( options(3) == 1) Gamma_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = exp(x_tj( RAM(r,4)-1, RAM(r,5)-1 ));
    }
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
    for(int j=0; j<n_j; j++){
      k = j * n_t;
      delta0_k1(k,0) = delta0_j(j);
    }

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
  //matrix<Type> tmp_tj;
  //vector<Type> ones_t( n_t );
  //ones_t.setOnes();
  //tmp_tj = (ones_t * mu_j.transpose());
  //xhat_tj = tmp_tj.array();
  //tmp_tj = delta_k.reshaped( n_t, n_j );
  //delta_tj = tmp_tj.array();

  // Apply GMRF
  array<Type> z_tj( n_t, n_j );
  // Option-1:  use full-rank GMRF
  if( options(0)==0 ){
    // Only compute Vinv_kk if Gamma_kk is full rank
    Eigen::SparseMatrix<Type> V_kk = Gamma_kk.transpose() * Gamma_kk;
    // Add diagonal if isTRUE(control$stabilize_Q)
    if( options(2) == 1 ){
      V_kk += I_kk * 1e-10;
    }
    matrix<Type> Vinv_kk = invertSparseMatrix( V_kk );
    Eigen::SparseMatrix<Type> Vinv2_kk = asSparseMatrix( Vinv_kk );
    Eigen::SparseMatrix<Type> Q_kk = IminusRho_kk.transpose() * Vinv2_kk * IminusRho_kk;
    
    // Eigen::SimplicialLDLT not working for some reason ... 
    //Eigen::SimplicialLDLT< Eigen::SparseMatrix<Type> > inverseV_kk;
    //inverseV_kk.compute(V_kk);
    //Eigen::SparseMatrix<Type> Q_kk = IminusRho_kk.transpose() * inverseV_kk.solve(IminusRho_kk);
    
    // Centered GMRF
    jnll_gmrf = GMRF(Q_kk)( x_tj - xhat_tj - delta_tj );
    z_tj = x_tj;
    REPORT( Q_kk );
  }
  // Option-2:  Rank-deficient (projection) method
  if( options(0)==1 ){
    jnll_gmrf = GMRF(I_kk)( x_tj );

    // Forward-format matrix
    matrix<Type> z_k1 = x_tj.reshaped( n_k, 1 );

    // (I-Rho)^{-1} * Gamma * Epsilon
    matrix<Type> z2_k1 = Gamma_kk * z_k1;
    matrix<Type> z3_k1 = inverseIminusRho_kk.solve(z2_k1);

    // Back-format vector
    z_tj = z3_k1.reshaped( n_t, n_j );
    
    // Add back mean and deviation
    z_tj += xhat_tj + delta_tj;
  }
  // Option-3:  use variance for full-rank component and projects to reduced-rank component
  // ALlows family = "fixed" for some variables, and rank-deficiency for other variables (e.g., determinstic composite variables)
  if( options(0)==2 ){
    DATA_IVECTOR( obs_idx );    // Full-rank component
    DATA_IVECTOR( unobs_idx );  // Reduced-rank component ... projecting from obs_idx to unobs_idx
    //error("not implemented yet");
    Eigen::SparseMatrix<Type> I_uu( unobs_idx.size(), unobs_idx.size() );
    I_uu.setIdentity();

    // Compute full covariance (potentially rank deficient)
    Eigen::SparseMatrix<Type> V_kk = Gamma_kk.transpose() * Gamma_kk;
    Eigen::SparseMatrix<Type> tmp_kk = inverseIminusRho_kk.solve(V_kk);
    Eigen::SparseMatrix<Type> tmp2_kk = tmp_kk.transpose();
    
    // Get Sigma components
    Eigen::SparseMatrix<Type> Sigma_kk = inverseIminusRho_kk.solve( tmp2_kk );
    Eigen::SparseMatrix<Type> Sigma_oo = get_submatrix( Sigma_kk, obs_idx, obs_idx );
    matrix<Type> V_oo = matrix<Type>(Sigma_oo);
    REPORT( Sigma_kk );

    // Extract sub-vectors for observed and unobserved components
    vector<Type> x_k = x_tj;
    vector<Type> dev_k = x_tj - xhat_tj - delta_tj;
    vector<Type> dev_o( obs_idx.size() );
    Eigen::SparseMatrix<Type> x_o1( obs_idx.size(), 1 );
    for( int index = 0; index < obs_idx.size(); index++ ){
      dev_o(index) = dev_k( obs_idx(index) );
      x_o1.coeffRef(index, 0) = x_k( obs_idx(index) );
    }
    vector<Type> x_u( unobs_idx.size() );
    Eigen::SparseMatrix<Type> x_u1( unobs_idx.size(), 1 );
    for( int u = 0; u < unobs_idx.size(); u++ ){
      x_u(u) = x_k( unobs_idx(u) );
      x_u1.coeffRef(u, 0) = x_k( unobs_idx(u) );
    }

    // Project residuals
    // mu_u = (V_uo %*% solve(V_oo) %*% x_o)[,1]
    // Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > inverseSigma_oo;
    // Using Eigen::SimplicialLDLT instead of Eigen::SparseLU because it's symmetric
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<Type> > inverseSigma_oo;
    inverseSigma_oo.compute(Sigma_oo);
    Eigen::SparseMatrix<Type> tmp_o1 = inverseSigma_oo.solve(x_o1);
    Eigen::SparseMatrix<Type> Sigma_uo = get_submatrix( Sigma_kk, unobs_idx, obs_idx );
    matrix<Type> mu_u1 = Sigma_uo * tmp_o1;
    
    // Get variance and Cholesky for remaining terms
    matrix<Type> xprime_u1( unobs_idx.size(), 1 );
    if( unobs_idx.size() > 0 ){
      Eigen::SparseMatrix<Type> Vprime_uu( unobs_idx.size(), unobs_idx.size() );
      Eigen::SparseMatrix<Type> Sigma_ou = Sigma_uo.transpose();
      Eigen::SparseMatrix<Type> tmp_ou = inverseSigma_oo.solve(Sigma_ou);
      Eigen::SparseMatrix<Type> Sigma_uu = get_submatrix( Sigma_kk, unobs_idx, unobs_idx );
      
      // CRASHING
      //Vprime_uu = Sigma_uu - (Sigma_uo * tmp_ou);   // CRASHES
      Vprime_uu = Sigma_uu;                // FINE
      //Vprime_uu = Sigma_uo * tmp_ou;     // FINE
      Vprime_uu -= (Sigma_uo * tmp_ou); 
      Vprime_uu += 1e-12 * I_uu;      // 1e-16 crashes
      
      // CONTINUE
      Eigen::SimplicialLLT< SparseMatrix<Type> > chol(Vprime_uu);
      SparseMatrix<Type> Lprime_uu = chol.matrixL();
      xprime_u1 = Lprime_uu * x_u1;
      REPORT( Lprime_uu );
    }
    
    // Add projected residuals + other comonents into linear predictor
    z_tj = x_tj;
    int u = 0;
    if( unobs_idx.size() > 0 ){
      for(int j=0; j<n_j; j++){
      for(int t=0; t<n_t; t++){
        k = j*n_t + t;
        if( (u < unobs_idx.size()) && (unobs_idx(u)==k) ){
        //if( (unobs_idx(u)==k) ){
          z_tj(t,j) = mu_u1(u,0) + xhat_tj(t,j) + delta_tj(t,j) + xprime_u1(u,0);
          u++;
        }
      }}
    }
    
    // Evaluate MVN density for full-rank component
    jnll_gmrf = MVNORM(V_oo)( dev_o );
    jnll_gmrf += GMRF(I_uu)( x_u );
  }

  // Option-4:  use full rank (some of which are fixed), 
  //            and project to zero-rank component (none of which are fixed and measured)
  // 
  // Given
  // x = (x_o, x_u)^T
  // where
  // x_o has V_oo that is full rank, and some are fixed ("observed")
  // x_u has V_uu that has no rank (V_uu = 0), and none are fixed ("unobserved" and projected to)
  //
  // NOTE:  x_u cannot include moderator variables
  //
  // Define
  // P = | P_oo, P_ou |
  //     | P_uo, P_uu  |
  //
  // V = | V_oo, V_ou |
  //     | V_uo, V_uu  |
  //
  // M = | I-P_oo,  P_ou   |  =  | M_oo,  M_ou |
  //     | P_uu_oo, I-P_uu |     | M_uo, M_uu  |
  //
  // Calculate
  // C = M_ou M_uu^-1
  // so
  // C^T = (M_uu^T)^-1 M_ou^T
  // and
  // Mtilda_oo = M_oo - M_ou M_uu^-1 M_uo
  // Vtilda_oo = V_oo + C V_uu C^T + C V_uo + V_ou C^T 
  //           = V_oo + C V_uo + V_ou C^T   (because V_uu = 0)
  // Q_oo = Mtilda_oo^T Vtilda_oo^-1 Mtilda_oo
  //
  // Then:
  // x_o ~ GMRF( Q_oo )
  // mu_u = -M_uu^-1 M_uo x_A (conditional krigging)
  //
  // And 
  // x_u = mu_u
  // Because 
  // x_u ~ MVN( mu_u, Q_uu^-1 )
  // And:
  // Q_uu = M_uu^T V_uu^-1 M_uu 
  // so 
  // Q_uu^-1 = 0 (because V_uu = 0)
  if( options(0)==3 ){
    DATA_IVECTOR( obs_idx );    // Full-rank component
    DATA_IVECTOR( unobs_idx );  // Zero-rank component ... projecting from obs_idx to unobs_idx
    Eigen::SparseMatrix<Type> Vtilda_oo;
    Eigen::SparseMatrix<Type> Mtilda_oo;
    vector<Type> dev_o( obs_idx.size() );
    z_tj = x_tj;
    if( unobs_idx.size() > 0 ){
      // Extract sub-vectors for observed and unobserved components
      vector<Type> dev_k = x_tj - xhat_tj - delta_tj;
      for( int o = 0; o < obs_idx.size(); o++ ){
        dev_o(o) = dev_k( obs_idx(o) );
      }
      // Extract V components
      Eigen::SparseMatrix<Type> V_kk = Gamma_kk.transpose() * Gamma_kk;
      Eigen::SparseMatrix<Type> V_oo = get_submatrix( V_kk, obs_idx, obs_idx );
      Eigen::SparseMatrix<Type> V_uo = get_submatrix( V_kk, unobs_idx, obs_idx );
      Eigen::SparseMatrix<Type> V_ou = get_submatrix( V_kk, obs_idx, unobs_idx );
      // Extract M components
      Eigen::SparseMatrix<Type> M_oo = get_submatrix( IminusRho_kk, obs_idx, obs_idx );
      Eigen::SparseMatrix<Type> M_uo = get_submatrix( IminusRho_kk, unobs_idx, obs_idx );
      Eigen::SparseMatrix<Type> M_ou = get_submatrix( IminusRho_kk, obs_idx, unobs_idx );
      Eigen::SparseMatrix<Type> M_uu = get_submatrix( IminusRho_kk, unobs_idx, unobs_idx );
      // Compute C
      Eigen::SparseMatrix<Type> Mt_ou = M_ou.transpose();
      Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > inverseMt_uu;
      inverseMt_uu.compute( M_uu.transpose().eval() );
      Eigen::SparseMatrix<Type> Ct = inverseMt_uu.solve(Mt_ou);
      // Mtilda_oo
      Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > inverseM_uu;
      inverseM_uu.compute(M_uu);
      Mtilda_oo = M_oo - M_ou * inverseM_uu.solve(M_uo);
      // Vtilda_oo
      Vtilda_oo = V_oo + Ct.transpose()*V_uo + V_ou*Ct;
      // Calculate devs
      matrix<Type> dev_u1 = -(inverseM_uu.solve(M_uo) * dev_o.matrix());
      REPORT( dev_u1 );
      // Add projected residuals + other comonents into linear predictor
      int u = 0;
      for(int j=0; j<n_j; j++){
      for(int t=0; t<n_t; t++){
        k = j*n_t + t;
        if( (u < unobs_idx.size()) && (unobs_idx(u)==k) ){
        //if( (unobs_idx(u)==k) ){
          z_tj(t,j) = dev_u1(u,0) + xhat_tj(t,j) + delta_tj(t,j);
          u++;
        }
      }}
    }else{
      dev_o = x_tj - xhat_tj - delta_tj;
      Vtilda_oo = Gamma_kk.transpose() * Gamma_kk;
      // Add diagonal if isTRUE(control$stabilize_Q)
      if( options(2) == 1 ){
        Vtilda_oo += I_kk * 1e-10;
      }
      Mtilda_oo = IminusRho_kk;
    }

    // Q_oo:  Eigen::SimplicialLDLT instead of Eigen::SparseLU because it's symmetric
    // SEEMS UNSTABLE
    //Eigen::SimplicialLDLT< Eigen::SparseMatrix<Type> > inverseVtilda_oo;
    //inverseVtilda_oo.compute(Vtilda_oo);
    //Eigen::SparseMatrix<Type> Q_oo = Mtilda_oo.transpose() * inverseVtilda_oo.solve(Mtilda_oo);

    // Same way as option(0) = 0
    matrix<Type> inverseVtilda_oo = invertSparseMatrix( Vtilda_oo );
    Eigen::SparseMatrix<Type> inverseVtilda2_oo = asSparseMatrix( inverseVtilda_oo );
    Eigen::SparseMatrix<Type> Q_oo = Mtilda_oo.transpose() * inverseVtilda2_oo * Mtilda_oo;

    // Get GMRF for data
    REPORT( Q_oo );
    //REPORT( dev_o );
    jnll_gmrf = GMRF( Q_oo )( dev_o );   
  }

  // Distribution for data
  // Simulates new data even for NA values, which can then be excluded during simulate.dsem
  array<Type> devresid_tj( n_t, n_j );
  array<Type> mu_tj( n_t, n_j );
  for(int t=0; t<n_t; t++){
  for(int j=0; j<n_j; j++){
    // Link function
    if( linkcode_j(j)==0 ){
      // identity link
      mu_tj(t,j) = z_tj(t,j);
    }
    if( linkcode_j(j)==1 ){
      // log link
      mu_tj(t,j) = exp(z_tj(t,j));
    }
    if( linkcode_j(j)==2 ){
      // logit link
      mu_tj(t,j) = invlogit(z_tj(t,j));
    }
    if( linkcode_j(j)==3 ){
      // cloglog link
      mu_tj(t,j) = Type(1.0) - exp( -1.0 * exp(z_tj(t,j)) );
    }

    // Likelihood
    if( familycode_j(j)==0 ){
      // familycode = 0 :  don't include likelihood
      SIMULATE{
        y_tj(t,j) = mu_tj(t,j);
      }
      devresid_tj(t,j) = 0;
    }
    if( familycode_j(j)==1 ){
      // familycode = 1 :  normal
      if(R_FINITE(asDouble(y_tj(t,j)))){
        loglik_tj(t,j) = dnorm( y_tj(t,j), mu_tj(t,j), sigma_z(sigmastart_j(j)), true );
      }
      SIMULATE{
        y_tj(t,j) = rnorm( mu_tj(t,j), sigma_z(sigmastart_j(j)) );
      }
      devresid_tj(t,j) = y_tj(t,j) - mu_tj(t,j);
    }
    if( familycode_j(j)==2 ){
      // familycode = 2 :  Bernoulli
      if(R_FINITE(asDouble(y_tj(t,j)))){
        loglik_tj(t,j) = dbinom( y_tj(t,j), Type(1.0), mu_tj(t,j), true );
      }
      SIMULATE{
        y_tj(t,j) = rbinom( Type(1), mu_tj(t,j) );
      }
      devresid_tj(t,j) = sign(y_tj(t,j) - mu_tj(t,j)) * pow(-2*(((1-y_tj(t,j))*log(1-mu_tj(t,j)) + y_tj(t,j)*log(mu_tj(t,j)))), 0.5);
    }
    if( familycode_j(j)==3 ){
      // familycode = 3 :  Poisson
      if(R_FINITE(asDouble(y_tj(t,j)))){
        loglik_tj(t,j) = dpois( y_tj(t,j), mu_tj(t,j), true );
      }
      SIMULATE{
        y_tj(t,j) = rpois( mu_tj(t,j) );
      }
      devresid_tj(t,j) = sign(y_tj(t,j) - mu_tj(t,j)) * pow(2*(y_tj(t,j)*log((Type(1e-10) + y_tj(t,j))/mu_tj(t,j)) - (y_tj(t,j)-mu_tj(t,j))), 0.5);
    }
    if( familycode_j(j)==4 ){
      // familycode = 4 :  Gamma:   shape = 1/CV^2; scale = mean*CV^2
      if(R_FINITE(asDouble(y_tj(t,j)))){
        loglik_tj(t,j) = dgamma( y_tj(t,j), pow(sigma_z(sigmastart_j(j)),-2), mu_tj(t,j)*pow(sigma_z(sigmastart_j(j)),2), true );
      }
      SIMULATE{
        y_tj(t,j) = rgamma( pow(sigma_z(sigmastart_j(j)),-2), mu_tj(t,j)*pow(sigma_z(sigmastart_j(j)),2) );
      }
      devresid_tj(t,j) = sign(y_tj(t,j) - mu_tj(t,j)) * pow(2 * ( (y_tj(t,j)-mu_tj(t,j))/mu_tj(t,j) - log(y_tj(t,j)/mu_tj(t,j)) ), 0.5);
    }
    if( familycode_j(j)==5 ){
      // familycode = 5 :  normal with known standard deviation
      if(R_FINITE(asDouble(y_tj(t,j)))){
        loglik_tj(t,j) = dnorm( y_tj(t,j), mu_tj(t,j), eps_tj(t,j), true );
      }
      SIMULATE{
        y_tj(t,j) = rnorm( mu_tj(t,j), eps_tj(t,j) );
      }
      devresid_tj(t,j) = NAN;
    }
    if( familycode_j(j)==6 ){
      // familycode = 6 :  lognormal
      if(R_FINITE(asDouble(y_tj(t,j)))){
        loglik_tj(t,j) = dlnorm( y_tj(t,j), log(mu_tj(t,j)), sigma_z(sigmastart_j(j)), true );
      }
      SIMULATE{
        y_tj(t,j) = exp(rnorm( log(mu_tj(t,j)), sigma_z(sigmastart_j(j)) ));
      }
      devresid_tj(t,j) = log(y_tj(t,j)) - log(mu_tj(t,j));
    }
    if( familycode_j(j)==7 ){
      // familycode = 7 :  tweedie
      if(R_FINITE(asDouble(y_tj(t,j)))){
        loglik_tj(t,j) = dtweedie( y_tj(t,j), mu_tj(t,j), exp(sigma_z(sigmastart_j(j))), 1.0 + invlogit(sigma_z(sigmastart_j(j)+1)), true );
      }
      SIMULATE{
        y_tj(t,j) = rtweedie( mu_tj(t,j), exp(sigma_z(sigmastart_j(j))), 1.0 + invlogit(sigma_z(sigmastart_j(j)+1)) );
      }
      devresid_tj(t,j) = devresid_tweedie( y_tj(t,j), mu_tj(t,j), 1.0 + invlogit(sigma_z(sigmastart_j(j)+1)) );
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
