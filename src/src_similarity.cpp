/* SIMILARITY MEASURES
 * 
 * (1) cpp_LR   : Linear Regression-Based Similarity
 * (2) cpp_DP   : Dot Product-Based Similarity
 * (3) cpp_HSIC : Hilbert-Schmidt Independence Criterion
 * (4) cpp_CKA  : Centered Kernel Alignment
 */

#include <RcppArmadillo.h>
#include "cpp_auxiliary.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat cpp_LR(arma::field<arma::mat> &config_list, bool centering){
  // STEP 0 : get parameters
  int N = config_list.n_elem;
  
  // STEP 1 : mean centering, extract Q, compute Frobenius norm
  arma::field<arma::mat> list_centered(N);
  arma::field<arma::mat> list_Q(N);
  arma::vec list_F(N,fill::zeros);
  arma::mat tmp_mat; 
  arma::mat tmp_Q;
  
  for (int n=0; n<N; n++){
    // centering each matrix (if necessary)
    if (centering){
      tmp_mat = Centering(config_list(n));  
    } else {
      tmp_mat = config_list(n);
    }

    // QR Decomposition
    tmp_Q = arma::orth(tmp_mat);
    
    // assign
    list_centered(n) = tmp_mat;
    list_Q(n) = tmp_Q;
    
    // compute Frobenius-norm
    list_F(n) = arma::norm(tmp_mat, "fro");
    
    // reset
    tmp_mat.reset();
    tmp_Q.reset();
  }
  
  // STEP 2 : iterate
  double dval = 0.0;
  arma::mat output(N,N,fill::ones);
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      if (i!=j){
        dval = arma::norm(arma::trans(list_Q(j))*list_centered(i), "fro")/list_F(i);
        output(i,j) = dval*dval;
      }
    }
  }
  
  // STEP 3 : return
  return(output);
}


// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat cpp_DP(arma::field<arma::mat> &config_list, bool centering){
  // STEP 0 : get parameters
  int N = config_list.n_elem;
  
  // STEP 1 : mean centering
  arma::field<arma::mat> list_centered(N);
  arma::mat tmp_mat;

  // centering or not
  if (centering){
    for (int n=0; n<N; n++){
      tmp_mat = Centering(config_list(n));
      list_centered(n) = tmp_mat;
      tmp_mat.reset();
    }
  } 
  
  // STEP 2 : iterate
  double dval = 0.0;
  arma::mat output(N,N,fill::zeros);
  for (int n=0; n<N; n++){
    if (centering){
      dval = arma::norm(arma::trans(list_centered(n))*list_centered(n), "fro");  
    } else {
      dval = arma::norm(arma::trans(config_list(n))*config_list(n), "fro");  
    }
    output(n,n) = dval*dval;
  }
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      if (centering){
        dval = arma::norm(arma::trans(list_centered(i))*list_centered(j), "fro");  
      } else {
        dval = arma::norm(arma::trans(config_list(i))*config_list(j), "fro");  
      }
      
      output(i,j) = dval*dval;
      output(j,i) = dval*dval;
    }
  }
  
  // STEP 3 : return
  return(output);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat cpp_HSIC(arma::field<arma::mat> &config_list, std::string par_kernel,
                   double par_sigma, bool use_auto){
  // STEP 0 : get some params
  int N = config_list.n_elem;
  
  // STEP 1 : compute Kernels
  arma::cube Kcube = Kernel3d(config_list, par_kernel, use_auto, par_sigma);
  
  // STEP 2 : compute HSIC
  arma::mat output(N,N,fill::zeros);
  for (int n=0; n<N; n++){
    output(n,n) = HSIC_pairwise(Kcube.slice(n), Kcube.slice(n));
  }
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      output(i,j) = HSIC_pairwise(Kcube.slice(i), Kcube.slice(j));
      output(j,i) = output(i,j);
    }
  }
  return(output);
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat cpp_CKA(arma::field<arma::mat> &config_list, std::string par_kernel,
                  double par_sigma, bool use_auto){
  // STEP 0 : get some params
  int N = config_list.n_elem;
  
  // STEP 1 : compute Kernels
  arma::cube Kcube = Kernel3d(config_list, par_kernel, use_auto, par_sigma);
  
  // STEP 2 : compute HSIC for selfs
  arma::vec vec_HSIC(N,fill::zeros);
  for (int n=0; n<N; n++){
    vec_HSIC(n) = HSIC_pairwise(Kcube.slice(n), Kcube.slice(n));
  }
  
  // STEP 3 : iterate
  arma::mat output(N,N,fill::ones);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      output(i,j) = HSIC_pairwise(Kcube.slice(i), Kcube.slice(j))/std::sqrt(vec_HSIC(i)*vec_HSIC(j));
      output(j,i) = output(i,j);
    }
  }
  return(output);
}