#include <RcppArmadillo.h>
#include "cpp_auxiliary.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Orthogonal Procrustes Analysis : return argmin |X-YQ|
arma::mat OrthogonalProcrustes(arma::mat X, arma::mat Y){
  arma::mat XtY = arma::trans(X)*Y;
  
  arma::mat U;
  arma::vec s;
  arma::mat V;
  
  arma::svd(U,s,V,XtY);
  return(U*V.t());
}

// Pseudo-Inverse A^{-1/2}
arma::mat PseudoSqrtInv(arma::mat X){
  int p = X.n_rows;
  
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, X);
  
  double thrval = 100.0*arma::datum::eps;
  arma::vec newval(p,fill::zeros);
  for (int i=0; i<p; i++){
    if (eigval(i) > thrval){
      newval(i) = 1.0/std::sqrt(eigval(i));
    } else {
      newval(i) = 0.0;
    }
  }
  
  return(eigvec*arma::diagmat(newval)*arma::trans(eigvec));
}



// Centering the matrix & the list of matrices
arma::mat Centering(arma::mat X){
  arma::rowvec y = arma::mean(X, 0);
  arma::mat Y = X.each_row() - y;
  return(Y);
}
arma::field<arma::mat> CenteringField(arma::field<arma::mat> config_list){
  int N = config_list.n_elem;
  arma::field<arma::mat> output(N);
  arma::mat tmp_mat;
  for (int n=0; n<N; n++){
    tmp_mat.reset();
    tmp_mat = Centering(config_list(n));
    output(n) = tmp_mat;
  }
  return(output);
}

// Kernel Computation
// [[Rcpp::export]]
arma::mat KernelLinear(arma::mat X){
  int N = X.n_rows;
  
  arma::mat output(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      output(i,j) = arma::dot(X.row(i), X.row(j));
      output(j,i) = output(i,j);
    }
  }
  for (int n=0; n<N; n++){
    output(n,n) = arma::dot(X.row(n), X.row(n));
  }
  return(output);
}
// [[Rcpp::export]]
arma::mat KernelRBF(arma::mat X, bool use_auto, double bandwidth){
  // parameters
  int n = X.n_rows;

  // compute pairwise distances
  arma::vec upper_dist(n*(n-1)/2, fill::zeros);
  arma::mat mat_dist(n,n,fill::zeros);
  
  int counter = 0;
  for (int i=0; i<(n-1); i++){
    for (int j=(i+1); j<n; j++){
      mat_dist(i,j) = arma::norm(X.row(i)-X.row(j), 2);
      mat_dist(j,i) = mat_dist(i,j);
      
      upper_dist(counter) = mat_dist(i,j);
      counter += 1;
    }
  }
  
  // bandwidth for exp(-|x-y|^2/(2*(sig^2)))
  double par_band = 0.0;
  if (use_auto==false){
    par_band = bandwidth;
  } else {
    par_band = arma::median(upper_dist);
  }
  
  // compute
  return(arma::exp(-(arma::pow(mat_dist, 2)/(2.0*par_band*par_band))));
}
arma::cube Kernel3d(arma::field<arma::mat> config_list, std::string kern_type, bool use_auto, double sigma){
  // params
  int N = config_list.n_elem;
  int m = config_list(0).n_rows;
  
  // prepare
  arma::cube Kcube(m,m,N,fill::zeros);
  arma::field<arma::mat> ctd_list(N);
  for (int n=0; n<N; n++){
    ctd_list(n) = Centering(config_list(n));
  }
  
  // case branching
  if (kern_type=="linear"){
    for (int n=0; n<N; n++){
      Kcube.slice(n) = KernelLinear(ctd_list(n));
    }
  } else {
    for (int n=0; n<N; n++){
      Kcube.slice(n) = KernelRBF(ctd_list(n), use_auto, sigma);
    }
  }
  return(Kcube);
}



double HSIC_pairwise(arma::mat K, arma::mat L){
  int n = K.n_rows;
  double nn = static_cast<double>(n);
  arma::mat H = arma::eye(n,n) - (1.0/nn)*arma::ones(n,n);
  
  double term1 = arma::trace(K*H*L*H);
  double term2 = (nn-1.0)*(nn-1.0);
  return(term1/term2);
}