/* DISSIMILARITY MEASURES
 * 
 * (1) cpp_GSR_procrustes : Williams et al. (2021) : rotation-invariant - procrustes
 * (2) cpp_GSR_kendall    : Williams et al. (2021) : rotation-invariant - kendall
 * (3) cpp_GSL            : Williams et al. (2021) : linear-invariant
 */

#include <RcppArmadillo.h>
#include "cpp_auxiliary.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// (1) cpp_GSR_procrustes ------------------------------------------------------
// [[Rcpp::export]]
arma::mat cpp_GSR_procrustes(arma::cube &data3d){
  // params
  int N = data3d.n_slices;
  int m = data3d.n_rows;
  int p = data3d.n_cols;
  
  // centering
  arma::mat Cm = arma::eye(m,m) - (1.0/static_cast<double>(m))*arma::ones(m,m);
  arma::cube array3d(m,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    array3d.slice(n) = Cm*data3d.slice(n);
  }
  
  // temporary variables
  arma::mat X(m,p,fill::zeros);
  arma::mat Y(m,p,fill::zeros);
  arma::mat Q(p,p,fill::zeros);
  
  // prepare an output
  arma::mat output(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    X = array3d.slice(i);
    for (int j=(i+1); j<N; j++){
      Y = array3d.slice(j);
      Q = OrthogonalProcrustes(X,Y);
      
      output(i,j) = arma::norm(X-(Y*Q),"fro");
      output(j,i) = output(i,j);
    }
  }
  return(output);
}

// (2) cpp_GSR_kendall ---------------------------------------------------------
// [[Rcpp::export]]
arma::mat cpp_GSR_kendall(arma::cube &data3d){
  // params
  int N = data3d.n_slices;
  int m = data3d.n_rows;
  int p = data3d.n_cols;
  
  // centering & compute Frobenius norm
  arma::mat Cm = arma::eye(m,m) - (1.0/static_cast<double>(m))*arma::ones(m,m);
  arma::vec Fnorm(N,fill::zeros);
  arma::cube array3d(m,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    array3d.slice(n) = Cm*data3d.slice(n);
    Fnorm(n) = arma::norm(array3d.slice(n), "fro");
  }
  
  // temporary variables
  arma::mat X(m,p,fill::zeros);
  arma::mat Y(m,p,fill::zeros);
  arma::mat Q(p,p,fill::zeros);
  
  // prepare an output
  double dval = 0.0;
  arma::mat output(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    X = array3d.slice(i);
    for (int j=(i+1); j<N; j++){
      Y = array3d.slice(j);
      Q = OrthogonalProcrustes(X,Y);
      
      dval = arma::dot(X, Y*Q)/(Fnorm(i)*Fnorm(j));
      
      output(i,j) = std::acos(dval);
      output(j,i) = output(i,j);
    }
  }
  return(output);
}

// (3) cpp_GSL -----------------------------------------------------------------
// [[Rcpp::export]]
arma::mat cpp_GSL(arma::cube &data3d, double alpha){
  // params
  int N = data3d.n_slices;
  int m = data3d.n_rows;
  int p = data3d.n_cols;
  
  // feature mapping & compute Frobenius norms
  arma::mat Ip(p,p,fill::eye);
  arma::mat ciq(p,p,fill::zeros);
  
  arma::vec Fnorm(N,fill::zeros);
  arma::mat Cm = arma::eye(m,m) - (1.0/static_cast<double>(m))*arma::ones(m,m);
  arma::cube array3d(m,p,N,fill::zeros);
  if (alpha >= (1.0 - arma::datum::eps)){
    for (int n=0; n<N; n++){
      array3d.slice(n) = Cm*data3d.slice(n);
      Fnorm(n) = arma::norm(array3d.slice(n), "fro");
    }
  } else {
    for (int n=0; n<N; n++){
      ciq = PseudoSqrtInv(arma::trans(data3d.slice(n))*Cm*data3d.slice(n));
      array3d.slice(n) = Cm*data3d.slice(n)*(alpha*Ip + (1.0-alpha)*ciq);
      Fnorm(n) = arma::norm(array3d.slice(n), "fro");
    }
  }
  
  // temporary variables
  arma::mat X(m,p,fill::zeros);
  arma::mat Y(m,p,fill::zeros);
  arma::mat Q(p,p,fill::zeros);
  
  // prepare an output
  double dval = 0.0;
  arma::mat output(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    X = array3d.slice(i);
    for (int j=(i+1); j<N; j++){
      Y = array3d.slice(j);
      Q = OrthogonalProcrustes(X,Y);
      
      dval = arma::dot(X, Y*Q)/(Fnorm(i)*Fnorm(j));
      
      output(i,j) = std::acos(dval);
      output(j,i) = output(i,j);
    }
  }
  return(output);
}