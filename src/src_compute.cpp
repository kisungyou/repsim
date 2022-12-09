/* GENERIC COMPUTATION ROUTINES EXPOSED TO R
 * 
 * src_orthobase : get an orthonormal base
 */

#include <RcppArmadillo.h>
#include "cpp_auxiliary.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat src_orthobase(arma::mat &X, bool centering){
  // case branching by centering
  if (centering){
    arma::mat Y = Centering(X);  
    arma::mat Z = arma::orth(Y);
    return(Z);
  } else {
    arma::mat Z = arma::orth(X);
    return(Z);
  }
}