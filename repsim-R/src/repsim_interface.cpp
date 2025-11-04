// Functions to be exposed
//  cpp_DotProduct
//  cpp_LinReg

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include "repsim_core.hpp"

using Eigen::MatrixXd;

// Convert an R list of numeric matrices to std::vector<repsim::Mat> (row-major)
static inline std::vector<repsim::Mat> list_to_eigens(Rcpp::List mats) {
  const int m = mats.size();
  std::vector<repsim::Mat> Xs;
  Xs.reserve(m);
  for (int i = 0; i < m; ++i) {
    // Simple and robust: convert SEXP -> Eigen::MatrixXd, then copy to row-major
    MatrixXd Xi = Rcpp::as<MatrixXd>(mats[i]);
    Xs.emplace_back(Xi);  // copy into repsim::Mat (row-major)
  }
  return Xs;
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_DotProduct(Rcpp::List mats){
  std::vector<repsim::Mat> Xs = list_to_eigens(mats);
  return(repsim::core_DotProduct(Xs));
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_LinReg(Rcpp::List mats){
  std::vector<repsim::Mat> Xs = list_to_eigens(mats);
  return(repsim::core_LinReg(Xs));
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_HSIC(Rcpp::List mats, std::string type, std::string estimator = "gretton"){
  std::vector<repsim::Mat> Xs = list_to_eigens(mats);
  return(repsim::core_HSIC(Xs, type, estimator));
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_CKA(Rcpp::List mats, std::string type, std::string estimator = "gretton"){
  std::vector<repsim::Mat> Xs = list_to_eigens(mats);
  return(repsim::core_CKA(Xs, type, estimator));
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_CCA(Rcpp::List mats, std::string type){
  std::vector<repsim::Mat> Xs = list_to_eigens(mats);
  return(repsim::core_CCA(Xs, type));
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_SVCCA(Rcpp::List mats, std::string type){
  std::vector<repsim::Mat> Xs = list_to_eigens(mats);
  return(repsim::core_SVCCA(Xs, type));
}

// [[Rcpp::export]]
Eigen::MatrixXd cpp_PWCCA(Rcpp::List mats){
  std::vector<repsim::Mat> Xs = list_to_eigens(mats);
  return(repsim::core_PWCCA(Xs));
}