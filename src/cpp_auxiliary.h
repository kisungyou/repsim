#ifndef _genshapes_CPP_AUXILIARY_H
#define _genshapes_CPP_AUXILIARY_H

#define ARMA_NO_DEBUG

// Orthogonal Procrustes Analysis : return argmin |X-YQ|
arma::mat OrthogonalProcrustes(arma::mat X, arma::mat Y);

// Pseudo-Inverse A^{-1/2}
arma::mat PseudoSqrtInv(arma::mat X);

// Centering the matrix
arma::mat Centering(arma::mat X);
arma::field<arma::mat> CenteringField(arma::field<arma::mat> config_list);

// Kernel Computations : Single
arma::mat  KernelLinear(arma::mat X);
arma::mat  KernelRBF(arma::mat X, bool use_auto, double bandwidth);

// Kernel Computation : all using {linear} or {rbf:auto/fixed bandwidth}
arma::cube Kernel3d(arma::field<arma::mat> config_list, 
                    std::string kern_type,
                    bool use_auto,
                    double sigma);

// HSIC given two kernel matrices
double HSIC_pairwise(arma::mat K, arma::mat L);
                    


#endif