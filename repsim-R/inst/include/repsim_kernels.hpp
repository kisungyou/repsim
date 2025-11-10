#pragma once
#include "repsim_util.hpp"
#include <string>

//  Kernel
//    kernel_linear ------ linear kernel
//    kernel_rbf --------- median heuristic for RBF kernel
//    kernel_rbf_mean ---- mean heuristic for RBF kernel
//    kernel_rbf_dual ---- dual median form RBF kernel
//    kernel_select ------ select kernel by string ========== update iteratively

namespace repsim
{

    // --------------------------------
    // [3] Kernel
    // --------------------------------
    inline Mat kernel_linear(const Mat &X)
    {
        Mat X_centered = util_centering(X);
        return (X_centered * X_centered.transpose());
    }

    inline Mat kernel_rbf(const Mat &X)
    {
        const Mat D = util_distmat(X);
        double sigma = util_upper_median(D);
        double sig2 = std::max(sigma * sigma, 1e-12);

        const Mat D2 = D.array().square().matrix();
        Mat K = (-D2.array() / (2.0 * sig2)).exp().matrix();
        K.diagonal().setOnes();
        return K;
    }

    inline Mat kernel_rbf_mean(const Mat &X)
    {
        const Mat D2 = util_distmat(X).array().square().matrix();
        const double mean_sq = util_upper_mean(D2);
        const double denom = std::max(mean_sq, 1e-12);

        Mat K = (-D2.array() / denom).exp().matrix();
        K.diagonal().setOnes();
        return K;
    }

    inline Mat kernel_rbf_dualmed(const Mat &X)
    {
        const Index N = X.rows();
        const Mat D = util_distmat(X);

        std::vector<double> row_medians;
        row_medians.reserve(N);
        for (Index i = 0; i < N; i++)
        {
            std::vector<double> row_vals;
            row_vals.reserve(N - 1);
            for (Index j = 0; j < N; j++)
            {
                if (i != j)
                    row_vals.emplace_back(D(i, j));
            }
            row_medians.emplace_back(util_median(row_vals));
        }

        Mat K = Mat::Zero(N, N);
        for (Index i = 0; i < (N - 1); i++)
        {
            for (Index j = (i + 1); j < N; j++)
            {
                double denom = 2.0 * row_medians[i] * row_medians[j];
                denom = std::max(denom, 1e-12);
                const double Dij = D(i, j);
                const double val = std::exp(-(Dij * Dij) / denom);
                K(i, j) = val;
                K(j, i) = val;
            }
        }
        K.diagonal().setOnes();
        return K;
    }

    // kernel_select
    inline Mat kernel_select(const Mat &X, const std::string &type)
    {
        if (type == "linear")
            return kernel_linear(X);
        else if (type == "rbf")
            return kernel_rbf(X);
        else if (type == "rbf_mean")
            return kernel_rbf_mean(X);
        else if (type == "rbf_dualmed")
            return kernel_rbf_dualmed(X);
        throw std::invalid_argument("kernel_select: unknown kernel type '" + type + "'");
    }

} // namespace repsim
