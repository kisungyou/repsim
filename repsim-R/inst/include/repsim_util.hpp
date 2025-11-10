#pragma once

#include <Eigen/Core>
#include <Eigen/SVD>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <string>
#include <cmath>
#include <vector>

//  Utilities
//    util_median --------- compute the median of a vector
//    util_ONB ------------ polar decomposition to get orthonormal basis
//    util_distmat -------- compute pairwise distance matrix
//    util_upper_median --- given D, compute median of all entries in the upper triangle
//    util_upper_mean ----- given D, compute mean of all entries in the upper triangle
//    util_nuclear_norm --- compute nuclear norm of a matrix
//    util_centering ------ return a mean centered version of input matrix
//    util_SVD_denoise ---- denoise by explained variance (centered input case)
//    util_HSIC_estimator - different types of estimators for HSIC

namespace repsim
{

    // --------------------------------
    // [1] Basic Typedefs
    // --------------------------------
    using Index = Eigen::Index;
    using Mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using Vec = Eigen::VectorXd;

    // --------------------------------
    // [2] Utilities
    // --------------------------------
    inline double util_median(std::vector<double> vals)
    {
        const size_t n = vals.size();
        if (n == 0)
            return std::numeric_limits<double>::quiet_NaN();

        auto mid = vals.begin() + n / 2;
        std::nth_element(vals.begin(), mid, vals.end());
        const double upper_mid = *mid;

        if (n % 2 == 1)
            return upper_mid;

        const auto lower_mid_it = std::max_element(vals.begin(), mid);
        const double lower_mid = *lower_mid_it;
        return 0.5 * (lower_mid + upper_mid);
    }

    inline Mat util_ONB(const Mat &X)
    {
        double rel_tol = 1e-12;

        Eigen::JacobiSVD<Mat> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const Eigen::VectorXd s = svd.singularValues();
        const double smax = (s.size() > 0) ? s.maxCoeff() : 0.0;

        const double thresh = rel_tol * smax;
        Eigen::Index r = 0;
        for (Eigen::Index k = 0; k < s.size(); ++k)
            if (s[k] > thresh)
                ++r;

        if (r == 0)
            return Mat::Zero(X.rows(), 0);
        return svd.matrixU().leftCols(r);
    }

    inline Mat util_distmat(const Mat &X)
    {
        Mat G = X * X.transpose();
        Vec d = G.diagonal();
        Mat D = (-2.0 * G);
        D.colwise() += d;
        D.rowwise() += d.transpose();
        D = D.cwiseMax(0.0);
        D = D.array().sqrt().matrix();
        return D;
    }

    inline double util_upper_median(const Mat &D)
    {
        const Index n = D.rows();
        std::vector<double> vals;
        vals.reserve(static_cast<size_t>(n) * (n - 1) / 2);
        for (Index i = 0; i < (n - 1); i++)
        {
            for (Index j = (i + 1); j < n; j++)
            {
                vals.emplace_back(D(i, j));
            }
        }
        if (vals.empty())
            return 0.0;
        return util_median(vals);
    }

    inline double util_upper_mean(const Mat &D)
    {
        const Index n = D.rows();
        if (n < 2)
            return 0.0;

        double sum = 0.0;
        size_t count = 0;
        for (Index i = 0; i < n - 1; ++i)
        {
            for (Index j = i + 1; j < n; ++j)
            {
                sum += D(i, j);
                ++count;
            }
        }
        return (count > 0) ? (sum / static_cast<double>(count)) : 0.0;
    }

    inline double util_nuclear_norm(const Mat &X)
    {
        Eigen::JacobiSVD<Mat> svd(X);
        return svd.singularValues().sum();
    }

    inline Mat util_centering(const Mat &X)
    {
        return X.rowwise() - X.colwise().mean();
    }

    inline Mat util_SVD_denoise(const Mat &X, double var_thresh = 0.99)
    {
        if (var_thresh < 0.0)
            var_thresh = 1e-9;
        if (var_thresh > 1.0)
            var_thresh = 1.0 - 1e-9;

        const Index n = X.rows(), p = X.cols();
        if (n == 0 || p == 0)
            return Mat::Zero(n, p);

        Eigen::JacobiSVD<Mat> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const Vec s = svd.singularValues();
        if (s.size() == 0)
            return Mat::Zero(n, p);

        const Vec s2 = s.array().square().matrix();
        const double total = s2.sum();
        if (total <= 0.0)
            return Mat::Zero(n, p);

        Index r = 1;
        double cum = s2[0];
        while (r < s2.size() && (cum / total) < var_thresh)
        {
            cum += s2[r];
            ++r;
        }

        const Mat U_r = svd.matrixU().leftCols(r);
        const Mat V_r = svd.matrixV().leftCols(r);
        const Mat Sr = s.head(r).asDiagonal();
        return U_r * Sr * V_r.transpose();
    }

    // --- HSIC helpers (kept in util so kernels/core can share) -------------------

    inline double util_HSIC_gretton(const Mat &Kx, const Mat &Ky)
    {
        const Index N = Kx.rows();
        const Mat H = Mat::Identity(N, N) - (Mat::Ones(N, N) / static_cast<double>(N));
        const double NN = static_cast<double>(N);
        const double denom = (NN - 1.0) * (NN - 1.0);
        return (Kx * H * Ky * H).trace() / denom;
    }

    inline double util_HSIC_song(const Mat &Kx_, const Mat &Ky_)
    {
        const Index n = Kx_.rows();
        if (n < 4)
            return util_HSIC_gretton(Kx_, Ky_);

        Mat Kx = Kx_;
        Mat Ky = Ky_;
        Kx.diagonal().setZero();
        Ky.diagonal().setZero();

        const double nd = static_cast<double>(n);
        const double denom_main = nd * (nd - 3.0);
        const double inv_nm1_nm2 = 1.0 / ((nd - 1.0) * (nd - 2.0));
        const double coef_cross = 2.0 / (nd - 2.0);

        const double term1 = (Kx.cwiseProduct(Ky)).sum();
        Vec ones = Vec::Ones(n);
        const Vec sK = Kx * ones;
        const Vec sL = Ky * ones;

        const double sumK = sK.sum();
        const double sumL = sL.sum();
        const double term2 = (sumK * sumL) * inv_nm1_nm2;
        const double term3 = coef_cross * sK.dot(sL);

        return (term1 + term2 - term3) / denom_main;
    }

    inline double util_HSIC_lange(const Mat &Kx, const Mat &Ky)
    {
        const Index n = Kx.rows();
        if (n < 4)
            return util_HSIC_gretton(Kx, Ky);

        const double nd = static_cast<double>(n);
        Mat H = Mat::Identity(n, n) - Mat::Ones(n, n) / nd;
        Mat Kc = H * Kx * H;
        Mat Lc = H * Ky * H;

        double s = 0.0;
        for (Index i = 1; i < n; ++i)
        {
            s += Kc.row(i).head(i).dot(Lc.row(i).head(i));
        }
        return (2.0 / (nd * (nd - 3.0))) * s;
    }

    inline double util_HSIC_estimator(const Mat &Kx, const Mat &Ky,
                                      const std::string &type = "gretton")
    {
        if (type == "gretton")
            return util_HSIC_gretton(Kx, Ky);
        if (type == "song")
            return util_HSIC_song(Kx, Ky);
        if (type == "lange")
            return util_HSIC_lange(Kx, Ky);
        throw std::invalid_argument("util_HSIC_estimator: unknown type '" + type + "'");
    }

} // namespace repsim
