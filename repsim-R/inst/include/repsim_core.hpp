// MAIN HEADER
//  Core Computation with multiple inputs
//    core_DotProduct ---- squared dot product between inputs
//    core_LinReg -------- R^2
//    core_HSIC ---------- HSIC
//    core_CKA ----------- CKA
//    core_CCA ----------- CCA
//    core_SVCCA --------- singular vector CCA
//    core_PWCCA --------- projection weighted CCA

#pragma once

// Re-export types & utils & kernels
#include "repsim_util.hpp"    // defines repsim::{Index, Mat, Vec, ... utils ...}
#include "repsim_kernels.hpp" // defines repsim::{kernel_* , kernel_select}

#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace repsim
{
  // --------------------------------
  inline Mat core_LinReg(const std::vector<Mat> &Xs)
  {
    // prep
    const Index M = static_cast<Index>(Xs.size());
    Mat output = Mat::Ones(M, M);

    // batch compute
    std::vector<Mat> Qs;
    Qs.reserve(M);
    std::vector<double> X_sqnorms;
    X_sqnorms.reserve(M);
    for (Index i = 0; i < M; i++)
    {
      Qs.emplace_back((util_ONB(util_centering(Xs[i])))); // get ONB of centered X
      X_sqnorms.emplace_back((Xs[i]).squaredNorm());
    }

    // compute off-diagonals
    for (Index i = 0; i < M; i++)
    {
      const Mat &now_Q = Qs[i];
      for (Index j = 0; j < M; j++)
      {
        const Mat &now_X = Xs[j];
        const double denom = X_sqnorms[j];
        if (denom > 0.0)
        {
          output(i, j) = (now_Q.transpose() * now_X).squaredNorm() / denom;
        }
        else
        {
          output(i, j) = 0.0;
        }
      }
    }

    // symmetrize
    for (Index i = 0; i < (M - 1); i++)
    {
      for (Index j = (i + 1); j < M; j++)
      {
        double val_ij = output(i, j);
        double val_ji = output(j, i);
        output(i, j) = (val_ij + val_ji) / 2.0;
        output(j, i) = output(i, j);
      }
    }

    // return the output
    return (output);
  }

  inline Mat core_DotProduct(const std::vector<Mat> &Xs)
  {
    // prep
    const Index M = static_cast<Index>(Xs.size());
    Mat output = Mat::Zero(M, M);

    // compute diagonal
    double val = 0.0;
    for (Index i = 0; i < M; i++)
    {
      val = (Xs[i].transpose() * Xs[i]).squaredNorm();
      output(i, i) = val;
    }

    // compute off-diagonal
    for (Index i = 0; i < (M - 1); i++)
    {
      for (Index j = (i + 1); j < M; j++)
      {
        val = (Xs[j].transpose() * Xs[i]).squaredNorm();
        output(i, j) = val;
        output(j, i) = val;
      }
    }

    // return the output
    return (output);
  }

  inline Mat core_HSIC(const std::vector<Mat> &Xs,
                       const std::string &type,
                       const std::string &estimator)
  {
    // prep
    const Index M = static_cast<Index>(Xs.size());
    Mat output = Mat::Zero(M, M);

    // precompute - kernels at batch
    std::vector<Mat> batch_K;
    batch_K.reserve(M);
    for (Index i = 0; i < M; i++)
    {
      batch_K.emplace_back(kernel_select(Xs[i], type));
    }

    // compute - diagonal
    for (Index m = 0; m < M; m++)
    {
      output(m, m) = util_HSIC_estimator(batch_K[m], batch_K[m], estimator);
    }
    for (Index i = 0; i < (M - 1); i++)
    {
      for (Index j = (i + 1); j < M; j++)
      {
        double val = util_HSIC_estimator(batch_K[i], batch_K[j], estimator);
        output(i, j) = val;
        output(j, i) = val;
      }
    }

    // return
    return (output);
  }

  inline Mat core_CKA(const std::vector<Mat> &Xs,
                      const std::string &kernel_type,
                      const std::string &estimator)
  {

    const Index M = static_cast<Index>(Xs.size());

    // Precompute kernels
    std::vector<Mat> K;
    K.reserve(M);
    for (Index i = 0; i < M; ++i)
      K.emplace_back(kernel_select(Xs[i], kernel_type));

    // Compute HSIC matrix with requested estimator
    Mat HSIC_mat = Mat::Zero(M, M);
    for (Index i = 0; i < M; ++i)
    {
      for (Index j = i; j < M; ++j)
      {
        double v = util_HSIC_estimator(K[i], K[j], estimator);
        HSIC_mat(i, j) = v;
        HSIC_mat(j, i) = v;
      }
    }

    // For normalization, use stable diagonals (prefer Gretton)
    Vec diag = Vec::Zero(M);
    for (Index i = 0; i < M; ++i)
    {
      double d = util_HSIC_gretton(K[i], K[i]); // force Gretton here
      // guard against tiny/negative numerical values
      diag[i] = std::max(d, 1e-12);
    }

    // Build CKA
    Mat out = Mat::Zero(M, M);
    out.diagonal().setOnes();
    for (Index i = 0; i < M - 1; ++i)
    {
      for (Index j = i + 1; j < M; ++j)
      {
        double denom = std::sqrt(diag[i] * diag[j]);
        double val = (denom > 0.0) ? (HSIC_mat(i, j) / denom) : 0.0;
        out(i, j) = val;
        out(j, i) = val;
      }
    }
    return out;
  }

  inline Mat core_CCA(const std::vector<Mat> &Xs,
                      const std::string &type)
  {
    const Index M = static_cast<Index>(Xs.size()); // number of reps
    std::vector<Mat> vec_Q;
    vec_Q.reserve(M); // batch compute ONBs
    std::vector<double> vec_dims;
    vec_dims.reserve(M);
    for (Index m = 0; m < M; m++)
    {
      Mat now_X = util_ONB(Xs[m]);                              // get the current matrix
      vec_dims.emplace_back(static_cast<double>(now_X.cols())); // current dim

      Mat X_center = util_centering(now_X);   // mean center
      vec_Q.emplace_back(util_ONB(X_center)); // get the ONB of centered matrix
    }

    Mat output = Mat::Zero(M, M);
    // compute the diagonals
    for (Index m = 0; m < M; m++)
    {
      const Mat Qx = vec_Q[m];
      if (type == "yanai")
      {
        output(m, m) = (Qx.transpose() * Qx).squaredNorm() / vec_dims[m];
      }
      else if (type == "pillai")
      {
        output(m, m) = util_nuclear_norm(Qx.transpose() * Qx) / vec_dims[m];
      }
      else
      {
        throw std::invalid_argument("core_CCA: unknown CCA type '" + type + "'");
      }
    }

    // compute the off-diagonals
    for (Index i = 0; i < (M - 1); i++)
    {
      const Mat Qx = vec_Q[i];
      for (Index j = (i + 1); j < M; j++)
      {
        const Mat Qy = vec_Q[j];
        if (type == "yanai")
        { // yanai
          double denom = std::min(vec_dims[i], vec_dims[j]);
          double val = (Qx.transpose() * Qy).squaredNorm() / denom;
          output(i, j) = val;
          output(j, i) = val;
        }
        else
        { // just pillai
          double denom = std::min(vec_dims[i], vec_dims[j]);
          double val = util_nuclear_norm(Qx.transpose() * Qy) / denom;
          output(i, j) = val;
          output(j, i) = val;
        }
      }
    }
    return (output);
  }

  inline Mat core_SVCCA(const std::vector<Mat> &Xs,
                        const std::string &type)
  {
    const Index M = static_cast<Index>(Xs.size()); // number of reps
    std::vector<Mat> vec_Q;
    vec_Q.reserve(M); // batch compute ONBs
    std::vector<double> vec_dims;
    vec_dims.reserve(M);
    for (Index m = 0; m < M; m++)
    {
      Mat now_X = util_ONB(Xs[m]);                              // get the current matrix
      vec_dims.emplace_back(static_cast<double>(now_X.cols())); // current dim

      Mat X_center = util_centering(now_X);       // mean center
      Mat X_denoise = util_SVD_denoise(X_center); // svd denoise
      vec_Q.emplace_back(util_ONB(X_denoise));    // get the ONB of denoised matrix
    }

    Mat output = Mat::Zero(M, M);
    // compute the diagonals
    for (Index m = 0; m < M; m++)
    {
      const Mat Qx = vec_Q[m];
      if (type == "yanai")
      {
        output(m, m) = (Qx.transpose() * Qx).squaredNorm() / vec_dims[m];
      }
      else if (type == "pillai")
      {
        output(m, m) = util_nuclear_norm(Qx.transpose() * Qx) / vec_dims[m];
      }
      else
      {
        throw std::invalid_argument("core_SVCCA: unknown CCA type '" + type + "'");
      }
    }

    // compute the off-diagonals
    for (Index i = 0; i < (M - 1); i++)
    {
      const Mat Qx = vec_Q[i];
      for (Index j = (i + 1); j < M; j++)
      {
        const Mat Qy = vec_Q[j];
        if (type == "yanai")
        { // yanai
          double denom = std::min(vec_dims[i], vec_dims[j]);
          double val = (Qx.transpose() * Qy).squaredNorm() / denom;
          output(i, j) = val;
          output(j, i) = val;
        }
        else
        { // just pillai
          double denom = std::min(vec_dims[i], vec_dims[j]);
          double val = util_nuclear_norm(Qx.transpose() * Qy) / denom;
          output(i, j) = val;
          output(j, i) = val;
        }
      }
    }
    return (output);
  }

  // --------------------------------
  // [5] PWCCA (Projection-Weighted CCA)
  // --------------------------------
  //
  // References/recipe:
  // - Center columns of X,Y (n x p), (n x q)
  // - Thin SVDs: Xc = Ux Sx Vx^T,  Yc = Uy Sy Vy^T
  // - Compute M = Ux^T Uy; SVD: M = L Σ R^T
  // - Canonical correlations: diag(Σ) (length r)
  // - Canonical variates (sample-space): UxQx = Ux * L,  UyQy = Uy * R
  // - Projection weights for X-side:
  //     C = (UxQx)^T Xc   [shape: r x p]
  //     w_i = sum_j |C_{i,j}|,  then normalize so sum_i w_i = 1
  // - PWCCA(X,Y) = sum_i w_i * Σ_i
  //
  // Notes:
  // - r = min(rank(Xc), rank(Yc)). Rank is estimated with relative tol to largest singular value.
  // - If total weight is zero (degenerate), fall back to uniform weights over r.
  // - Returns similarity in [0,1]. A distance variant is (1 - PWCCA).

  inline double core_PWCCA_pair(const Mat &X, const Mat &Y)
  {
    // 0) parameter
    double rel_tol = 1e-12;

    // 1) Center columns
    const Mat Xc = util_centering(X);
    const Mat Yc = util_centering(Y);

    const Index n = Xc.rows();
    if (n == 0 || Xc.cols() == 0 || Yc.cols() == 0)
    {
      return 0.0;
    }

    // 2) Thin SVDs
    Eigen::JacobiSVD<Mat> svdX(Xc, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::JacobiSVD<Mat> svdY(Yc, Eigen::ComputeThinU | Eigen::ComputeThinV);

    const Vec sX = svdX.singularValues();
    const Vec sY = svdY.singularValues();
    const double sXmax = (sX.size() > 0 ? sX.maxCoeff() : 0.0);
    const double sYmax = (sY.size() > 0 ? sY.maxCoeff() : 0.0);
    if (sXmax <= 0.0 || sYmax <= 0.0)
    {
      return 0.0;
    }

    // Rank thresholds
    const double tX = rel_tol * sXmax;
    const double tY = rel_tol * sYmax;

    Index rX = 0, rY = 0;
    for (Index i = 0; i < sX.size(); ++i)
      if (sX[i] > tX)
        ++rX;
    for (Index i = 0; i < sY.size(); ++i)
      if (sY[i] > tY)
        ++rY;
    const Index r = std::max<Index>(1, std::min(rX, rY)); // ensure >=1

    // Truncate U to ranks
    const Mat Ux = svdX.matrixU().leftCols(rX); // n x rX
    const Mat Uy = svdY.matrixU().leftCols(rY); // n x rY

    // 3) CCA via SVD on Ux^T Uy
    //    M is rX x rY => SVD: L (rX x r), Σ (r x r), R (rY x r)
    const Mat M = Ux.transpose() * Uy;
    Eigen::JacobiSVD<Mat> svdM(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Mat L = svdM.matrixU().leftCols(r);      // rX x r
    const Vec sig = svdM.singularValues().head(r); // canonical correlations (length r)
    const Mat R = svdM.matrixV().leftCols(r);      // rY x r

    // 4) Sample-space canonical variates for X: UxQx = Ux * L  (n x r)
    const Mat UxQx = Ux * L;

    // 5) Projection coefficients of neurons (columns of Xc) onto UxQx
    //    C = (UxQx)^T * Xc  (r x p)
    const Mat C = UxQx.transpose() * Xc;

    // 6) Weights: rowwise sum of absolute coefficients, then normalize
    Vec w = Vec::Zero(r);
    for (Index i = 0; i < r; ++i)
    {
      w[i] = C.row(i).cwiseAbs().sum();
    }
    double wsum = w.sum();
    if (wsum <= 0.0 || !std::isfinite(wsum))
    {
      // fallback: uniform weights
      w.setOnes();
      wsum = static_cast<double>(r);
    }
    w.array() /= wsum;

    // 7) Weighted sum of canonical correlations
    double pwcca = 0.0;
    for (Index i = 0; i < r; ++i)
    {
      // guard for numerical leakage of corr > 1
      const double rho = std::min(1.0, std::max(0.0, sig[i]));
      pwcca += w[i] * rho;
    }
    return pwcca;
  }

  // Batched PWCCA over a vector of matrices (same n rows each)
  inline Mat core_PWCCA(const std::vector<Mat> &Xs)
  {
    const Index M = static_cast<Index>(Xs.size());
    Mat out = Mat::Zero(M, M);
    // Diagonals = 1 by convention (similarity with itself)
    out.diagonal().setOnes();

    for (Index i = 0; i < M; ++i)
    {
      for (Index j = i + 1; j < M; ++j)
      {
        const double sij = core_PWCCA_pair(Xs[i], Xs[j]);
        // PWCCA is symmetric (if you always weight from the same side,
        // you *can* get asymmetry; common convention is X->Y weights.
        // Here we symmetrize by averaging both directions.)
        const double sji = core_PWCCA_pair(Xs[j], Xs[i]);
        const double s = 0.5 * (sij + sji);
        out(i, j) = s;
        out(j, i) = s;
      }
    }
    return out;
  }
} // namespace repsim