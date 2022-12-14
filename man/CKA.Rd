% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_CKA.R
\name{CKA}
\alias{CKA}
\title{Centered Kernel Alignment}
\usage{
CKA(dlist, ...)
}
\arguments{
\item{dlist}{a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.}

\item{...}{optional parameters including\describe{
\item{kernel}{(case-insensitive) name of the kernel to be used; either \code{"linear"} or \code{"rbf"} (default: \code{"rbf"}).}
\item{sigma}{bandwidth parameter for the RBF kernel when \code{kernel="rbf"}. If \code{NULL}, it automatically selects 
a bandwidth as the median of pairwise distances. Otherwise, a number should be provided (default: \code{NULL}).
}
}}
}
\value{
an \eqn{(N\times N)} matrix of pairwise similarity measures.
}
\description{
Given two centered representation matrices \eqn{X \in \mathbf{R}^{m\times p}} and \eqn{Y \in \mathbf{R}^{m\times q}}, 
the Centered Kernel Alignment (CKA) is defined using the HSIC as 
\deqn{\textrm{CKA}(X,Y) = \frac{\textrm{HSIC}(X,Y)}{\sqrt{
\textrm{HSIC}(X,X)\cdot \textrm{HSIC}(Y,Y)
}}}
where \eqn{K} and \eqn{L} are \eqn{(m\times m)} kernel matrices of \eqn{X} and \eqn{Y}, respectively and \eqn{H} is a centering matrix. 
Its definition is similar to the relationship between covariance and correlation where the latter is normalized version of the former. 
Currently, our implementation supports linear and RBF kernels. For a linear kernel, CKA is reduced to the RV coefficient.
}
\examples{
\donttest{
## load the 'twoclass' data
data("twoclass", package="repsim")

# compare linear & rbf kernels
CKAlin <- CKA(twoclass$data, kernel="linear")
CKArbf <- CKA(twoclass$data, kernel="rbf")

## visualize
graphics.off()
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2), pty="s")
fields::imagePlot(CKAlin, main="CKA-linear", xaxt="n", yaxt="n", horizontal=TRUE)
fields::imagePlot(CKArbf, main="CKA-rbf",    xaxt="n", yaxt="n", horizontal=TRUE)
par(opar)
}

}
\references{
\insertRef{cortes_2012_AlgorithmsLearningKernels}{repsim}

\insertRef{kornblith_2019_SimilarityNeuralNetwork}{repsim}
}
\concept{sim}
