% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_HSIC.R
\name{HSIC}
\alias{HSIC}
\title{Hilbert-Schmidt Independence Criterion}
\usage{
HSIC(dlist, ...)
}
\arguments{
\item{dlist}{a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.}

\item{...}{optional parameters including\describe{
\item{kernel}{(case-insensitive) name of the kernel to be used; either \code{"linear"} or \code{"rbf"} (default: \code{"rbf"}).}
\item{sigma}{bandwidth parameter for the RBF kernel when \code{kernel="rbf"}. If \code{NULL}, it automatically selects 
a bandwidth as the median of pairwise distances. Otherwise, a number should be provided (default: \code{NULL}).}
}}
}
\value{
an \eqn{(N\times N)} matrix of pairwise similarity measures.
}
\description{
Given two centered representation matrices \eqn{X \in \mathbf{R}^{m\times p}} and \eqn{Y \in \mathbf{R}^{m\times q}}, 
the Hilbert-Schmidt Independence Criterion (HSIC) is defined as 
\deqn{\textrm{HSIC}(X,Y) = \frac{1}{(m-1)^2} \textrm{tr}(KHLH),}
where \eqn{K} and \eqn{L} are \eqn{(m\times m)} kernel matrices of \eqn{X} and \eqn{Y}, respectively and \eqn{H} is a centering matrix. 
Currently, our implementation supports linear and RBF kernels.
}
\examples{
\donttest{
## load the 'twoclass' data
data("twoclass", package="repsim")

## compare linear & rbf kernels
HSIClin <- HSIC(twoclass$data, kernel="linear")
HSICrbf <- HSIC(twoclass$data, kernel="rbf")

## visualize
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2), pty="s")
fields::imagePlot(HSIClin, main="HSIC-linear", xaxt="n", yaxt="n", horizontal=TRUE)
fields::imagePlot(HSICrbf, main="HSIC-rbf",    xaxt="n", yaxt="n", horizontal=TRUE)
par(opar)
}

}
\references{
\insertRef{gretton_2005_MeasuringStatisticalDependence}{repsim}
}
\concept{sim}
