% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dis_GSP.R
\name{GSP}
\alias{GSP}
\title{Generalized Shape Metric invariant to Permutation}
\usage{
GSP(dlist, ...)
}
\arguments{
\item{dlist}{a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.}

\item{...}{optional parameters including\describe{
\item{centering}{a logical to apply centering (default:\code{TRUE}).}
\item{ndim}{an integer-valued dimension for the common feature space mapping (default: 2).}
}}
}
\value{
an \eqn{(N\times N)} matrix of pairwise dissimilarity measures.
}
\description{
\insertCite{williams_2021_GeneralizedShapeMetrics;textual}{repsim} proposed 
a series of algorithms named as generalized shape metric. This function implements 
the \bold{permutation-invariant} metric.
}
\examples{
\donttest{
## load the 'twoclass' data
data("twoclass", package="repsim")

## seperate data
data1 = twoclass$data

## copy the data by shuffling columns
data2 = vector("list", length=length(data1))
for (i in 1:length(data1)){
  i_data = data1[[i]]
  i_idx  = sample(1:base::ncol(i_data), base::ncol(i_data), replace=FALSE)
  data2[[i]] = i_data[,i_idx]
}

## compute GSP for the two data with default options
gsp1 <- GSP(data1)
gsp2 <- GSP(data2)

## visualize
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2), pty="s")
fields::imagePlot(gsp1, xaxt="n", yaxt="n", horizontal=TRUE, main="GSP-data")
fields::imagePlot(gsp2, xaxt="n", yaxt="n", horizontal=TRUE, main="GSP-shuffled")
par(opar)
}

}
\references{
\insertAllCited{}
}
\concept{dis}
