#' Data : Two Class
#' 
#' This is a generated data that are perturbed versions of two model shapes, 
#' \emph{spiral} and \emph{cassini}. The data contains 25 perturbed versions 
#' of each shape embedded in 5-dimensional space.
#' 
#' @usage data(twoclass)
#' 
#' @format a named list containing\describe{
#' \item{data}{a list of \eqn{50} landmark representation matrices in \eqn{\mathbf{R}^{100\times 5}}.}
#' \item{label}{a length-\eqn{50} vector of class labels; either \code{"spiral"} or \code{"cassini"}.}
#' }
#' 
#' @examples 
#' # load the data
#' data("twoclass", package="repsim")
#' 
#' # take exemplary images from each class
#' im1 <- twoclass$data[[1]]
#' im2 <- twoclass$data[[50]]
#' 
#' # apply PCA for low-dimensional embedding
#' embed1 <- im1%*%eigen(stats::cov(im1))$vectors[,1:2]
#' embed2 <- im2%*%eigen(stats::cov(im2))$vectors[,1:2]
#' 
#' # visualize 2-dimensional embeddings
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(embed1, xaxt='n', yaxt='n', xlab='', ylab='', pch=19, cex=0.5, main="spiral")
#' plot(embed2, xaxt='n', yaxt='n', xlab='', ylab='', pch=19, cex=0.5, main="cassini")
#' par(opar)
#' 
#' @concept data
"twoclass"

# library(mlbench)
# shape1 = mlbench.1spiral(100)
# shape2 = mlbench.cassini(100)$x
# 
# data  = vector("list", length=50L)
# label = rep(c("spiral","cassini"), each=25)
# 
# for (i in 1:25){
#   tmp1 = shape1 + matrix(rnorm(100*2, sd=0.1), ncol=2)
#   tmp1 = cbind(tmp1, array(0,c(100,3)))
#   data[[i]] = runif(1, min=0.5, max=5)*(tmp1%*%qr.Q(qr(matrix(rnorm(25),ncol=5))))
# }
# for (i in 26:50){
#   tmp2 = shape2 + matrix(rnorm(100*2, sd=0.1), ncol=2)
#   tmp2 = cbind(tmp2, array(0,c(100,3)))
#   data[[i]] = runif(1, min=0.5, max=5)*(tmp2%*%qr.Q(qr(matrix(rnorm(25),ncol=5))))
# }
# 
# twoclass = list(data=data, label=label)
