#' Data : Gorilla Skull
#' 
#' This is 29 male and 30 female gorillas' skull landmark data where each 
#' individual is represented as 8-ad/landmarks in 2 dimensions. This is a 
#' re-arranged version of the data from \pkg{shapes} package where 
#' each representation is mean-centered and normalized by its Frobenius norm.
#' 
#' @usage data(gorilla)
#' 
#' @format a named list containing\describe{
#' \item{data}{a list of \eqn{59} landmark representation matrices in \eqn{\mathbf{R}^{8\times 2}}.}
#' \item{label}{a length-\eqn{59} vector of class labels; either \code{"male"} or \code{"female"}.}
#' }
#' 
#' @examples 
#' # load the data
#' data("gorilla", package="repsim")
#' 
#' # show 2-dimensional landmarks
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3), pty="s")
#' for (i in 1:3){
#'   plot(gorilla$data[[i]], cex=0.5, main=paste0("male #",i), 
#'        xlab="DIM 1", ylab="DIM 2", xaxt="n", yaxt="n")
#'   lines(gorilla$data[[i]])
#' }
#' for (j in 30:32){
#'   plot(gorilla$data[[j]], cex=0.5, main=paste0("female #",j-29), 
#'        xlab="DIM 1", ylab="DIM 2", xaxt="n", yaxt="n")
#'   lines(gorilla$data[[j]])
#' }
#' par(opar)
#' 
#' 
#' @references 
#' Reno PL, Meindl RS, McCollum MA, Lovejoy CO (2003). "Sexual dimorphism in Australopithecus afarensis was similar to that of modern humans." \emph{Proceedings of the National Academy of Sciences}, 100(16):9404–9409.
#' 
#' @concept data
"gorilla"