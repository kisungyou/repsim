#' Canonical Correlation Analysis 
#' 
#' The method of Canonical Correlation Analysis (CCA) is one of the most fundamental 
#' methods to measure association between two matrices. This function implements 
#' two goodness of fit measures. One is the mean squared CCA correlation, which is 
#' also known as \emph{Yanai's GCD measure} \insertCite{ramsay_1984_MatrixCorrelation}{repsim}. 
#' Another is \emph{Pillai's trace}, which is a popular statistic in MANOVA.
#' 
#' @param dlist a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.
#' @param ... optional parameters including\describe{
#' \item{centering}{a logical to apply centering (default:\code{TRUE}).}
#' \item{method}{type of summary statistics of the goodness of fit for CCA, either \code{"pillai"} or \code{"yanai"} (default: \code{"yanai"}).}
#' }
#' 
#' @return an \eqn{(N\times N)} matrix of pairwise similarity measures.
#' 
#' @examples 
#' \donttest{
#' ## load the 'twoclass' data
#' data("twoclass", package="repsim")
#' 
#' # compare Yanai's method (default) with and without centering
#' Y1 <- CCA(twoclass$data, centering=TRUE)
#' Y2 <- CCA(twoclass$data, centering=FALSE)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' fields::imagePlot(Y1, xaxt="n", yaxt="n", horizontal=TRUE, main="CCA-centering")
#' fields::imagePlot(Y2, xaxt="n", yaxt="n", horizontal=TRUE, main="CCA-no centering")
#' par(opar)
#' }
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @concept sim
#' @export
CCA <- function(dlist, ...){
  # ---------------------------------------------------------------
  # PREP
  # Inputs : explicit
  aux_check_configlist("CCA", dlist)
  
  # Inputs : implicit
  params = list(...)
  pnames = names(params)
  
  if ("centering"%in%pnames){
    par_centering = as.logical(params$centering)
  } else {
    par_centering = TRUE
  }
  if ("method"%in%pnames){
    par_measure = match.arg(params$method, c("yanai","pillai"))
  } else {
    par_measure = "yanai"
  }
  
  # ---------------------------------------------------------------
  # COMPUTE AND RETURN
  # prep
  N = length(dlist)
  output = array(0,c(N,N))
  
  # centering and compute orthonormal bases
  vec_Q = vector("list", length=N)
  for (n in 1:N){
    # Kornblith paper wants X(X^T X)^{-1/2} for extracting orthonormal bases.
    # vec_Q[[n]] = qr.Q(qr(as.matrix(base::scale(dlist[[n]], center=TRUE, scale=FALSE))))
    vec_Q[[n]] = src_orthobase(dlist[[n]], par_centering)
  }
  
  # iterate : diagonal
  for (n in 1:N){
    if (identical(par_measure,"yanai")){
      output[n,n] = CCA_yanai(vec_Q[[n]], vec_Q[[n]])
    } else {
      output[n,n] = CCA_pillai(vec_Q[[n]], vec_Q[[n]])
    }
  }
  
  # iterate : off-diagonal
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (identical(par_measure, "yanai")){
        output[i,j] <- CCA_yanai(vec_Q[[i]], vec_Q[[j]])
      } else {
        output[i,j] <- CCA_pillai(vec_Q[[i]], vec_Q[[j]])
      }
      output[j,i] <- output[i,j]
    }
  }
  
  # return
  return(output)
}

# auxiliary functions for CCA ---------------------------------------------
#' @keywords internal
#' @noRd
CCA_yanai <- function(Qx, Qy){
  p = min(base::ncol(Qx), base::ncol(Qy))
  return((base::norm(t(Qy)%*%Qx,type="F")^2)/p)
}
#' @keywords internal
#' @noRd
CCA_pillai <- function(Qx, Qy){
  p = min(base::ncol(Qx), base::ncol(Qy))
  svdQ = base::svd(t(Qy)%*%Qx)
  return(base::sum(svdQ$d)/p)
}