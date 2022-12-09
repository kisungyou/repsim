#' Centered Kernel Alignment
#' 
#' Given two centered representation matrices \eqn{X \in \mathbf{R}^{m\times p}} and \eqn{Y \in \mathbf{R}^{m\times q}}, 
#' the Centered Kernel Alignment (CKA) is defined using the HSIC as 
#' \deqn{\textrm{CKA}(X,Y) = \frac{\textrm{HSIC}(X,Y)}{\sqrt{
#' \textrm{HSIC}(X,X)\cdot \textrm{HSIC}(Y,Y)
#' }}}
#' where \eqn{K} and \eqn{L} are \eqn{(m\times m)} kernel matrices of \eqn{X} and \eqn{Y}, respectively and \eqn{H} is a centering matrix. 
#' Its definition is similar to the relationship between covariance and correlation where the latter is normalized version of the former. 
#' Currently, our implementation supports linear and RBF kernels. For a linear kernel, CKA is reduced to the RV coefficient. 
#' 
#' @param dlist a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.
#' @param ... optional parameters including\describe{
#' \item{kernel}{(case-insensitive) name of the kernel to be used; either \code{"linear"} or \code{"rbf"} (default: \code{"rbf"}).}
#' \item{sigma}{bandwidth parameter for the RBF kernel when \code{kernel="rbf"}. If \code{NULL}, it automatically selects 
#' a bandwidth as the median of pairwise distances. Otherwise, a number should be provided (default: \code{NULL}).
#' }
#' }
#' 
#' @return an \eqn{(N\times N)} matrix of pairwise similarity measures.
#' 
#' @examples 
#' \donttest{
#' ## load the 'twoclass' data
#' data("twoclass", package="repsim")
#' 
#' # compare linear & rbf kernels
#' CKAlin <- CKA(twoclass$data, kernel="linear")
#' CKArbf <- CKA(twoclass$data, kernel="rbf")
#' 
#' ## visualize
#' graphics.off()
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' fields::imagePlot(CKAlin, main="CKA-linear", xaxt="n", yaxt="n", horizontal=TRUE)
#' fields::imagePlot(CKArbf, main="CKA-rbf",    xaxt="n", yaxt="n", horizontal=TRUE)
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{cortes_2012_AlgorithmsLearningKernels}{repsim}
#' 
#' \insertRef{kornblith_2019_SimilarityNeuralNetwork}{repsim}
#' 
#' @concept sim
#' @export
CKA <- function(dlist, ...){
  # ---------------------------------------------------------------
  # PREP
  # Inputs : explicit
  aux_check_configlist("CKA", dlist)
  
  # Inputs : implicit
  params = list(...)
  pnames = names(params)
  
  if ("kernel"%in%pnames){
    par_kernel = params$kernel
    if (!is.character(par_kernel)){
      stop("* CKA : 'kernel' should be a character.")
    }
    par_kernel = match.arg(tolower(par_kernel), c("linear","rbf"))
  } else {
    par_kernel = "rbf"
  }
  if ("sigma"%in%pnames){
    if (is.null(params$sigma)&&(length(params$sigma)==1)){
      par_sigma = 1
      rbf_auto  = TRUE
    } else {
      par_sigma = as.double(params$sigma)
      rbf_auto  = FALSE
      if (par_sigma <= .Machine$double.eps){
        stop("* CKA : provide a positive number for the kernel bandwidth.")
      } 
    }
  } else {
    par_sigma = 1
    rbf_auto  = TRUE
  }
  
  # ---------------------------------------------------------------
  # COMPUTE & RETURN
  return(cpp_CKA(dlist, par_kernel, par_sigma, rbf_auto))
}
