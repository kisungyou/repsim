#' Hilbert-Schmidt Independence Criterion
#' 
#' Given two centered representation matrices \eqn{X \in \mathbf{R}^{m\times p}} and \eqn{Y \in \mathbf{R}^{m\times q}}, 
#' the Hilbert-Schmidt Independence Criterion (HSIC) is defined as 
#' \deqn{\textrm{HSIC}(X,Y) = \frac{1}{(m-1)^2} \textrm{tr}(KHLH),}
#' where \eqn{K} and \eqn{L} are \eqn{(m\times m)} kernel matrices of \eqn{X} and \eqn{Y}, respectively and \eqn{H} is a centering matrix. 
#' Currently, our implementation supports linear and RBF kernels. 
#' 
#' @param dlist a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.
#' @param ... optional parameters including\describe{
#' \item{kernel}{(case-insensitive) name of the kernel to be used; either \code{"linear"} or \code{"rbf"} (default: \code{"rbf"}).}
#' \item{sigma}{bandwidth parameter for the RBF kernel when \code{kernel="rbf"}. If \code{NULL}, it automatically selects 
#' a bandwidth as the median of pairwise distances. Otherwise, a number should be provided (default: \code{NULL}).}
#' }
#' 
#' @return an \eqn{(N\times N)} matrix of pairwise similarity measures.
#' 
#' @examples 
#' \donttest{
#' ## load the 'twoclass' data
#' data("twoclass", package="repsim")
#' 
#' ## compare linear & rbf kernels
#' HSIClin <- HSIC(twoclass$data, kernel="linear")
#' HSICrbf <- HSIC(twoclass$data, kernel="rbf")
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' fields::imagePlot(HSIClin, main="HSIC-linear", xaxt="n", yaxt="n", horizontal=TRUE)
#' fields::imagePlot(HSICrbf, main="HSIC-rbf",    xaxt="n", yaxt="n", horizontal=TRUE)
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{gretton_2005_MeasuringStatisticalDependence}{repsim}
#' 
#' @concept sim
#' @export
HSIC <- function(dlist, ...){
  # ---------------------------------------------------------------
  # PREP
  # Inputs : explicit
  aux_check_configlist("HSIC", dlist)
  
  # Inputs : implicit
  params = list(...)
  pnames = names(params)
  
  if ("kernel"%in%pnames){
    par_kernel = params$kernel
    if (!is.character(par_kernel)){
      stop("* HSIC : 'kernel' should be a character.")
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
        stop("* HSIC : provide a positive number for the kernel bandwidth.")
      } 
    }
  } else {
    par_sigma = 1
    rbf_auto  = TRUE
  }

  # ---------------------------------------------------------------
  # COMPUTE & RETURN
  return(cpp_HSIC(dlist, par_kernel, par_sigma, rbf_auto))
}