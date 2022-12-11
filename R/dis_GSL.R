#' Generalized Shape Metric invariant to Linear Transform
#' 
#' \insertCite{williams_2021_GeneralizedShapeMetrics;textual}{repsim} proposed 
#' a series of algorithms named as generalized shape metric. This function implements 
#' the \bold{linear-invariant} metric. 
#' 
#' @param dlist a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.
#' @param ... optional parameters including\describe{
#' \item{alpha}{a regularization parameter in \eqn{[0,1]} (default: 0.5).}
#' \item{centering}{a logical to apply centering (default:\code{FALSE}).}
#' \item{ndim}{an integer-valued dimension for the common feature space mapping (default: 2).}
#' }
#' 
#' @return an \eqn{(N\times N)} matrix of pairwise similarity measures.
#' 
#' @examples 
#' \donttest{
#' ## load the 'twoclass' data
#' data("twoclass", package="repsim")
#' 
#' ## compare different levels of regularizations
#' alpha1 <- GSL(twoclass$data, alpha=0)
#' alpha2 <- GSL(twoclass$data, alpha=0.5)
#' alpha2 <- GSL(twoclass$data, alpha=1)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' fields::imagePlot(alpha1, xaxt="n", yaxt="n", horizontal=TRUE, main="GSL-alpha=0")
#' fields::imagePlot(alpha2, xaxt="n", yaxt="n", horizontal=TRUE, main="GSL-alpha=0.5")
#' fields::imagePlot(alpha2, xaxt="n", yaxt="n", horizontal=TRUE, main="GSL-alpha=1")
#' par(opar)
#' }
#' 
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @concept dis
#' @export
GSL <- function(dlist, ...){
  # ---------------------------------------------------------------
  # PREP
  # Inputs : explicit
  aux_check_configlist("GSL", dlist)
  
  # Inputs : implicit
  params = list(...)
  pnames = names(params)
  
  if ("ndim"%in%pnames){
    par_ndim = round(params$ndim)
  } else {
    par_ndim = 2
  }
  if ("alpha"%in%pnames){
    par_alpha = as.double(params$alpha)
  } else {
    par_alpha = 0.5
  }
  if ((par_alpha < 0)||(par_alpha > 1)){
    stop("* GSL : please set the 'alpha' parameter in [0,1].")
  }
  if ("centering"%in%pnames){
    par_centering = as.logical(params$centering)
  } else {
    par_centering = FALSE
  }
  
  max_dim = max(unlist(lapply(dlist, base::ncol)))
  if ((par_ndim > max_dim)||(par_ndim<2)){
    stop(paste0("* GSL : please set the 'ndim' parameter to be in [2,",max_dim,"]."))
  }
  
  # ---------------------------------------------------------------
  # COMPUTE AND RETURN
  # common space maping : (m x p x N) : N reps, m reference pts, p is ndim.
  common_data = aux_common_mapping(dlist, par_ndim, par_centering)
  
  # compute
  output = cpp_GSL(common_data, par_alpha)
  
  # return
  return(output)
}