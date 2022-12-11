#' Generalized Shape Metric invariant to Rotation
#' 
#' \insertCite{williams_2021_GeneralizedShapeMetrics;textual}{repsim} proposed 
#' a series of algorithms named as generalized shape metric. This function implements 
#' the \bold{rotation-invariant} metric. 
#' 
#' @param dlist a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.
#' @param ... optional parameters including\describe{
#' \item{centering}{a logical to apply centering (default:\code{TRUE}).}
#' \item{ndim}{an integer-valued dimension for the common feature space mapping (default: 2).}
#' \item{use.angular}{a logical; \code{FALSE} returns a Procrustes size-and-shape distance with reflections \eqn{d_1} 
#' and \code{TRUE} to return Kendall's shape space distance \eqn{\theta_1} (default: \code{FALSE}).}
#' }
#' 
#' @return an \eqn{(N\times N)} matrix of pairwise similarity measures.
#' 
#' @examples 
#' \donttest{
#' ## load the 'twoclass' data
#' data("twoclass", package="repsim")
#' 
#' ## compare Kendall and Procrustes distances
#' gsrK <- GSR(twoclass$data, use.angular=TRUE)
#' gsrP <- GSR(twoclass$data, use.angular=FALSE)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' fields::imagePlot(gsrK, xaxt="n", yaxt="n", horizontal=TRUE, main="GSR-Kendall")
#' fields::imagePlot(gsrP, xaxt="n", yaxt="n", horizontal=TRUE, main="GSR-Procrustes")
#' par(opar)
#' }
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @concept dis
#' @export
GSR <- function(dlist, ...){
  # ---------------------------------------------------------------
  # PREP
  # Inputs : explicit
  aux_check_configlist("GSR", dlist)
  
  # Inputs : implicit
  params = list(...)
  pnames = names(params)
  
  if ("use.angular"%in%pnames){
    par_angular = as.logical(params$use.angular)
  } else {
    par_angular = FALSE
  }
  if ("centering"%in%pnames){
    par_centering = as.logical(params$centering)
  } else {
    par_centering = TRUE
  }
  if ("ndim"%in%pnames){
    par_ndim = round(params$ndim)
  } else {
    par_ndim = 2
  }
  max_dim = max(unlist(lapply(dlist, base::ncol)))
  if ((par_ndim > max_dim)||(par_ndim<2)){
    stop(paste0("* GSR : set the 'ndim' parameter to be in [2,",max_dim,"]."))
  }
  
  # ---------------------------------------------------------------
  # COMPUTE AND RETURN
  # common space maping : (m x p x N) : N reps, m reference pts, p is ndim.
  common_data = aux_common_mapping(dlist, par_ndim, par_centering)
  
  # branching
  if (par_angular){
    output = cpp_GSR_kendall(common_data)
  } else {
    output = cpp_GSR_procrustes(common_data)
  }
  
  # return
  return(output)
}
