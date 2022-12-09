#' Dot Product-Based Similarity
#' 
#' Given two centered representation matrices \eqn{X \in \mathbf{R}^{m\times p}} and \eqn{Y \in \mathbf{R}^{m\times q}}, 
#' the dot product-based similarity is measured as
#' \deqn{s(X,Y) = \|X^\top Y \|_F^2,}
#' which measures cross-variable similarities. 
#' 
#' @param dlist a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.
#' @param ... optional parameters including\describe{
#' \item{centering}{a logical to apply centering (default:\code{TRUE}).}
#' }
#' 
#' @return an \eqn{(N\times N)} matrix of pairwise similarity measures.
#' 
#' @examples
#' \donttest{
#' ## load the 'twoclass' data
#' data("twoclass", package="repsim")
#' 
#' ## compute the similarity with and without centering
#' DP1 <- DP(twoclass$data, centering=TRUE)
#' DP2 <- DP(twoclass$data, centering=FALSE)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' fields::imagePlot(DP1, xaxt="n", yaxt="n", horizontal=TRUE, main="DP-centering")
#' fields::imagePlot(DP2, xaxt="n", yaxt="n", horizontal=TRUE, main="DP-no centering")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{kornblith_2019_SimilarityNeuralNetwork}{repsim}
#' 
#' @concept sim
#' @export
DP <- function(dlist, ...){
  # ---------------------------------------------------------------
  # PREP
  # Inputs : explicit
  aux_check_configlist("DP", dlist)
  
  # Inputs : implicit
  params = list(...)
  pnames = names(params)
  
  if ("centering"%in%pnames){
    par_centering = as.logical(params$centering)
  } else {
    par_centering = TRUE
  }
  
  # ---------------------------------------------------------------
  # COMPUTE AND RETURN
  return(cpp_DP(dlist, par_centering))
}