#' Linear Regression-Based Similarity
#' 
#' Given two centered representation matrices \eqn{X \in \mathbf{R}^{m\times p}} and \eqn{Y \in \mathbf{R}^{m\times q}}, 
#' the goodness-of-fit of multivariate linear regression is reported using the \eqn{R^2} statistic. 
#' This measure is asymmetric and bounded in \eqn{[0,1]}. 
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
#' LR1 <- LR(twoclass$data, centering=TRUE)
#' LR2 <- LR(twoclass$data, centering=FALSE)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' fields::imagePlot(LR1, xaxt="n", yaxt="n", horizontal=TRUE, main="LR-centering")
#' fields::imagePlot(LR2, xaxt="n", yaxt="n", horizontal=TRUE, main="LR-no centering")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{kornblith_2019_SimilarityNeuralNetwork}{repsim}
#' 
#' @concept sim
#' @export
LR <- function(dlist, ...){
  # ---------------------------------------------------------------
  # PREP
  # Inputs : explicit
  aux_check_configlist("LR", dlist)

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
  return(cpp_LR(dlist, par_centering))
}
