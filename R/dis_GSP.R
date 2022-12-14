#' Generalized Shape Metric invariant to Permutation
#' 
#' \insertCite{williams_2021_GeneralizedShapeMetrics;textual}{repsim} proposed 
#' a series of algorithms named as generalized shape metric. This function implements 
#' the \bold{permutation-invariant} metric. 
#' 
#' @param dlist a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.
#' @param ... optional parameters including\describe{
#' \item{centering}{a logical to apply centering (default:\code{TRUE}).}
#' \item{ndim}{an integer-valued dimension for the common feature space mapping (default: 2).}
#' }
#' 
#' @return an \eqn{(N\times N)} matrix of pairwise dissimilarity measures.
#' 
#' @examples 
#' \donttest{
#' ## load the 'twoclass' data
#' data("twoclass", package="repsim")
#' 
#' ## seperate data
#' data1 = twoclass$data
#' 
#' ## copy the data by shuffling columns
#' data2 = vector("list", length=length(data1))
#' for (i in 1:length(data1)){
#'   i_data = data1[[i]]
#'   i_idx  = sample(1:base::ncol(i_data), base::ncol(i_data), replace=FALSE)
#'   data2[[i]] = i_data[,i_idx]
#' }
#' 
#' ## compute GSP for the two data with default options
#' gsp1 <- GSP(data1)
#' gsp2 <- GSP(data2)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' fields::imagePlot(gsp1, xaxt="n", yaxt="n", horizontal=TRUE, main="GSP-data")
#' fields::imagePlot(gsp2, xaxt="n", yaxt="n", horizontal=TRUE, main="GSP-shuffled")
#' par(opar)
#' }
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @concept dis
#' @export
GSP <- function(dlist, ...){
  # ---------------------------------------------------------------
  # PREP
  # Inputs : explicit
  aux_check_configlist("GSP", dlist)

  # Inputs : implicit
  params = list(...)
  pnames = names(params)
  
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
    stop(paste0("* GSP : set the 'ndim' parameter to be in [2,",max_dim,"]."))
  }  
  
  # ---------------------------------------------------------------
  # COMPUTE AND RETURN
  # Common Space Mapping : (m x p x N) : N reps, m reference pts, p is ndim.
  common_data = aux_common_mapping(dlist, par_ndim, par_centering)
  
  # iterate
  N = length(dlist)
  output = array(0,c(N,N))
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      output[i,j] <- output[j,i] <- GSP_pairwise(common_data[,,i], common_data[,,j])
    }
  }
  
  # return
  return(output)
}


# auxiliary function ------------------------------------------------------
#' @keywords internal
#' @noRd
GSP_pairwise <- function(X, Y){
  cost_mat = t(Y)%*%X
  opt_Q    = lpSolve::lp.assign(cost_mat, direction="max")$solution
  return(base::norm(X-(Y%*%opt_Q), "F"))
}
