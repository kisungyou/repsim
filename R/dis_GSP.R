#' Generalized Shape Metric invariant to Permutation
#' 
#' @param dlist a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.
#' @param ... optional parameters including\describe{
#' \item{centering}{a logical to apply centering (default:\code{FALSE}).}
#' \item{ndim}{an integer-valued dimension for the common feature space mapping (default: 2).}
#' }
#' 
#' @references 
#' \insertRef{williams_2021_GeneralizedShapeMetrics}{repsim}
#' 
#' @concept dis
#' @export
GSP <- function(dlist){
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
    par_centering = FALSE
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
  # COMPUTATION
}