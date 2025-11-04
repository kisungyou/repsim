#' List of kernel functions
#' 
#' Returns a character vector of available kernel methods implemented in the package.
#' 
#' @return A character vector of kernel names.
#' 
#' @export
repsim_kernels <- function(){
  all_kernels = c("linear","rbf","rbf_mean","rbf_dualmed")
  return(all_kernels)
}