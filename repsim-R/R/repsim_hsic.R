#' List of HSIC estimators
#' 
#' Returns a character vector of available HSIC estimators implemented in 
#' the package.
#' 
#' @return A character vector of estimator names.
#' 
#' @references 
#' \insertRef{gretton_2005_MeasuringStatisticalDependence}{repsim}
#' 
#' \insertRef{song_2007_SupervisedFeatureSelection}{repsim}
#' 
#' \insertRef{lange_2023_DeepNetworksPaths}{repsim}
#' 
#' @export
repsim_hsic <- function(){
  all_ests = c("gretton", "song", "lange")
  return(all_ests)
}