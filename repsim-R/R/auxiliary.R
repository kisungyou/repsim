# aux_checker -------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_checker <- function(mats, fnname){
  # it should be a list
  if (!is.list(mats)){
    stop(paste0("* ", fnname, ": 'mats' should be a list of matrices"))
  }
  
  # iterate over mats
  for (it in seq_along(mats)){
    # set the current matrix
    now_mat = mats[[it]]
    
    # it should be a matrix
    if (!is.matrix(now_mat)){
      stop(paste0("* ", fnname, ": element #", it, " in 'mats' is not a matrix"))
    }
    
    # it should be numeric with finite values only
    if (!is.numeric(now_mat) || any(!is.finite(now_mat))){
      stop(paste0("* ", fnname, ": element #", it, " in 'mats' should be a numeric matrix with finite values only"))
    }
  }
}
