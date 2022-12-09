# auxliary functions
# aux_check_configlist : valid list of configurations.


# CHECKERS ----------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_check_configlist <- function(f.name, config.list){
  # first, it should be a list
  if (!is.list(config.list)){
    stop(paste0("* ",f.name," : an input 'dlist' is not a list."))
  }
  # second, all should be matrices
  check_mat = unlist(lapply(config.list, base::is.matrix))
  if (any(check_mat==FALSE)){
    stop(paste0("* ",f.name," : the provided input 'dlist' contains a non-matrix element."))
  }
  # thrid, all should have equal number of rows
  check_nrow = unlist(lapply(config.list, base::nrow))
  if (length(unique(check_nrow))>1){
    stop(paste0("* ",f.name," : all elements in the 'dlist' must have same number of rows/observations."))
  }
}