# auxliary functions
# (01) aux_check_configlist -- : valid list of configurations.
# (02) aux_quick_pca --------- : apply PCA to the matrix
# (03) aux_common_mapping ---- : Williams et al. (2021) strategy for common space mapping
#                                PCA or zero padding after centering



# (01) aux_check_configlist -----------------------------------------------
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

# (02) aux_quick_pca ------------------------------------------------------
#' @keywords internal
#' @noRd
aux_quick_pca <- function(mat, ndim){
  mat_svds = RSpectra::svds(as.matrix(scale(mat, center = TRUE, scale=FALSE)), ndim)
  return(mat%*%mat_svds$v)
}

# (03) aux_common_mapping -------------------------------------------------
#' @keywords internal
#' @noRd
aux_common_mapping <- function(config_list, tgt_dim, centering=TRUE){
  N = length(config_list)
  m = base::nrow(config_list[[1]])
  
  output = array(0,c(m,tgt_dim,N))
  if (centering){
    for (n in 1:N){
      now_mat = as.matrix(base::scale(config_list[[n]], center=TRUE, scale=FALSE))
      now_dim = base::ncol(now_mat)
      if (now_dim < tgt_dim){
        output[,1:now_dim,n] = now_mat
      } else if (now_dim==tgt_dim){
        output[,,n] = now_mat
      } else {
        output[,,n] = aux_quick_pca(now_mat, tgt_dim)
      }
    }
  } else {
    for (n in 1:N){
      now_mat = config_list[[n]]
      now_dim = base::ncol(now_mat)
      if (now_dim < tgt_dim){
        output[,1:now_dim,n] = now_mat
      } else if (now_dim==tgt_dim){
        output[,,n] = now_mat
      } else {
        output[,,n] = aux_quick_pca(now_mat, tgt_dim)
      }
    }
  }
  return(output)
}

