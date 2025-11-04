#' Hilbert Schmidt Independence Criterion
#' 
#' Compute pairwise HSIC similarities between multiple representations using a
#' chosen kernel and estimator.
#' 
#' @param mats A list of length \eqn{M} containing data matrices of size 
#' \eqn{(n_\mathrm{samples},\, p_k)}. All matrices must share the same number of rows for matching samples.
#' @param kernel_type Character scalar indicating the kernel. Defaults to \code{"rbf"} (if \code{NULL}). 
#' See \code{\link{repsim_kernels}} for a list of available kernels.
#' @param estimator Character scalar indicating the HSIC estimator. Defaults to \code{"gretton"} (if \code{NULL}). 
#' See \code{\link{repsim_hsic}} for a list of available estimators.
#'
#' @return An \eqn{M \times M} matrix of HSIC values.
#' 
#' @examples
#' \donttest{
#' # --------------------------------------------------
#' # Use "iris" and "USArrests" datasets
#' #   1. apply scaling to reduce the effect of scales
#' #   2. add white noise to create multiple representations
#' #   3. generate 10 perturbations per each dataset
#' # --------------------------------------------------
#' # prepare the prototype
#' set.seed(1)
#' X = as.matrix(scale(as.matrix(iris[sample(1:150, 50, replace=FALSE),1:4])))
#' Y = as.matrix(scale(as.matrix(USArrests)))
#' n = nrow(X)
#' p_X = ncol(X)
#' p_Y = ncol(Y)
#' 
#' # generate 10 of each by perturbation
#' mats = vector("list", length=20L)
#' for (i in 1:10){
#'   mats[[i]] = X + matrix(rnorm(n*p_X, sd=1), nrow=n)
#' }
#' for (j in 11:20){
#'   mats[[j]] = Y + matrix(rnorm(n*p_Y, sd=1), nrow=n)
#' }
#' 
#' # compute similarity with rbf kernel and different estimators
#' hsic1 = hsic(mats, estimator="gretton")
#' hsic2 = hsic(mats, estimator="song")
#' hsic3 = hsic(mats, estimator="lange")
#' 
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' labs <- paste0("rep ",1:20)
#' par(mfrow=c(1,3), pty="s")
#' 
#' image(hsic1[,20:1], axes=FALSE, main="HSIC (Gretton)")
#' axis(1, seq(0, 1, length.out=20), labels = labs, las = 2)
#' axis(2, at = seq(0, 1, length.out=20), labels = labs[20:1], las = 2)
#' 
#' image(hsic2[,20:1], axes=FALSE, main="HSIC (Song)")
#' axis(1, seq(0, 1, length.out=20), labels = labs, las = 2)
#' axis(2, at = seq(0, 1, length.out=20), labels = labs[20:1], las = 2)
#' 
#' image(hsic3[,20:1], axes=FALSE, main="HSIC (Lange)")
#' axis(1, seq(0, 1, length.out=20), labels = labs, las = 2)
#' axis(2, at = seq(0, 1, length.out=20), labels = labs[20:1], las = 2)
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{gretton_2005_MeasuringStatisticalDependence}{repsim}
#' 
#' @export
hsic <- function(mats, kernel_type = NULL, estimator = NULL){
  # check the data
  aux_checker(mats, "hsic")
  
  # check the kernel type
  if ((is.null(kernel_type))||(length(kernel_type) < 1)){
    par_kernel = "rbf"
  } else {
    all_kernels = repsim_kernels()
    par_kernel  = tolower(kernel_type)
    if (!(par_kernel %in% all_kernels)){
      stop("* hsic : 'kernel_type' is invalid.")
    }
  }
  
  # check the estimator
  if ((is.null(estimator)||(length(estimator) < 1))){
    par_estimator = "gretton"
  } else {
    all_ests = repsim_hsic()
    par_estimator = tolower(estimator)
    if (!(par_estimator %in% all_ests)){
      stop("* hsic : 'estimator' is invalid.")
    }
  }
  
  # pass onto C++ part
  return(cpp_HSIC(mats, par_kernel, par_estimator))
}

# # my test
# set.seed(1)
# X = as.matrix(scale(as.matrix(iris[sample(1:150, 50, replace=FALSE),1:4])))
# Y = as.matrix(scale(as.matrix(USArrests)))
# n = nrow(X)
# p_X = ncol(X)
# p_Y = ncol(Y)
# 
# mats = vector("list", length=20L)
# for (i in 1:10){
#   mats[[i]] = X + matrix(rnorm(n*p_X, sd=1), nrow=n)
# }
# for (j in 11:20){
#   mats[[j]] = Y + matrix(rnorm(n*p_Y, sd=1), nrow=n)
# }
# 
# par(mfrow=c(2,2), pty="s")
# image(hsic(mats, kernel_type= "linear"), main="Linear")
# image(hsic(mats, kernel_type="rbf"), main="RBF")
# image(hsic(mats, kernel_type="rbf_mean"), main="RBF-mean")
# image(hsic(mats, kernel_type="rbf_dualmed"), main="RBF-dual")
