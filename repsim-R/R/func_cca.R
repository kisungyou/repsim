#' Canonical Correlation Analysis
#' 
#' Compute pairwise CCA-based similarities between multiple representations,
#' summarized by either Yanai’s GCD measure \insertCite{ramsay_1984_MatrixCorrelation}{repsim} or Pillai’s trace statistic \insertCite{raghu_2017_SVCCASingularVector}{repsim}.
#'
#' @param mats A list of length \eqn{M} containing data matrices of size
#'   \eqn{(n_\mathrm{samples},\, p_k)}. All matrices must share the same number
#'   of rows for matching samples.
#' @param summary_type Character scalar indicating the CCA summary statistic.
#'   One of \code{"yanai"} or \code{"pillai"}. Defaults to \code{"yanai"} if
#'   \code{NULL}.
#'
#' @return An \eqn{M \times M} symmetric matrix of CCA summary similarities.
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
#' # compute two similarities
#' cca_gcd = cca(mats, summary_type="yanai")
#' cca_trace = cca(mats, summary_type="pillai")
#' 
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' labs <- paste0("rep ",1:20)
#' par(pty="s", mfrow=c(1,2))
#' 
#' image(cca_gcd[,20:1], axes=FALSE, main="CCA:Yanai's GCD")
#' axis(1, seq(0, 1, length.out=20), labels = labs, las = 2)
#' axis(2, at = seq(0, 1, length.out=20), labels = labs[20:1], las = 2)
#' 
#' image(cca_trace[,20:1], axes=FALSE, main="CCA:Pillai's Trace")
#' axis(1, seq(0, 1, length.out=20), labels = labs, las = 2)
#' axis(2, at = seq(0, 1, length.out=20), labels = labs[20:1], las = 2)
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{golub_1995_CanonicalCorrelationsMatrix}{repsim}
#' 
#' \insertAllCited{}
#' 
#' @export
cca <- function(mats, summary_type = NULL){
  # check the data
  aux_checker(mats, "cca")
  
  # check the summary type
  if ((is.null(summary_type))&&(length(summary_type) < 1)){
    par_summary = "yanai"
  } else {
    par_summary = match.arg(tolower(summary_type), c("yanai","pillai"))
  }
  
  # pass onto C++ part
  return(cpp_CCA(mats, par_summary))
}

# # my test
# set.seed(1)
# X = as.matrix(scale(as.matrix(iris[sample(1:150, 50, replace=FALSE),1:4])))
# Y = as.matrix(scale(as.matrix(USArrests)))
# n = nrow(X)
# p_X = ncol(X)
# p_Y = ncol(Y)
# this
# mats = vector("list", length=20L)
# for (i in 1:10){
#   mats[[i]] = X + matrix(rnorm(n*p_X, sd=1), nrow=n)
# }
# for (j in 11:20){
#   mats[[j]] = Y + matrix(rnorm(n*p_Y, sd=1), nrow=n)
# }
# 
# par(mfrow=c(1,3), pty="s")
# image(cca(mats, summary_type = "Pillai"), main="Pillai")
# image(cca(mats, summary_type = "Yanai"), main="Yanai")
# image(cca(mats), main="Default: Yanai")
#