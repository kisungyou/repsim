#' Singular Vector Canonical Correlation Analysis
#' 
#' Compute pairwise singular vector CCA (SVCCA) similarities between multiple representations.
#' SVCCA first mean-centers and denoises each representation via SVD, 
#' retaining components explaining a high fraction of variance at 99% threshold. Then, 
#' CCA is applied to the reduced representations, and the similarity is summarized 
#' with either Yanai’s GCD or Pillai’s trace.
#' 
#' @param mats A list of length \eqn{M} containing data matrices of size
#'   \eqn{(n_\mathrm{samples},\, p_k)}. All matrices must share the same number
#'   of rows for matching samples.
#' @param summary_type Character scalar indicating the CCA summary statistic.
#'   One of \code{"yanai"} or \code{"pillai"}. Defaults to \code{"yanai"} if
#'   \code{NULL}.
#'
#' @return An \eqn{M \times M} symmetric matrix of SVCCA summary similarities.
#' 
#' @references 
#' \insertRef{raghu_2017_SVCCASingularVector}{repsim}
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
#' svcca_gcd = svcca(mats, summary_type="yanai")
#' svcca_trace = svcca(mats, summary_type="pillai")
#' 
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' labs <- paste0("rep ",1:20)
#' par(pty="s", mfrow=c(1,2))
#' 
#' image(svcca_gcd[,20:1], axes=FALSE, main="SVCCA:Yanai's GCD")
#' axis(1, seq(0, 1, length.out=20), labels = labs, las = 2)
#' axis(2, at = seq(0, 1, length.out=20), labels = labs[20:1], las = 2)
#' 
#' image(svcca_trace[,20:1], axes=FALSE, main="SVCCA:Pillai's Trace")
#' axis(1, seq(0, 1, length.out=20), labels = labs, las = 2)
#' axis(2, at = seq(0, 1, length.out=20), labels = labs[20:1], las = 2)
#' par(opar)
#' }
#' 
#' @seealso \code{\link{cca}}
#' 
#' @export
svcca <- function(mats, summary_type = NULL){
  # check the data
  aux_checker(mats, "svcca")
  
  # check the summary type
  if ((is.null(summary_type))&&(length(summary_type) < 1)){
    par_summary = "yanai"
  } else {
    par_summary = match.arg(tolower(summary_type), c("yanai","pillai"))
  }
  
  # pass onto C++ part
  return(cpp_SVCCA(mats, par_summary))
}