#' Singular Vector Canonical Correlation Analysis 
#' 
#' Since a standard version of Canonical Correlation Analysis (CCA) often suffers from 
#' perturbation, denoising based on singular value decomposition is applied as a preprocessing step. 
#' This is coined as the method of Singular Vector CCA (SVCCA), which measures similarity 
#' between representation matrices \eqn{[0,1]}. Denoising is applied by the parameter \code{proportion}, 
#' which truncates the singular vector components by the amount of variance explained.
#' 
#' @param dlist a length-\eqn{N} list of representation matrices. Note that all matrices must have same number of rows that are assumed to be matched observations across different representations.
#' @param ... optional parameters including\describe{
#' \item{centering}{a logical to apply centering (default:\code{TRUE}).}
#' \item{method}{type of summary statistics of the goodness of fit for CCA, either \code{"pillai"} or \code{"yanai"} (default: \code{"pillai"}).}
#' \item{proportion}{threshold of the amount of variance explained for denoising (default: \eqn{0.95}).}
#' }
#' 
#' @return an \eqn{(N\times N)} matrix of pairwise similarity measures.
#' 
#' @examples 
#' \donttest{
#' ## load the 'twoclass' data
#' data("twoclass", package="repsim")
#' 
#' # compare CCA and SVCCA with two levels of threshold
#' runCC <- CCA(twoclass$data)
#' runS1 <- SVCCA(twoclass$data, proportion=0.50)
#' runS2 <- SVCCA(twoclass$data, proportion=0.95)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' fields::imagePlot(runCC, xaxt="n", yaxt="n", horizontal=TRUE, main="CCA")
#' fields::imagePlot(runS1, xaxt="n", yaxt="n", horizontal=TRUE, main="SVCCA thr=50%")
#' fields::imagePlot(runS2, xaxt="n", yaxt="n", horizontal=TRUE, main="SVCCA thr=95%")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{raghu_2017_SVCCASingularVector}{repsim}
#' 
#' @concept sim
#' @export
SVCCA <- function(dlist, ...){
  # ---------------------------------------------------------------
  # PREP
  # Inputs : explicit
  aux_check_configlist("SVCCA", dlist)
  
  # Inputs : implicit
  params = list(...)
  pnames = names(params)
  
  if ("centering"%in%pnames){
    par_center = as.logical(params$centering)
  } else {
    par_center = TRUE
  }
  if ("proportion"%in%pnames){
    par_ratio = as.double(params$proportion)
    par_eps   = 1e-6
    if ((par_ratio <= par_eps)||(par_ratio >= (1-par_eps))){
      stop("* SVCCA : 'proportion' parameter should be in (0,1) as it is ratio of variance explained.")
    }
  } else {
    par_ratio = 0.95
  }
  if ("method"%in%pnames){
    tmp_method = tolower(as.character(params$method))
    par_method = match.arg(tmp_method, c("yanai","pillai"))
  } else {
    par_method = "pillai"
  }
  
  # ---------------------------------------------------------------
  # COMPUTE
  # preps
  N = length(dlist)
  output = array(0,c(N,N))
  vec_Q = vector("list", length=N)
  for (n in 1:N){
    vec_Q[[n]] = SVCCA_extractQ(dlist[[n]], par_center, par_ratio)
  }
  
  # iterate : diagonal ------------------------------- use mean : pillai
  for (n in 1:N){
    if (identical(par_method, "yanai")){
      output[n,n] = CCA_yanai(vec_Q[[n]], vec_Q[[n]])
    } else {
      output[n,n] = CCA_pillai(vec_Q[[n]], vec_Q[[n]])      
    }
  }
  
  # iterate : off-diagonal
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (identical(par_method, "pillai")){
        output[i,j] <- output[j,i] <- CCA_pillai(vec_Q[[i]], vec_Q[[j]])  
      } else {
        output[i,j] <- output[j,i] <- CCA_yanai(vec_Q[[i]], vec_Q[[j]])  
      }
    }
  }
  
  # return
  return(output)
}



# auxiliary for SVCCA -----------------------------------------------------
#' @keywords internal
#' @noRd
SVCCA_extractQ <- function(X, centering, threshold){
  # run SVD + (contingent) centering
  if (centering){
    svdX = base::svd(as.matrix(base::scale(X, center=TRUE, scale=FALSE)))
  } else {
    svdX = base::svd(X)
  }
  
  # thresholding
  svdd = svdX$d
  csum = base::cumsum(svdd^2)/base::sum(svdd^2)
  ncom = max(1, min(which(csum>=threshold)))
  
  # reconstruction
  if (ncom < 2){
    recon = base::outer(as.vector(svdX$u[,1]), as.vector(svdX$v[,1]))*(svdX$d[1])
  } else {
    recon = svdX$u[,1:ncom]%*%diag(svdX$d[1:ncom])%*%t(svdX$v[,1:ncom])
  }
  
  # compute Q and return
  matQ = src_orthobase(recon, FALSE)
  return(matQ)
}
