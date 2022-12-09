#' Representational Similarity Metrics
#' 
#' not written yet
#' 
#' smaller value means closeness => metric => distance
#'                               => 0 means close
#'                               => 0 means similar
#'                               => 1 means dissimilar
#'                               => dissimilarity (metric/distance)
#'                               
#' larger  value means closeness => 0 means not close
#'                               => 0 means not similar 
#'                               => 0 means dissimilar
#'                               => 1 means similar
#'                               => similarity
#' 
#' 
#' similarity
#' 
#' @docType package
#' @noRd
#' @import Rdpack
#' @importFrom stats cov
#' @importFrom fields imagePlot
#' @importFrom Rcpp evalCpp
#' @useDynLib repsim
NULL
# pack <- "repsim"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))