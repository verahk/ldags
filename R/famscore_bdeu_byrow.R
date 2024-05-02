



#' Compute BDEU-score for each parent config
#'
#' @param x a matrix with observed counts
#' @param ess
#' @param r
#' @param q
#' @param s
#'  scaling factor for the hyperparam `alpha`
#'
#' @return an array with score for each parent config
#' @example 
#' n <- matrix(1:10, 10, 2)
#' famscore_bdeu_by_row(n, 1)
famscore_bdeu_byrow <- function(x, ess, r = ncol(x), q = nrow(x), s = 1) {
  alpha <- ess*s/(r*q) 
  lgamma(r*alpha) - r*lgamma(alpha) + rowSums(lgamma(alpha + x)) - lgamma(r*alpha + rowSums(x))
}
