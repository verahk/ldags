



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
#' 
#' dag <- rbind(Z = c(0, 1, 1), 
#'              X = c(0, 0, 1),
#'              Y = c(0, 0, 0))
#'
#' counts1 <- list()
#' counts1$Z <- round(rgamma(2, 10))
#' counts1$X <- matrix(round(rgamma(4, 10)), nrow = 2)
#' counts1$Y <- matrix(round(rgamma(6, 10)), nrow = 2)
#' 
#' s1 <- list(1, rep(1, 2), c(2, 1, 1))
#' coun


famscore_bdeu_1row <- function(x, ess, r = length(x), q = 1, s = 1) {
  ralpha <- ess*s/q
  alpha <- ralpha/r
  lgamma(ralpha) - r*lgamma(alpha) + sum(lgamma(alpha + x)) - lgamma(ralpha + sum(x))
}

famscore_bdeu_byrow <- function(x, ess, r = ncol(x), q = nrow(x), s = 1) {
  ralpha <- ess*s/q
  alpha <- ralpha/r
  lgamma(ralpha) - r*lgamma(alpha) + rowSums(lgamma(alpha + x)) - lgamma(ralpha + rowSums(x))
}
