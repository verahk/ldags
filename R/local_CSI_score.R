

#' Compute local score in a labeled DAG from data
#' 
#' Compute sample counts over the relevant subset of variables and compute score.
#' If more than one parent return the score of the optimal local CSI-structure.
#' 
#' @param j node to be scored
#' @param parentsnodes parents of this node
#' @param n number of nodes
#' @param param object of `scoreparameters`
#'
#' @return score of node `j` given `parentnodes` and optimal CSI-structure 
local_CSI_score <- function(data, nlev, j, parentnodes, ess = 1, kappa = .1) {
  
  npar <- length(parentnodes) 
  if (npar == 0) {
    r <- nlev[j]
    counts <- matrix(tabulate(data[, j] +1, nlev[j]), ncol = nlev[j])
    famscore_bdeu_byrow(counts, ess, r, 1)
  } else if (npar == 1) {
    r <- nlev[j]
    q <- nlev[parentnodes]
    counts <- tabulate(data[, c(parentnodes, j)]%*%c(1, q) +1, q*r)
    dim(counts) <- c(q, r)
    sum(famscore_bdeu_byrow(counts, ess, r, q))
  } else {

    # optimal CSI-structure
    res <- optimize_CSI_structure(data, nlev, j, parentnodes, ess = ess, kappa = kappa)
    
    # return penalized score
    r <- nlev[j]
    q <- prod(nlev[parentnodes])
    sum(res$scores) + (q-length(res$partition))*(r-1)*log(kappa)
  }
}

