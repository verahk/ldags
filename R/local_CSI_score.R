

#' Compute local score in a labeled DAG from data
#' 
#' Compute sample counts over the relevant subset of variables and compute score.
#' If more than one parent return the score of the optimal local CSI-structure.

#' @rdname 
#' @param data 
#' @param nlev 
#' @param j 
#' @param parentnodes 
#' @param ess 
#' @param kappa 
#'
#' @return score of node `j` given `parentnodes` and optimal CSI-structure 
#' @export
#'
#' @examples
local_CSI_score <- function(data, nlev, j, parentnodes, ess, kappa) {
  npar <- length(parentnodes)
  if (npar == 0) {
    r <- nlev[j]
    counts <- matrix(tabulate(data[, j] +1, nlev[j]), ncol = nlev[j])
    score  <- famscore_bdeu_byrow(counts, ess, r, 1)
    new_local_CSI_score(j, parentnodes, score, counts)
  } else if (npar == 1) {
    r <- nlev[j]
    q <- nlev[parentnodes]
    counts <- tabulate(data[, c(parentnodes, j)]%*%c(1, q) +1, q*r)
    dim(counts) <- c(q, r)
    score <- sum(famscore_bdeu_byrow(counts, ess, r, q))
    new_local_CSI_score(j, parentnodes, score, counts)
  } else {
    
    # optimal CSI-structure
    res    <- optimize_CSI_structure(data, nlev, j, parentnodes, ess = ess, kappa = kappa)
    
    # compute penalized score
    r <- nlev[j]
    q <- prod(nlev[parentnodes])
    score <- sum(res$scores) + (q-length(res$partition))*(r-1)*log(kappa)
    
    # return score object
    new_local_CSI_score(j, parentnodes, score, res$counts, res$partition)
  }
}

#' @rdname local_CSI_score
#' @param score 
#' @param counts 
#' @param partition 
#' @return a 
new_local_CSI_score <- function(j, parentnodes, score, counts, partition = NULL) {
  list(node = j,
       parentnodes = parentnodes, 
       score = score, 
       counts = counts,
       partition = partition)
}