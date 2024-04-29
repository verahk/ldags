

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
famscore_csi_from_data <- function(data, nlev, j, parentnodes, ess = 1, kappa = .1) {
  
  npa <- length(parentnodes) 
  if (npa == 0) {
    r <- nlev[j]
    counts <- tabulate(data[, j] +1, r)
    dim(counts) <- c(1, r)
    return(sum(famscore_bdeu_byrow(counts, ess, r, 1)))
  } else {
    
    cp  <- cumprod(nlev[parentnodes])
    stride <- c(1, cp)
    
    r <- nlev[j]
    q <- cp[npa]
    counts <- tabulate(data[, c(parentnodes, j)]%*%stride+1, q*r)
    dim(counts) <- c(q, r)
    
    if (FALSE && npa == 1) {
      # do not search for local CSIs if only one parent
      return(sum(famscore_bdeu_byrow(counts, ess, r, q)))
    } else {
      tmp <- optimize_local_csi_part(counts, nlev[parentnodes], ess)
      score <- sum(tmp$scores) + r*(q-length(tmp$scores))*log(kappa)
    }
  }
  return(score)
}



