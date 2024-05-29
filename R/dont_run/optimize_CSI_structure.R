
#' Optimize local CSI structure
#'
#' @param data (integer matrix)
#' @param nlev (integer vector)
#' @param j (integer)
#' @param parentnodes (integer vector)
#' @param ess (numeric constant)
#'  imaginary sample size for the parameter prior.
#' @param kappa (numeric constant)
#'  penalization factor for CSI-complexity. Search is terminated when the maximal 
#'  improvement in the score from fusing two parts is less than `-(r-1)*log(kappa)`, 
#'  where `r-1` is the reduction in the number of free parameters in the reduced CPT.
#' @return a list with the optimal partition and the associated scores
#' @export
#'
#' @examples
optimize_partition_from_data <- function(data, nlev, j, parentnodes, method = "LDAG", ess = 1, kappa = .1, verbose = FALSE) {
  
  r <- nlev[j]
  q <- prod(nlev[parentnodes])
  
  # compute marginal counts
  stride <- c(1, cumprod(nlev[parentnodes]))
  counts <- tabulate(data[, c(parentnodes, j)]%*%stride+1, q*r)
  dim(counts) <- c(q, r)
  
  # optimize CSI-structure
  if (method == "LDAG") {
    
    P <- as.list(seq.int(0, q-1)) # init 
    levels <- lapply(nlev[parentnodes]-1, seq.int, from = 0)
    
    optimize_CSI_structure_greedy(counts, levels, nlev[parentnodes], P, ess = ess, lkappa = log(kappa), verbose = verbose)
  } else if (method == "tree") {
    
  }
  
}

optimize_partition <- function(counts, lev, nlev = lengths(lev), params, verbose = F) {
  
  r <- nlev[j]
  q <- prod(nlev[parentnodes])
  
  # compute marginal counts
  stride <- c(1, cumprod(nlev[parentnodes]))
  counts <- tabulate(data[, c(parentnodes, j)]%*%stride+1, q*r)
  dim(counts) <- c(q, r)
  
  # optimize CSI-structure
  if (method == "LDAG") {
    
    P <- as.list(seq.int(0, q-1)) # init 
    levels <- lapply(nlev[parentnodes]-1, seq.int, from = 0)
    
    optimize_CSI_structure_greedy(counts, levels, nlev[parentnodes], P, ess = ess, lkappa = log(kappa), verbose = verbose)
  } else if (method == "tree") {
    
  }
  
}

