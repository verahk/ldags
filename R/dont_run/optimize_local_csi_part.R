
#' Compute optimal local CSI-structure of an local distribtuion
#'
#' @param counts (integer matrix)
#'  conditional frequency table 
#'  Each row correspond to one parent-config, each column to a outcome of the node.
#' @param nlev (integer vector)
#'  cardinality of parent variables
#' @param ess (integer)
#'  equivalent sample size
#'  
#' @return a list that contains: 
#' - `partitioing:` a vector of length `prod(nlev)` that defines a partitioning 
#'   of the parent-outcome space.
#' - `partsize`: number of outcomes in each region
#' - `scores`: the contributions of each region to the local score
#'
#' @examples
#' counts <- cbind(c(0, 5, 10, 100), 1)
#' counts/rowSums(counts)
#' nlev <- c(2, 2)
#' ess <- 1
#' optimize_csi_part(counts, nlev, ess)
#' 
#' 
#' counts <- structure(c(89L, 107L, 83L, 148L, 160L, 258L, 74L, 81L), dim = c(4L, 2L))
#' nlev <- c(2, 2)




optimize_local_csi_part <- function(counts, nlev, kappa, ess = 1, verbose = F) {
  q   <- nrow(counts)
  r   <- ncol(counts)
  
  # init
  S  <- uS <- letters[seq_len(q)] # seq_len(q)   # partition and list with regions
  nS <- rep(1, q)
  scores <- famscore_bdeu_byrow(counts, ess, r, q, s = nS)
  while (length(uS) > 2) {
    
    max_diff <- kappa
    for (i in seq.int(1, length(uS)-1)) {
      for (j in seq.int(i+1, length(uS))) {
        
        # ad hoc check for regularity:
        if (length(uS) < max(q/nlev +1)) {
          
        }
        tmp_counts <- counts[i,, drop = FALSE] + counts[j, , drop = FALSE]
        tmp_score <- famscore_bdeu_byrow(tmp_counts, ess, r, q, nS[i] + nS[j])
        diff <- tmp_score - (scores[i] + scores[j])
        if (diff > max_diff) {
          max_diff <- diff
          collapse <- c(i, j)
          score_collapsed  <- tmp_score
          counts_collapsed <- tmp_counts
        }
      }
    }
    if (max_diff > 0) {
      
      i <- collapse[1]
      j <- collapse[2]
      
      nS[i]       <- nS[i] + nS[j]
      counts[i, ] <- counts_collapsed
      scores[i]   <- score_collapsed
      
      indx <- S %in% uS[collapse]
      S[indx] <- uS[collapse][1]
      
      if (verbose) {
        cat("\nscore improvement:", max_diff, "collapse regions:", uS[collapse], "new partition:", S)
      }
      
      
      uS <- uS[-j]
      nS <- nS[-j]
      counts <- counts[-j, , drop = FALSE]
      scores  <- scores[-j]
      
    } else {
      break
    }
  }
  
  return(list(partitioning = S,
              partsize = nS,
              scores = scores))
}