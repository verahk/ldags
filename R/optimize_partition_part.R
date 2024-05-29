
#' Title
#'
#' @param counts 
#' @param levels 
#' @param nlev 
#' @param P 
#' @param ess 
#' @param lkappa 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
#' levels <- list(0:1, 0:1)
#' counts <- cbind(c(10, 100, 100, 100), rep(10, 4))
#' optimize_partition_part(counts, levels, ess = 1, verbose = T)
optimize_partition_part <- function(counts, levels,  ess, P = as.list(1:nrow(counts)-1), lkappa = 0, verbose = FALSE){
  
  q   <- nrow(counts)
  r   <- ncol(counts)

  # compute seach part's contribution to the local score
  part_size   <- lengths(P)
  part_counts <- rowsum(counts, unlist_partition(P, part_size), reorder = T)
  part_scores <- famscore_bdeu_byrow(part_counts, ess, r, q, part_size)
  
  # find the two current parts that increase the score most if collapsed
  keep_climb <- FALSE
  if (length(P) > 2) {
    max_diff <- -(r-1)*lkappa
    nlev <- lengths(levels)
    stopifnot(all(nlev == 2))
    stride <- c(1, cumprod(nlev[-length(nlev)]))
    
    for (i in seq.int(1, length(P)-1)) {
      for (j in seq.int(i+1, length(P))) {
        
        # candidate partition
        PP  <- replace(P, i, list(c(P[[i]], P[[j]])))
        PP[j] <- NULL
        
        # compute score of collapsed parts
        tmp_count <- part_counts[i,, drop = FALSE] + part_counts[j, , drop = FALSE]
        tmp_score <- famscore_bdeu_byrow(tmp_count, ess, r, q, length(PP[[i]]))
        
        # compare score with score of non-collapsed regions
        diff <- tmp_score - (part_scores[i] + part_scores[j])
        
        # compare score with current best partition
        if (diff > max_diff) {
          
          # check if collapsing part i and j result in regular, CSI-consistent partition
          if (!is_valid_partition(PP, levels, nlev, stride, verbose = F))  next
          
          if (verbose) cat("\n (i, j):", c(i, j), "diff:", diff, "partition:", unlist_partition(PP))
          max_diff <- diff
          best_partition <- PP
          keep_climb <- TRUE
        }
      }
    }
  }
  
  if (keep_climb) {
    if (verbose) cat("\nbest partition:", unlist_partition(best_partition))
    optimize_partition_part(counts, levels,  ess, best_partition,  lkappa, verbose)
  } else {
    return(list(partition = P, counts = part_counts, scores = part_scores))
  }
}
