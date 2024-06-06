
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
    best_diff <- -(r-1)*lkappa
    nlev <- lengths(levels)
    stopifnot(all(nlev == 2))
    stride <- c(1, cumprod(nlev[-length(nlev)]))
    joint <- seq_len(prod(nlev))-1 
    
    for (i in seq.int(1, length(P)-1)) {
      for (j in seq.int(i+1, length(P))) {
        
        # check if collapsed part can be represented by labels
        new_part <- c(P[[i]], P[[j]])
        if (!is_ldag_consistent_part(new_part, nlev, stride, joint)) next 
        
        # candidate partition
        # PP  <- replace(P, i, list(new_part))
        # PP[j] <- NULL
        PP  <- c(P[-c(i, j)], list(new_part))
        
        # compute score of collapsed parts
        tmp_count <- part_counts[i,, drop = FALSE] + part_counts[j, , drop = FALSE]
        tmp_score <- famscore_bdeu_1row(tmp_count, ess, r, q, length(new_part))
        
        # compare score with score of non-collapsed regions
        diff <- tmp_score - (part_scores[i] + part_scores[j])
        
        # compare score with current best partition
        if (diff > best_diff) {
          
          # check if collapsing part i and j result in regular partition
          if (!all(is_regular(PP,  nlev, stride)))  next
          best_diff <- diff
          best_partition <- PP
          best_parts <- c(i, j)
          keep_climb <- TRUE
        }
      }
    }
  }
  
  if (keep_climb) {
    if (verbose) cat(sprintf("\ndiff = %s, parts = %s, new partition: %s",
                             best_diff,
                             paste(best_parts, collapse = ","),
                             paste(unlist_partition(best_partition), collapse = " ")))
    
    optimize_partition_part(counts, levels,  ess, best_partition,  lkappa, verbose)
  } else {
    return(list(partition = P, counts = part_counts, scores = part_scores))
  }
}
