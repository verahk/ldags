

#' Optimize a local LDAG structure for a frequency table
#' 
#' @rdname optimize_partition
#' @examples
#' # ldag
#' levels <- list(0:1, 0:1)
#' counts <- cbind(c(10, 20, 30, 40), rep(10, 4))
#' optimize_partition_ldag(counts, levels, ess = 1, verbose = T)
#' 
#' # mixed cardinality
#' levels <- list(0:1, 0:2)
#' counts <- cbind(10, 1, c(10, 10, 10**2, 10**2, 10**3, 10**3))
#' res <- optimize_partition_ldag(counts, levels, ess = ess, verbose = T)
#' cbind(expand.grid(levels), n = counts, part = unlist_partition(res$partition))
optimize_partition_ldag <- function(counts, levels, ess, 
                                    P = as.list(1:nrow(counts)-1), 
                                    labels = rep(list(integer()), length(levels)), 
                                    verbose = FALSE){
  
  r <- ncol(counts)
  q <- nrow(counts)
  
  nlev <- lengths(levels)
  stride <- c(1, cumprod(nlev[-length(nlev)]))

  # compute seach part's contribution to the local score
  partition   <- unlist_partition(P)
  part_size   <- tabulate(partition)
  part_counts <- rowsum(counts, partition, reorder = T)
  part_scores <- famscore_bdeu_byrow(part_counts, ess, r, q, part_size)
  
  # search for label that improves score the most 
  keep_climb <- FALSE
  best_diff <- 0
  for (i in seq_along(nlev)) {
    
    # find candidate contexts to add to label
    len <- length(labels[[i]])
    if (len == (q/nlev[[i]]-1)) { 
      # restrict labels to not include all co-parent contexts 
      next 
    } else {
      # enumerate contexts defined by the co-parents of node i
      contexts <- bida:::expand_grid_fast(levels[-i])%*%stride[-i]
    
      # remove contexts already in label on edge from i
      contexts <- contexts[match(contexts, labels[[i]], nomatch = 0) == 0]
    }

    add_to_rows <- 1+ levels[[i]]*stride[i] 
    for (context in contexts) {
      
      rows  <- context + add_to_rows    # rows in the CPT satisfied by label
      parts <- unique(partition[rows])  # corresponding parts in current partition
      
      # if context is consistent with one single region, add label 
      if (length(parts) == 1) {
        labels[[i]] <- append(labels[[i]], context)
        next
      }
      
      # compare score of collapsed vs separated parts
      tmp_counts <- matrix(colSums(part_counts[parts,]), nrow = 1)
      tmp_size   <- sum(part_size[parts])
      tmp_score  <- famscore_bdeu_byrow(tmp_counts, ess, r, q, tmp_size)
      diff  <-  tmp_score - sum(part_scores[parts])
      
      if (diff > best_diff) {
        best_diff <- diff
        best_lab  <- list(node = i, context = context, parts = parts)
        keep_climb <- TRUE
      }
    }
  }
    
  if (keep_climb) {
    if (verbose) cat("\nbest label:", unlist(best_lab)[1:2], "\n")
    
    # update labels and partition
    labels[[best_lab$node]] <- append(labels[[best_lab$node]], best_lab$context)
    P <- c(P[-best_lab$parts], list(unlist(P[best_lab$parts])))
    
    return(optimize_partition_ldag(counts, levels, ess, P, labels, verbose))
    
  } else {
    return(list(partition = P,
                labels = labels,
                counts = part_counts, 
                scores = part_scores))
  }
}
