

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
#' res <- optimize_partition_ldag(counts, levels, ess = 1, verbose = T)
#' cbind(expand.grid(levels), n = counts, part = unlist_partition(res$partition))
optimize_partition_ldag <- function(counts, levels, ess, 
                                    P = as.list(1:nrow(counts)-1), 
                                    labels = rep(list(integer()), length(levels)), 
                                    conf = bida:::expand_grid_fast(levels, nlev),
                                    verbose = FALSE){
  
  r <- ncol(counts)
  q <- nrow(counts)
  
  nlev <- lengths(levels)
  stride <- c(1, cumprod(nlev[-length(nlev)]))
  conf_enum <- seq_len(nrow(conf))-1 
  #conf_enum_par <- conf_enum - sweep(conf, 2, stride, "*")
  
  # compute seach part's contribution to the local score
  partition   <- unlist_partition(P)
  part_size   <- tabulate(partition)
  part_counts <- rowsum(counts, partition, reorder = T)
  part_scores <- famscore_bdeu_byrow(part_counts, ess, r, q, part_size)
  
  # search for label that improves score the most 
  keep_climb <- FALSE
  best_diff <- 0
  for (i in seq_along(nlev)) {
    
    # check if all co-parent configs is already added to label
    len <- length(labels[[i]])
    if (len >= (q/nlev[[i]]-1)) { 
      next 
    }
    
    # find candidate  labels /  co-parent configs to add to label
    # - enumerate these contexts by the rows in the CPT where node i is zero
    contexts <- conf_enum[conf[, i] == 0]
    
    # remove contexts already in label on edge from i
    contexts <- contexts[match(contexts, labels[[i]], nomatch = 0) == 0]
    
    # find the parts satisfied by each candidate label
    rows  <- outer(contexts, levels[[i]]*stride[i]+1, "+")
    parts <- partition[rows] # each row is the parts collapsed by corresp label
    dim(parts) <- dim(rows)
    
    for (j in 1:nrow(parts)) {
      
      collapse <- parts[j, ]
      if (anyDuplicated(collapse) > 0) {
        collapse <- unique(collapse)
        # if context is consistent with one single region, add label 
        if (length(collapse) == 1) {
          labels[[i]] <- append(labels[[i]], contexts[j])
          next
        }
      } 
      
      # compare score of collapsed vs separated parts
      tmp_counts <- colSums(part_counts[collapse,])
      tmp_size   <- sum(part_size[collapse])
      tmp_score  <- famscore_bdeu_1row(tmp_counts, ess, r, q, tmp_size)
      diff  <-  tmp_score - sum(part_scores[collapse])
      
      if (diff > best_diff) {
        best_diff <- diff
        best_lab  <- list(node = i, context = contexts[j], parts = collapse)
        keep_climb <- TRUE
      }
    }
  }
    
  if (keep_climb) {
    if (verbose) cat("\nbest label:", unlist(best_lab), "\n")
    
    # update labels and partition
    labels[[best_lab$node]] <- append(labels[[best_lab$node]], best_lab$context)
    P <- c(P[-best_lab$parts], list(unlist(P[best_lab$parts])))
    
    return(optimize_partition_ldag(counts, levels, ess, P, labels, conf, verbose))
    
  } else {
    return(list(partition = P,
                labels = labels,
                counts = part_counts, 
                scores = part_scores))
  }
}
