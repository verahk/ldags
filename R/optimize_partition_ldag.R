

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
                                    conf = as.matrix(expand.grid(levels)),
                                    verbose = FALSE){
  
  r <- ncol(counts)
  q <- nrow(counts)
  
  nlev <- lengths(levels)
  stride <- c(1, cumprod(nlev[-length(nlev)]))
  conf_enum <- seq_len(nrow(conf))-1 
  #conf_enum_par <- conf_enum - sweep(conf, 2, stride, "*")
  
  # compute seach part's contribution to the local score
  partition   <- unlist_partition(P)
  part_size   <- lengths(P)
  part_counts <- rowsum(counts, partition, reorder = T)
  part_scores <- famscore_bdeu_byrow(part_counts, ess, r, q, part_size)
  
  
  # search for label that improves score the most 
  keep_climb <- FALSE
  best_diff <- 0
  best_lab <- list()
  
  # check if all co-parent configs is already added to label
  #lens <- lengths(labels)
  #for (i in seq_along(nlev)[lens < (q/nlev-1)]) {
  for (i in seq_along(nlev)) { 
    # find candidate  labels /  co-parent configs to add to label
    # - enumerate these contexts by the rows in the CPT where node i is zero
    contexts <- conf_enum[conf[, i] == 0]
    
    # remove contexts already in label on edge from i
    contexts <- contexts[match(contexts, labels[[i]], nomatch = 0) == 0]
    
    # find the parts satisfied by each candidate label
    rows  <- outer(contexts, levels[[i]]*stride[i]+1, "+")
    parts <- split(partition[rows], .row(dim(rows)))
    #apply(rows, 1, function(x) cbind(conf[x, ], p = partition[x]), simplify = F)
    
    # add labels that are satisfied by current labeling
    redundant_labels <- vapply(parts, function(x) length(unique(x)), integer(1)) == 1
    
    if (any(redundant_labels)) {
      if (verbose) cat("\nAdd implicit labels: node", i, "contexts = ", paste(contexts[redundant_labels], collapse = ","))
      labels[[i]] <- append(labels[[i]], contexts[redundant_labels])
      if (all(redundant_labels)) next 
      parts <- parts[!redundant_labels]
      contexts <- contexts[!redundant_labels]
    }
    
    if (! (length(labels[[i]]) < q/nlev[i]-1) ) next
    
    parts <- lapply(parts, sort)
    dups  <- duplicated(parts)
    for (j in seq_along(parts)[!dups]) {
      
      collapse <- parts[[j]]
      
      # compare score of collapsed vs separated parts
      tmp_counts <- colSums(part_counts[collapse,])
      tmp_size   <- sum(part_size[collapse])
      tmp_score  <- famscore_bdeu_1row(tmp_counts, ess, r, q, tmp_size)
      diff  <-  tmp_score - sum(part_scores[collapse])
      
      if (diff > best_diff) {
        
        # check if regular 
        PP <- c(P[-collapse], list(unlist(P[collapse])))
        if (!all(is_regular(PP,  nlev, stride)))  next 
        
        if (any(dups)) {
          # find all labels implied by collapsing the two parts
          matches <- vapply(parts, function(x) all(x == collapse), logical(1))
          if (!(length(labels[[i]]) + sum(matches) < q/nlev[i])) {
            next 
          } else {
            context <- contexts[matches]
          }
        } else {
          context <- contexts[j]
        }
        
        best_partition <- PP
        best_diff <- diff
        best_lab  <- list(node = i, context = context, parts = collapse)
        keep_climb <- TRUE
      }
    }
  }
  

  
  if (keep_climb) {
  
    # update labels and partition
    labels[[best_lab$node]] <- c(labels[[best_lab$node]], best_lab$context)
    P <- best_partition 
    
    if (verbose) cat(sprintf("\ndiff = %s, node %s, context = %s, parts = %s, new partition: %s",
                             best_diff,
                             best_lab$node, 
                             paste(best_lab$context, collapse = ","),
                             paste(best_lab$parts, collapse = ","),
                             paste(unlist_partition(P), collapse = " ")))
   
    
    return(optimize_partition_ldag(counts, levels, ess, P, labels, conf, verbose))
    
  } else {
    return(list(partition = P,
                labels = labels,
                counts = part_counts, 
                scores = part_scores))
  }
}

