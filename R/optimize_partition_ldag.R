

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
    parts <- array(partition[rows], dim(rows))
    # check 
    # apply(rows, 1, function(x) cbind(conf[x, ], p = partition[x]), simplify = F)
    
    # find labels that implies no new independencies
    redundant_labels <- rowSums(parts == parts[, 1]) == ncol(parts)
    
    # add redundant labels to current set of labels
    if (any(redundant_labels)) {
      if (verbose) cat("\nAdd implicit labels: node", i, "contexts = ", paste(contexts[redundant_labels], collapse = ","))
      labels[[i]] <- append(labels[[i]], contexts[redundant_labels])
      if (all(redundant_labels)) next 
      parts <- parts[!redundant_labels, , drop = F]
      contexts <- contexts[!redundant_labels]
    }
    
    # if no additional label can be added to edge from i, go to next node
    if (! (length(labels[[i]]) < q/nlev[i]-1) ) next
    
   
    # find labels that implies the same independencies
    dups <- find_duplicated_rows(parts)
    #parts <- lapply(parts, sort)
    #dups  <- duplicated(parts)
    
    for (j in seq_along(contexts)[!dups]) {
      
      collapse <- parts[j, ]
      
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
          # candidates <- which(dups)[rowSums(parts[dups,, drop = F]) == sum(collapse)]
          # matches <- vector("logical", length(candidates))
          # for (jj in candidates) {
          #   matches[j] <- all(match(parts[jj, ], collapse, 0L) > 0)
          # } 
          
          matches <- apply(parts, 1, 
                           function(x) all(match(x, collapse, 0L) > 0))
          if (!(length(labels[[i]]) + sum(matches) < q/nlev[i])) {
            next 
          } else {
            #context <- contexts[c(j, candidates[matches])]
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


find_duplicated_rows <- function(x){
  n <- nrow(x)
  dups <- vector("logical", n)
  sums <- rowSums(x)
  if (anyDuplicated(sums) == 0) return(dups)
  seqn <- seq_len(n)
  for (i in seqn[-n]) {
    if (!dups[i]) {
      tab <- x[i, ]
      candidates <- seqn[replace((!dups & (sums == sums[i])), i, FALSE)]
      for (j in candidates) {
        dups[j] <- all(match(x[j, ], tab, 0L) > 0)
      }
    }
  }
  return(dups)
}


if (FALSE) {
  all_elem_equal <- function(y) vapply(y, function(x) all(x == x[1]), logical(1))
  all_elem_equal2 <- function(y) vapply(y, function(x) length(unique(x)), integer(1)) == 1
  all_elem_equal3 <- function(m) rowSums(m == m[, 1]) == 1
  
  tmp <- replicate(5, sample(1:3, 10, T), simplify = F)
  m <- do.call(rbind, tmp) 

  microbenchmark::microbenchmark(all_elem_equal(tmp),
                                 all_elem_equal2(tmp),
                                 all_elem_equal3(m))

  find_duplicated_rows <- function(x){
    n <- nrow(x)
    dups <- vector("logical", n)
    seqn <- seq_len(n)
    sums <- rowSums(x)
    
    for (i in seqn[-n]) {
      if (!dups[i]) {
        tab <- x[i, ]
        candidates <- seqn[replace((!dups & (sums == sums[i])), i, FALSE)]
        for (j in candidates) {
          dups[j] <- all(match(x[j, ], tab, 0L) > 0)
        }
      }
    }
    return(dups)
  }
  
  find_duplicated_rows2 <- function(x) {
    tmp <- lapply(asplit(x, 1), sort)
    duplicated(tmp)
  }
  
  
  tmp <- replicate(5, sample(1:10, 3, T), simplify = F)
  m <- do.call(rbind, tmp) 
  microbenchmark::microbenchmark(find_duplicated_rows(m),
                                 find_duplicated_rows2(m))
}

