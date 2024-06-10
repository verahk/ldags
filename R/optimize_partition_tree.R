
#' Optimize a local tree structure for a frequency table
#'
#' Optimize a regression-tree partitioning of the conditional outcome space. 
#' 
#' @rdname optimize_partition
#' @examples
#' 
#' # binary tree
#' levels <- list(0:1, 0:1)
#' counts <- cbind(c(10, 100, 100, 100), 1)
#' ess <- 1
#' res <- optimize_partition_tree(counts, levels, ess = ess, verbose = TRUE)
#' res$partition 
#' cbind(counts, get_parts(res$partition))
#' 
#' 
#' # mixed cardinality
#' levels <- list(0:1, 0:2)
#' counts <- cbind(10, 1, c(10, 10, 10**2, 10**2, 10**3, 10**3))
#' res <- optimize_partition_tree(counts, levels, ess = ess, verbose = T)
#' cbind(expand.grid(levels), n = counts, part = get_parts(res$partition))
#' 
#' # compare returned scores with score of partition
#' gr <- get_parts(res$partition)
#' tmp <- rowsum(counts, gr)
#' scores <- famscore_bdeu_byrow(tmp, ess, r = ncol(counts), q = nrow(counts), s = lengths(res$partition))
#' all.equal(scores, res$scores)
optimize_partition_tree <- function(counts, levels, ess, verbose = verbose) {
  

  nlev <- lengths(levels)
  stride <- c(1, cumprod(nlev[-length(nlev)]))
  conf  <- as.matrix(expand.grid(levels)) 
  r <- ncol(counts)
  q <- nrow(counts)
  
  score <- famscore_bdeu_byrow(matrix(colSums(counts), nrow = 1), r, q = 1, s = 1, ess)
  tree  <- grow_tree(counts, conf, score, stride, ess, r, q, verbose = verbose)
  
  # return partition and scores of each part / leaf
  list(partition = unname(list_leaves(tree, "outcomes")),
       scores = unname(unlist(list_leaves(tree, "scores"))))


}


list_leaves <- function(tree, name) {
  if (is.null(tree$branches)) {
    tree[name]
  } else {
    unlist(lapply(tree$branches, list_leaves, name = name), recursive = FALSE)
  }
}


grow_tree <- function(counts, conf, score, stride, ess, r = ncol(counts), q = nrow(counts), verbose = F) {
  
  
  keep_splitting <- FALSE
  if (nrow(counts) > 1) {
    # find best split
    best_split <- find_best_split(counts, conf, ess, r, q, best_score = score)
    keep_splitting <- length(best_split) > 0
  }
  if (keep_splitting) {
    # grow a new tree for each value of split variable
    if (verbose) cat(sprintf("\nSplitvariable: %s, Splitvalue: %s, Scores: %s",
                             best_split$var, best_split$vals, best_split$scores))
    
   
    best_split$branches <- vector("list", length(best_split$vals))
    for (k in best_split$vals+1) {
      indx <- conf[, best_split$var] == k-1
      best_split$branches[[k]] <- grow_tree(counts[indx, , drop = F], 
                                            conf[indx, , drop = F], 
                                            best_split$scores[k], 
                                            stride, ess, r, q, verbose)
    }
  } else {
    # return leaf node
    best_split <- list(outcomes = conf%*%stride,   # enumerate outcomes in part
                       scores = score)            # local-local score
  }
  return(best_split) 
}

find_best_split <- function(counts, conf, ess, r, q, best_score = 0) {
  best_split <- list()
  for (i in 1:ncol(conf)) {
    vals <- unique(conf[, i])
    if (length(vals) > 1) {
      tmp <- unname(rowsum(counts, conf[, i]))
      scores <- famscore_bdeu_byrow(tmp, ess, r, q, s = tabulate(conf[, i]+1))
      score  <- sum(scores)  
      if (score > best_score) {
        best_score <- score
        best_split$var <- i
        best_split$vals <- vals
        best_split$scores <- scores
      }
    }
  }
  return(best_split)
}
