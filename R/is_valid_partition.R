

#' Test a partition for CSI consistency and regularity
#'
#' CSI-consistency: No part in $P$ includes 
#' Regularity: All levels 
#'
#' @param P (list of integer vectors) 
#'  a partition of the outcome space defined by `levels`.
#'  enumerating the joint outcomes from `0` to `prod(lengths(nlev))-1`. 
#' @param levels (list of integer vectors)
#'  a list with outcomes of each variable.
#' @param check_regularity (logical constant)
#'  if TRUE, test if parition is regular, conditional on it beeing CSI-consistent.
#'  if FALSE, only CSI-consistency is tested.
#' @param verbose (logical constant)
#'  if TRUE
#' @return a logical constant 
#' @export
#'
#' @examples
#' nlev <- c(2, 2)
#' levels <- lapply(nlev-1, seq.int, from = 0)
#' 
#' # CSI-consistent partition
#' P <- list(c(0, 1, 3), 2)
#' is_valid_partition(P, levels, verbose = T) == TRUE
#' 
#' # not CSI-consistent partition
#' P <- list(c(0, 3), 2, 1)
#' is_valid_partition(P, levels, verbose = T) == FALSE
#' 
#' # CSI-consistent, but not regular
#' P <- list(c(0, 1), c(2, 3))
#' is_valid_partition(P, levels, verbose = T) == FALSE
#' P <- list(c(0, 2), c(1, 3))
#' is_valid_partition(P, levels, verbose = T) == FALSE
#' 
#' # Different cardinality
#' nlev <- c(2, 3, 4)
#' levels <- lapply(nlev-1, seq.int, from = 0)
#' P <- c(list(c(0, 2, 4)), list(c(1, 3, 5)), list(23:18), as.list(6:17))
#' is_valid_partition(P, levels, verbose = T)

is_valid_partition <- function(P, levels, nlev = lengths(levels), check_regularity = TRUE, verbose = FALSE) {
  n    <- length(levels)
  seqn <- seq_len(n)
  stride <- c(1, cumprod(nlev[-n]))
 
  # for tracking which partitions that includes all levels of each variable
  includes_all_lev <- matrix(FALSE, nrow = length(P), ncol = n)
  
  for (j in seq_along(P)[lengths(P) > 1]) {
    
    p <- P[[j]]
    K <- list()
    
    # configuration of nodes in the first row of current context
    vals <- (min(p)%/%stride)%%nlev
    
    for (i in which(vals == 0)) {
      
      # matrix where each row corresponds to a label of the edge i -> child
      labels <- bida:::expand_grid_fast(levels[-i], nlev[-i])
      
      
      # matrix where each row corresponds to equality constraints
      rows <- outer(c(labels%*%stride[-i]), levels[[i]]*stride[i], "+")
      
      # check which labels that is consistent with the part p
      indx <- apply(rows, 1, function(x) !any(match(x, p, 0L) == 0))
      if (any(indx)) {
        K[[i]] <- c(rows[indx, ])
        includes_all_lev[j, i] <- T
      }
    }
    
    # stop if current part can not be produced by labels
    K <- unique(unlist(K))
    if (! (length(K) == length(p) && all(sort(K) == sort(p))) ) {
      if (verbose) cat("P is not CSI-consistent\n")
      return(FALSE)
    }
  }
  
  # check regularity 
  if (all(lengths(P) <= min(nlev)) &&  any(colMeans(includes_all_lev) == 1)) {
    if (verbose) cat("P is not regular\n")
    return(FALSE)
  } else {
    return(TRUE)
  }
}



# test with 3 categories ----
if (FALSE) {
  nlev <- c(3, 3)
  levels <- lapply(nlev-1, seq.int, from = 0)
  P <- c(list(0:1), as.list(3:8))
  
  is_valid_partition(P, levels)
  
}
