

#' Test a partition for CSI consistency  and regularity
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
#' is_valid_partition(P, levels, verbose = T)
#' 
#' # not CSI-consistent partition
#' P <- list(c(0, 3), 2, 1)
#' is_valid_partition(P, levels, verbose = T)
#' 
#' # CSI-consistent, but not regular
#' P <- list(c(0, 1), c(2, 3))
#' is_valid_partition(P, levels, verbose = T)
#' P <- list(c(0, 2), c(1, 3))
#' is_valid_partition(P, levels, verbose = T)
is_valid_partition <- function(P, levels, check_regularity = TRUE, verbose = FALSE) {
  n    <- length(levels)
  nlev <- lengths(levels)
  stride <- c(1, cumprod(nlev[-length(n)]))
  seqn <- seq_along(stride)
  lens    <- lengths(P)
  
  # for tracking which partitions that includes all levels of each variable
  includes_all_lev <- matrix(FALSE, nrow = length(P), ncol = n)

  for (j in seq_along(P)[lens > 0]) {
    p <- P[[j]]
    K <- p[[1]]
    for (pp in p) {
      
      # configuration of nodes in the current context 
      vals <- (pp%/%stride)%%nlev
      
      # the corresponding rows, setting each value to 0
      first_row <- pp-stride*vals
      for (i in seqn) {
        
        # list rows where the config of the co-parents of i 
        # are consistent with pp 
        # - a label on the edge from i to j would imply that all these rows 
        #   belongs to the same part p
        rows <- first_row[i] + stride[i]*levels[[i]]
        
        if (all(rows %in% p)) {
          K <- unique(c(K, rows))
          includes_all_lev[j, i] <- TRUE
        }
      }
    }
    
    if ( !(length(K) == length(p)) || !all(p == K) ) {
      if (verbose) cat("P is not CSI-consistent\n")
      return(FALSE)
    } 
  }
  
  if (check_regularity && all(lens >= min(nlev)) && any(colMeans(includes_all_lev) == 1)) {
    if (verbose) cat("P is not regular\n")
    return(FALSE)
  } else {
    return(TRUE)
  }
}
