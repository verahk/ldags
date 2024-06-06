

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
#' stride <- c(1, 2)
#' 
#' # CSI-consistent, regular partition
#' P <- list(c(0, 1, 3), 2)
#' is_ldag_consistent(P,  nlev, stride)
#' is_regular(P, nlev, stride)
#' 
#' # not CSI-consistent partition
#' P <- list(c(0, 3), 2, 1)
#' is_ldag_consistent(P,  nlev, stride)
#' is_regular(P, nlev, stride)
#' 
#' # CSI-consistent, but not regular
#' P <- list(c(0, 1), c(2, 3))
#' is_ldag_consistent(P,  nlev, stride)
#' is_regular(P, nlev, stride)
#' 
#' P <- list(c(0, 2), c(1, 3))
#' is_ldag_consistent(P,  nlev, stride)
#' is_regular(P, nlev, stride)
#' 
#' # Different cardinality


is_ldag_consistent <- function(P, nlev, stride = c(1, cumprod(nlev[-length(nlev)]))) {
  nparts <- length(P)
  joint <- seq_len(prod(nlev))-1
  
  is_consistent <- TRUE
  j <- 0
  while(is_consistent && j < nparts) {
    j <- j+1
    is_consistent <- is_ldag_consistent_part(P[[j]], nlev, stride, joint)
  }
  return(is_consistent)
}

is_ldag_consistent_part <- function(p, nlev, stride, joint = seq_len(prod(nlev))-1) {
  if (length(p) == 1) return(TRUE)
  p <- p+1                        # refer to rows, starting at 1
  K <- logical(length(joint))
  for (i in seq_along(nlev)) {
    
    # possible labels on the edge i -> child
    labels <- joint[(joint%/%stride[i])%%nlev[i] == 0]
    
    # rows by each label (by row)
    rows  <- outer(labels, stride[i]*(seq_len(nlev[i])-1)+1, "+")
    
    # check which labels could produce p
    indx <- rowSums(array(match(rows, p, 0L), dim(rows)) > 0) == nlev[i]
    K[c(rows[indx, ])] <- TRUE
  }
  
  #print(list(p-1, K))
  all(K[p])
}


is_regular <- function(P, nlev, stride = c(1, cumprod(nlev[-length(nlev)]))) {
  is_regular <- !logical(length(nlev))
  if (min(lengths(P)) >= min(nlev)) {
    partition <- unlist_partition(P)
    joint     <- seq_len(prod(nlev))-1
    for (i in seq_along(nlev)) {
      contexts <- joint-stride[i]*(joint%/%stride[i])%%nlev[i]
      parts    <- split(partition, contexts)
      is_regular[i] <- any(! vapply(parts, function(x) all(x == x[1]), logical(1)))
    }
  }
  return(is_regular)
}

