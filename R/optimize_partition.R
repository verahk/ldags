#' Optimize partition over parent outcome space
#' 
#' @param counts (integer matrix)
#'  a frequency table.
#' @param levels (list of integer vectors)
#'  levels of each conditioning variable in `counts`, 
#'  such that `expand.grid(levels)` gives the joint configurations corresponding
#'  to each row in `counts`.
#' @param ess (numeric)
#'  imaginary sample size for the bdeu-prior.
#' @param method (character)
#' @param make_regular (logical)
#'  if `FALSE` (default) the optimized partition is returned, also if it implies
#'  conditional independencies (is not "regular"). If `TRUE` each part of the 
#'  optimized partition is divided by the outcomes of the relevant variables.
#' @param verbose (logical)
#' @return a list with named elements: 
#' - `partition`: the partition implied by the tree
#' - `scores`: each part's contribution to the local score
#' -  additional output from the different optimization procedures.n
#' @export
optimize_partition <- function(counts, levels, ess, method, make_regular = F, verbose = FALSE){
  method <- match.arg(method, c("tree", "ldag", "part"))
  res <- switch(method, 
               "tree" = optimize_partition_tree(counts, levels, ess, verbose = verbose),
               "ldag" = optimize_partition_ldag(counts, levels, ess, verbose = verbose),
               "part" = optimize_partition_part(counts, levels, ess, verbose = verbose))
  
  if (make_regular) {
    # ensure that partition is regular
    tmp <- make_partition_regular(res$partition, lengths(levels))
    if (!length(tmp) == length(res$partition)) {
      res$partition <- tmp
      parts <- unlist_partition(tmp)
      scores <- famscore_bdeu_byrow(rowsum(counts, parts), ess, r = ncol(counts), q = nrow(counts), s = lengths(tmp))
      res$scores <- scores
    }
  }
  return(res)
}

#' @rdname optimize_partition
optimize_partition_from_cpt <- function(cpt, ess, method, make_regular = T, return_cpt = TRUE) {
  n <- length(cpt$levels)
  if (n > 2) {
    res <- optimize_partition(cpt$counts, cpt$levels[-n], ess, method, make_regular, verbose = FALSE)
    if (return_cpt) {
      cpt$partition <- res$partition
      cpt$counts <- res$counts
      cpt$score  <- sum(res$scores)
      return(cpt)
    } else {
      res$partition
    }
  } else if (return_cpt) {
    return(cpt)
  }
}