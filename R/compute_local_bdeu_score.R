

#' Compute local bdeu score 
#'
#' @param data 
#' @param levels 
#' @param nlev 
#' @param j (integer)
#'  column index of node
#' @param parentnodes (integer vector)
#'  column indicies of parent nodes. A vector of length 0 (e.g. `integer(0)`) 
#'  translates to no parents.
#' @param ess (numeric constant)
#'  equivalent sample size.
#' @param method (character)
#'  algorithm to learn local structure. See \link{optimize_partition}. 
#'  If `NULL` (default) no local structure is learned. Ignored if `length(parentnodes) < 2`.
#' @param ... additional parameters sent to [optimize_partition()].
#'
#' @return a numeric constant, the Bdeu-score. If a local structure is inferred,
#'  this is returned as an named attribute (`structure`).
#' @export
#'
#' @examples


compute_local_bdeu_score <- function(data, levels, nlev = lengths(levels), j, parentnodes, ess = 1, method = NULL, ...) {
  counts <- compute_freq_table(data, nlev, j, parentnodes)
  if (length(parentnodes) > 1 && !is.null(method)) {
    struct <- optimize_partition(counts, levels[parentnodes], ess, method = method, ...)
    score  <- sum(struct$scores)
    attr(score, "structure") <- struct
    return(score)
  } else {
    sum(famscore_bdeu_byrow(counts, ess))
  }
}

compute_local_bdeu_score_from_cpt <- function(cpt, ess) {
  if (is.null(cpt$partition)) {
    sum(famscore_bdeu_byrow(counts, ess))
  } else {
    nlev <- lengths(cpt$levels)
    n <- length(nlev)
    s <- lengths(cpt$partition)
    sum(famscore_bdeu_byrow(cpt$counts, ess = ess, r = nlev[n], q = prod(nlev[-n]), s = s))
  }
}



