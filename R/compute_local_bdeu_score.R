

#' Title
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
#' 
#' levels <- rep(list(0:1), 3)
#' nlev   <- lengths(levels)
#' data   <- sapply(levels, sample, size = 10, replace = T)
#' 
#' # no parent
#' compute_local_bdeu_score(data, levels, nlev, 1, integer(0))
#' 
#' # single parent
#' compute_local_bdeu_score(data, levels, nlev, 1, 2)
#' 
#' # multiple parents
#' compute_local_bdeu_score(data, levels, nlev, 1, 2:3, method = "ldag")

compute_local_bdeu_score <- function(data, levels, nlev = lengths(levels), j, parentnodes, ess = 1, method = NULL, ...) {
  counts <- compute_freq_table(data, nlev, j, parentnodes)
  if (length(parentnodes) > 1 && !is.null(method)) {
    struct <- optimize_partition(counts, levels[parentnodes], ess, method = method, ...)
    score  <- sum(struct$scores)
    attr(score, "structure") <- struct
    return(score)
  } else {
    compute_local_bdeu_score_from_counts(counts, ess)
  }
}

compute_local_bdeu_score_from_counts <- function(counts, ess, r = ncol(counts), q = nrow(counts), s = 1) {
  sum(famscore_bdeu_byrow(counts, ess, r, q, s))
}

compute_local_bdeu_score_from_cpt <- function(cpt, ess) {
  if (is.null(cpt$partition)) {
    compute_local_bdeu_score_from_counts(cpt$counts, ess = ess)
  } else {
    nlev <- lengths(cpt$levels)
    n <- length(nlev)
    s <- lengths(cpt$partition)
    compute_local_bdeu_score_from_counts(cpt$counts, ess = ess, r = nlev[n], q = prod(nlev[-n]), s = s)
  }
}



