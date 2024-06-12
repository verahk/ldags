

#' Title
#'
#' @param counts 
#' @param ess 
#' @param r 
#' @param q 
#' @param s 
#'
#' @return
#' @export
#'
#' @examples
#' 
compute_local_bdeu_score <- function(data, levels, nlev = lengths(levels), j, parentnodes, ess, method = NULL, ...) {
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

compute_local_bdeu_score_from_cpt <- function(cpt, ess = 1) {
  if (is.null(cpt$partition)) {
    compute_local_bdeu_score_from_counts(cpt$counts, ess = ess)
  } else {
    nlev <- lengths(cpt$levels)
    n <- length(nlev)
    s <- lengths(cpt$partition)
    compute_local_bdeu_score_from_counts(cpt$counts, ess = ess, r = nlev[n], q = prod(nlev[-n]), s = s)
  }
}



