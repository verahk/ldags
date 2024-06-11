

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
compute_local_bdeu_score <- function(counts, ess, r = ncol(counts), q = nrow(counts), s = 1) {
  sum(famscore_bdeu_byrow(counts, ess, r, q, s))
}

compute_local_bdeu_score_from_cpt <- function(cpt, ess = 1) {
  if (is.null(cpt$partition)) {
    compute_local_bdeu_score(cpt$counts, ess = ess, r = ncol(cpt$counts), q = nrow(cpt$counts), s = 1)
  } else {
    nlev <- lengths(cpt$levels)
    n <- length(nlev)
    s <- lengths(cpt$partition)
    compute_local_bdeu_score(cpt$counts, ess = ess, r = nlev[n], q = prod(nlev[-n]), s = s)
  }
}

compute_local_bdeu_score_from_data <- function(data, levels, nlev = lengths(levels), j, parentnodes, ess, struct = NULL, regular = FALSE, lookup = NULL) {
  counts <- compute_freq_table(data, nlev, j, parentnodes, lookup = NULL)
  if (length(parentnodes) > 1 && !is.null(struct)) {
    struct <- optimize_partition(counts, levels[parentnodes], ess, method = struct, regular = regular)
    score  <- sum(struct$scores)
    attr(score, "structure") <- struct
    return(score)
  } else {
    compute_local_bdeu_score(counts, ess, nlev[j], prod(nlev[parentnodes]), s = 1)
  }
}

