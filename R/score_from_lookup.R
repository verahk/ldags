#' Compute local score in a labeled DAG from data
#' 
#' Compute sample counts over the relevant subset of variables and compute score.
#' If more than one parent return the score of the optimal local CSI-structure.
#' 
#' @param data 
#' @param nlev 
#' @param parentnodes 
#' @param method 
#' @param table 
#' @param j node to be scored
#'
#' @return score of node `j` given `parentnodes` and optimal CSI-structure 
score_from_lookup <- function(data, levels, nlev = lengths(levels), j, parentnodes, ess, method, regular, lookup) {
  
  if (length(parentnodes) > 1) {
    
    # get score from look up table
    parentnodes <- sort(parentnodes)
    id <- paste0(c(j, parentnodes), collapse = "|")
    score <- lookup[[method]][[id]]
    
    # if not already computed, compute score and store in lookup-table
    if (is.null(score)) {
      
      score <- compute_local_bdeu_score(data, levels, nlev, j, parentnodes, ess, 
                                        method = method, regular = regular)
      lookup[[method]][[id]] <- score
    }
    
    return(score)
  } else {
    compute_local_bdeu_score(data, levels, nlev, j, parentnodes, ess)
  }
}




