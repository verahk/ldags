#' Compute local score in a labeled DAG from data
#' 
#' Compute sample counts over the relevant subset of variables and compute score.
#' If more than one parent return the score of the optimal local CSI-structure.
#' 
#' @param j node to be scored
#' @param parentsnodes parents of this node
#' @param n number of nodes
#' @param param object of `scoreparameters`
#'
#' @return score of node `j` given `parentnodes` and optimal CSI-structure 
score_from_lookup <- function(data, nlev, j, parentnodes, env) {
  
  if (length(parentnodes) == 0) {
    paid <- "NA"
  } else {
    parentnodes <- sort(parentnodes)
    paid <- paste0(parentnodes, collapse = "|")
  }
 
  scoreobj <- env$scores[[j]][[paid]]
  
  if (is.null(scoreobj)) {
    scoreobj  <- local_CSI_score(data, nlev, j, parentnodes, env$params$ess, env$params$kappa)
    env$scores[[j]][[paid]] <- scoreobj
  }
  
  return(scoreobj$score)
}

init_lookup_scoretable <- function(n, params) {
  rlang::new_environment(list(scores = vector("list", n),
                                     params = params))
}

## test ----
if (FALSE) {

  nlev <- rep(2, 5)
  params <- list(ess = 1, kappa = .5)
  data <- replicate(length(nlev), sample(0:1, size = 100, T))

  
  # add values in wrapper function
  test <- function(data, nlev, j, parentnodes, env) {
    score_from_lookup(data, nlev, j, parentnodes, env)
  }
  
  tab <- init_lookup_scoretable(length(nlev), params)
  test(data, nlev, 2, 1, env = tab) 
  tab$scores[[2]][["1"]]

}



