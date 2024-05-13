

new_CPT <- function(j, parentnodes, counts = NULL, score = NULL, partition = NULL, params = NULL) {
  list(node = j, 
       parentnodes = parentnodes, 
       score = score,
       counts = counts, 
       partition = partition,
       params = params)
}

#' Title
#'
#' @param data 
#' @param nlev 
#' @param j 
#' @param parentnodes 
#' @param params 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' nlev <- rep(2, 5)
#' lev  <- lapply(nlev-1, seq.int, from = 0)
#' data <- sapply(lev, sample, size = 10, replace = T)
#' j <- 1
#' parentnodes <- 2:4
CPT_from_data <- function(data, nlev, j, parentnodes, params, verbose = FALSE) {
  
  npar <- length(parentnodes) 
  if (npar == 0) {
    counts <- matrix(tabulate(data[, node] +1, nlev[j]), ncol = nlev[j])
    new_CPT(j, parentnodes, counts = counts, params = params)
  } else if (npar == 1) {
    r <- nlev[j]
    q <- nlev[parentnodes]
    counts <- tabulate(data[, c(parentnodes, j)]%*%c(1, q) +1, q*r)
    dim(counts) <- c(q, r)
    new_CPT(j, parentnodes, counts = counts, params = params)
  } else {
    r <- nlev[j]
    q <- prod(nlev[parentnodes])
    
    # compute marginal counts
    stride <- c(1, cumprod(nlev[parentnodes]))
    counts <- tabulate(data[, c(parentnodes, j)]%*%stride+1, q*r)
    dim(counts) <- c(q, r)
    
    if (optimize_CSI_structure) {
      
      # optimize CSI-structure
      P <- as.list(seq.int(0, q-1))
      levels <- lapply(nlev[parentnodes]-1, seq.int, from = 0)
      tmp <- optimize_CSI_structure_greedy(counts, levels, nlev, P, ess = params$ess, lkappa = log(params$kappa), verbose = verbose)
      
      # return params over optimized structure
      new_CPT(j, parentnodes, counts = tmp$counts, score = tmp$score, partition = tmp$partition, params = params)
    } else {
      new_CPT(j, parentnodes, counts = counts, params = params)
    }
  }
}


score_CPT <- function(x) {
  if (is.null(x$score)) {
    sum(famscore_bdeu_byrow(x$counts, x$params$ess))
  } else {
    x$score
  }
}

mean_CPT <- function(x) {
  if (is.null(x$partition)) {
    alpha <- x$ess/length(x$counts)
  } else {
    alpha <- x$ess/length(unlist(x$partition))*lengths(x$partition)
  }
  upd_counts <- alpha + x$counts
  
  return((upd_counts)/rowSums(upd_counts))
}

