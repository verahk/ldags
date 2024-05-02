


rand_cpt_array <- function(dag, labels, nlev, i, alpha, name_dims = TRUE) {
  
  if (is.character(i)) i <- match(i, colnames(dag))
  
  pa <- which(dag[, i] == 1)
  r  <- nlev[i]
 
  if (length(alpha) == 1) {
     alpha <- rep(alpha, r)
  }
  
  if (length(pa) == 0) {
    p <- bida:::rDirichlet(1, alpha, r)
  } else {
    if (is.null(labels[[i]])) {
      q  <- prod(nlev[pa])
      p  <- vapply(seq_len(q), function(x) bida:::rDirichlet(1, alpha, r), numeric(r))
      
    } else {
      stopifnot(length(labels[[i]]) == length(pa))
      P <- labels_to_partition(labels[[i]], nlev[pa])
      p <- vapply(seq_along(P), 
                  function(x) bida:::rDirichlet(1, alpha, r), numeric(r))[, unlist_partition(P)]
    } 
  }
 
  dim(p) <- nlev[c(i, pa)]
  if (name_dims) {
    scope <- colnames(dag)[c(i, pa)]
    stopifnot(!is.null(scope))
    dimnames(p) <- setNames(lapply(dim(p)-1, seq.int, from = 0), scope)
  }
  
  return(p)
}
