


famscore_bb_from_data <- function(j, parentnodes, data, levels, ess, kappa) {
  
  nlev   <- vapply(levels, max, integer(1)) +1
  r <- nlev[j]
  q <- prod(nlev[parentnodes])
  
  stride <- c(1, cumprod(nlev[parentnodes]))
  counts <- tabulate(data[, c(parentnodes, j)]%*%stride + 1, q*r)
  dim(counts) <- c(q, r)
  
  # run branch and bound
  P <- as.list(seq_len(q))
  B <- -Inf
  res <- opt_part_bb(counts, levels, P, 1, B, ess = ess, lkappa = log(kappa))
}

opt_part_bb <- function(counts, levels, P, f, B, ess = 1, lkappa = log(.5)) {
  
  r <- ncol(counts)
  q <- nrow(counts)
  
  # compute score
  part_sizes  <- lengths(P)
  part_counts <- rowsum(counts, rep.int(seq_along(P), part_sizes))
  loglik <- sum(famscore_bdeu_byrow(part_counts, ess, r, q, part_sizes))
  score  <- loglik + (r-1)*(q-length(P))*lkappa
  
  if (score > B) {
    if (is_CSI_consistent(P, levels)) {
      B <- score
    }
  } else {
    # collapsing more rows will neither increase log-lik nor reduce penalty
    return(list(B = B, P = P, f= f))
  }
  
  # maximum score of further collapsing rows = max loglik + minimum penalty
  # maxScore <- loglik + (r-1)*(q-length(P))*lkappa   
  # if (maxScore < B) {
  #   return(list(B = B, P = P, f= f))
  # }
  
  
  for (i in seq_len(f)) {
    PP  <- replace(P, i, list(c(P[[i]], P[[f+1]])))
    PP[f+1] <- NULL
    res <- opt_part_bb(counts, levels, PP, f, B, ess, lkappa)
    B <- res$B
  }
  
  opt_part_bb(counts, levels, P, f+1, B, ess, lkappa)
}


