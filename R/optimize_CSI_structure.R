
#' Optimize local CSI structure
#'
#' @param data (integer matrix)
#' @param nlev (integer vector)
#' @param j (integer)
#' @param parentnodes (integer vector)
#' @param ess (numeric constant)
#'  imaginary sample size for the parameter prior.
#' @param kappa (numeric constant)
#'  penalization factor for CSI-complexity. Search is terminated when the maximal 
#'  improvement in the score from fusing two parts is less than `-(r-1)*log(kappa)`, 
#'  where `r-1` is the reduction in the number of free parameters in the reduced CPT.
#' @return a list with the optimal partition and the associated scores
#' @export
#'
#' @examples
optimize_CSI_structure <- function(data, nlev, j, parentnodes, ess = 1, kappa = .1, verbose = FALSE) {
  
  r <- nlev[j]
  q <- prod(nlev[parentnodes])
  
  # compute marginal counts
  stride <- c(1, cumprod(nlev[parentnodes]))
  counts <- tabulate(data[, c(parentnodes, j)]%*%stride+1, q*r)
  dim(counts) <- c(q, r)
  
  # optimize CSI-structure
  P <- as.list(seq.int(0, q-1))
  levels <- lapply(nlev[parentnodes]-1, seq.int, from = 0)
  optimize_CSI_structure_greedy(counts, levels, nlev[parentnodes], P, ess = ess, lkappa = log(kappa), verbose = verbose)
}


#' @rdname optimize_CSI_structure
#' @param counts (integer matrix)
#'  observed frequency table.
#' @param levels (list of integer vectors)
#'  levels of variables.
#' @param P (list of integer vectors)
#'  current partition.
#' @param lkappa (numeric constant)
#' @param verbose (logic)
#'  if TRUE, print intermediate results from search
#' @examples
#' levels <- list(0:1, 0:1)
#' counts <- cbind(c(10, 20, 30, 40), rep(10, 4))
#' P <- as.list(seq_len(nrow(counts))-1)
#' optimize_CSI_structure_greedy(counts, levels, P, ess = 1, lkappa = 0, verbose = T)
optimize_CSI_structure_greedy <- function(counts, levels, nlev = lengths(levels), P, ess, lkappa, verbose = FALSE){
  
  q   <- nrow(counts)
  r   <- ncol(counts)
  
  # compute seach part's contribution to the local score
  part_size   <- lengths(P)
  part_counts <- rowsum(counts, unlist_partition(P, part_size), reorder = T)
  part_scores <- famscore_bdeu_byrow(part_counts, ess, r, q, part_size)
  
  # find the two current parts that increase the score most if collapsed
  keep_climb <- FALSE
  if (length(P) > 2) {
    max_diff <- -(r-1)*lkappa
    for (i in seq.int(1, length(P)-1)) {
      for (j in seq.int(i+1, length(P))) {
        
        # candidate partition
        PP  <- replace(P, i, list(c(P[[i]], P[[j]])))
        PP[j] <- NULL
   
        # compute score of collapsed parts
        tmp_count <- part_counts[i,, drop = FALSE] + part_counts[j, , drop = FALSE]
        tmp_score <- famscore_bdeu_byrow(tmp_count, ess, r, q, length(PP[[i]]))
        
        # compare score with score of non-collapsed regions
        diff <- tmp_score - (part_scores[i] + part_scores[j])
        
        # compare score with current best partition
        if (diff > max_diff) {
          
          # check if collapsing part i and j result in regular, CSI-consistent partition
          if (!is_valid_partition(PP, levels, nlev, verbose = F))  next
          
          if (verbose) cat("\n (i, j):", c(i, j), "diff:", diff, "partition:", unlist_partition(PP))
          max_diff <- diff
          best_partition <- PP
          keep_climb <- TRUE
        }
      }
    }
  }
    
  if (keep_climb) {
    if (verbose) cat("\nbest partition:", unlist_partition(best_partition))
    optimize_CSI_structure_greedy(counts, levels, nlev, best_partition, ess, lkappa, verbose)
  } else {
    return(list(partition = P, counts = part_counts, scores = part_scores))
  }
}

# alternative: re-use counts and scores
optimize_CSI_structure_greedy_2 <- function(P, counts, scores, levels, nlev = lengths(levels), ess, lkappa, verbose = FALSE){
  
  r   <- ncol(counts)
  q   <- prod(nlev)
  
  # find the two current parts that increase the score most if collapsed
  keep_climb <- FALSE
  if (length(P) > 2) {
    max_diff <- -(r-1)*lkappa
    for (i in seq.int(1, length(P)-1)) {
      for (j in seq.int(i+1, length(P))) {
        
        # candidate partition
        PP  <- replace(P, i, list(c(P[[i]], P[[j]])))
        PP[j] <- NULL
        
        
        # compute score of collapsed parts
        tmp_count <- counts[i,, drop = FALSE] + counts[j, , drop = FALSE]
        tmp_score <- famscore_bdeu_byrow(tmp_count, ess, r, q, length(PP[[i]]))
        
        # compare score with score of non-collapsed regions
        diff <- tmp_score - (scores[i] + scores[j])
        
        # compare score with current best partition
        if (diff > max_diff) {

          
          # check if collapsing part i and j result in regular, CSI-consistent partition
          if (!is_valid_partition(PP, levels, nlev, verbose = F))  next
          if (verbose) cat("\n (i, j):", c(i, j), "diff:", diff, "partition:", unlist_partition(PP)) 
   
          max_diff <- diff
          istar <- i
          jstar <- j
          scorestar <- tmp_score
          keep_climb <- TRUE
        }
      }
    }
  }
  
  if (keep_climb) {
    i <- istar
    j <- jstar
    counts[i, ] <- counts[i, ] + counts[j, ]
    counts <- counts[-j, ]
    scores[i] <- scorestar
    scores <- scores[-j]
    P[[i]] <- unlist(P[c(i, j)])
    P[[j]] <- NULL
    if (verbose) cat("\nbest partition:", unlist_partition(P))
    optimize_CSI_structure_greedy_2(P, counts, scores, levels, nlev = lengths(levels), ess, lkappa, verbose)
  } else {
    return(list(partition = P, counts = counts, scores = scores))
  }
}


### test ----
if (FALSE) {
  levels <- rep(list(0:2), 3)
  nlev <- lengths(levels)
  q <- prod(lengths(levels))
  P <- as.list(seq_len(q)-1)
  r <- 2
  ess <- 1
  kappa <- 1
  
  ##  same proportion, different counts
  counts <- matrix(seq(1, q, length.out = q), nrow = q, ncol = 2)
  counts 
  scores <- famscore_bdeu_byrow(counts, ess, r, q, 1)
  
  optimize_CSI_structure_greedy(counts, levels, nlev, P, 1, 0, verbose = T)
  optimize_CSI_structure_greedy_2(P, counts, scores, levels, nlev, 1, 0, verbose = T)
  
  microbenchmark::microbenchmark(optimize_CSI_structure_greedy(counts, levels, nlev, P, 1, 0, verbose = F), 
                                 optimize_CSI_structure_greedy_2(P, counts, scores, levels, nlev, 1, 0, verbose = F), 
                                 times = 20)
  
  profvis::profvis(
    optimize_CSI_structure_greedy(counts, levels, nlev, P, 1, 0, verbose = F)
  )
  profvis::profvis(
    optimize_CSI_structure_greedy_2(P, counts, scores, levels, nlev, 1, 0, verbose = F)
  )
  
  res <- optimize_CSI_structure_greedy(counts, levels, nlev, P, 1, 0, verbose = T)
  res <- optimize_CSI_structure_greedy_2(P, counts, scores, levels, nlev, 1, 0, verbose = T)
  # regions with high support are combined first. 
  

  ## different proportion, same counts
  prop <- c(.5, .4, .3, .2, .1, .05, .025, .01)
  prop <- cbind(prop, 1-prop)
  res <- optimize_CSI_structure_greedy(prop*10, levels, nlev, P, 1, 0, verbose = T)
  res <- optimize_CSI_structure_greedy(prop*100, levels, nlev, P, 1, 0, verbose = T)
  res <- optimize_CSI_structure_greedy(prop*1000, levels, nlev, P, 1, 0, verbose = T)
  # larger sample size, less willing to fuse parts
  
}

if (FALSE) {
  
  j <- 10 
  parentnodes <- 1:6 #match(bn[[j]]$parents, names(bn))
  profvis::profvis(optimize_CSI_structure(data, nlev, 10, parentnodes, verbose = T))
}
