

example_bn <- function(name) {
  if (name == "collider") {
    dag <- rbind(Z = c(0, 0, 1),
                 X = c(0, 0, 1),
                 Y = c(0, 0, 0))
    colnames(dag) <- rownames(dag)
    
   
    levels <- rep(list(0:1), 3)
    #names(levels) <- rownames(dag)
    nlev <- lengths(levels)
    n <- length(nlev)
    
    cpts <- list() 
    cpts$Z <- array(c(.5, .5), 2, levels["Z"])
    cpts$X <- array(c(.5, .5), 2, levels[c("X")])
    
    # encode some "weak" context-specific independence
    cpts$Y <- array(c(rep(c(.4, .6), 2), .5, .5,  .3, .7), 
                    c(2, 2, 2), 
                    levels[c("Y", "Z", "X")])
    cpts
  } else if (name == "minimal_parent") {
    dag <- rbind(L  = c(0, 1, 1, 0, 1),
                 Z1 = c(0, 0, 0, 1, 0),
                 Z2 = c(0, 0, 0, 1, 0),
                 X  = c(0, 0, 0, 0, 1),
                 Y  = rep(0, 5))
    colnames(dag) <- rownames(dag)
    
  } else if (name == "LDAG10") {
    n <- 10
    nlev <- rep(2, n)
    levels <- lapply(nlev-1, seq.int, from = 0)
    dag <- matrix(0, n, n)
    colnames(dag) <- rownames(dag) <- paste0("X", seq_len(n))
    
    
    dag[1, 2:5] <- 1
    dag[2, c(3, 7)] <- 1
    dag[3, c(4, 7)] <- 1
    dag[4, c(5, 7, 8)] <- 1
    dag[5, c(9)] <- 1
    dag[6, c(7, 10)] <- 1
    dag[7, c(8, 10)] <- 1
    dag[8, c(5, 9, 10)] <- 1
    dag[9, c(10)] <- 1
    
    labels <- matrix(list(), n, n)
    labels[[2, 3]] <- rbind(0)
    labels[[1, 4]] <- rbind(1)
    labels[[4, 5]] <- rbind(c(0, NA))
    labels[[8, 5]] <- rbind(c(0, NA))
    labels[[2, 7]] <- rbind(c(1, 1, 0))
    labels[[3, 7]] <- rbind(c(0, 1, 1), c(1, NA, 1))
    labels[[4, 7]] <- rbind(c(1, 1, NA))
    labels[[6, 7]] <- rbind(c(1, 1, NA))
    labels[[5, 9]] <- rbind(1)
    labels[[7, 10]] <- rbind(c(1, NA, NA))
    labels[[8, 10]] <- rbind(c(1, NA, NA))
    labels[[9, 10]] <- rbind(c(1, NA, NA))
    
    partitions <- vector("list", ncol(dag))
    for (i in 1:ncol(dag)) {
      pa   <- which(dag[, i] == 1)
      if (length(pa) > 1) {
        tmp <- labels[pa, i]
        if (any(lengths(tmp) > 0)) {
          partitions[[i]] <- labels_to_partition(tmp, levels[pa], type = "outcome_vectors")
        }
      }
    }
    nlev <- rep(2, ncol(dag))
    bn <- rand_bn(dag, partitions, nlev, alpha = 1)
    return(bn)
    
    # test
    # i <- 7
    # pa  <- which(dag[, i] == 1)
    # tmp <- matrix(aperm(bn[[i]]$prob, c(1:length(pa)+1, 1)), ncol = nlev[i])
    # cbind(expand.grid(levels[pa]), get_parts(partitions[[i]]), tmp)
  }
  
  # return bn.fit object
  g <- bnlearn::empty.graph(colnames(dag))
  bnlearn::amat(g) <- dag
  bn <- bnlearn::custom.fit(g, cpts)
  return(bn)
}

