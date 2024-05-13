

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
    
    
  } else if (FALSE) {
    dag <- rbind(Z1 = c(0, 0, 1, 0),
                 Z2 = c(0, 0, 1, 1),
                 X  = c(0, 0, 0, 1),
                 Y  = rep(0, 4))
    colnames(dag) <- rownames(dag)
    return(dag)
  }
  
  # return bn.fit object
  g <- bnlearn::empty.graph(colnames(dag))
  bnlearn::amat(g) <- dag
  bn <- bnlearn::custom.fit(g, cpts)
  return(bn)
}

