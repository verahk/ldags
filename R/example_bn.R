

example_bn <- function() {
  dag <- rbind(Z = c(0, 0, 1),
               X = c(0, 0, 1),
               Y = c(0, 0, 0))
  colnames(dag) <- rownames(dag)
  
  cpts <- setNames(vector("list", ncol(dag)), colnames(dag))
  levels <- rep(list(0:1), 3)
  names(levels) <- names(cpts)
  nlev <- lengths(levels)
  n <- length(nlev)
  
  
  cpts$Z <- array(c(.5, .5), 2, levels["Z"])
  cpts$X <- array(c(.5, .5), 2, levels[c("X")])
  
  # encode some "weak" context-specific independence
  cpts$Y <- array(c(rep(c(.4, .6), 2), .5, .5,  .3, .7), 
                  c(2, 2, 2), 
                  levels[c("Y", "Z", "X")])
  cpts
  g <- bnlearn::empty.graph(names(cpts))
  bnlearn::amat(g) <- dag
  
  # return bn.fit object
  return(bnlearn::custom.fit(g, cpts))
}

