test_that("multiplication works", {
  

  dag  <- matrix(0, 3, 3)
  dag[upper.tri(dag)] <- 1
  colnames(dag) <- rownames(dag) <- c("z", "x", "y")
  
  dage <- dag
  dage["x", "y"] 
  partitions  <- list(y = list(c(0, 2), 1, 3))
  bn <- rand_bn(dag, partitions = partitions)  
  
  compute_scores <- function(data, dag, method, regular) {
    lapply(seqn, 
           function(j) compute_local_bdeu_score(data, levels, nlev, j, seqn[dag[, j] == 1],
                                                method = method, regular = regular))
  }
   
 data <- bida:::sample_data_from_bn(bn, 10**4)
 tree <- compute_scores(data, dag, "tree", TRUE) 
    
})
