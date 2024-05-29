
new_CPT <- function(node, 
                    parentnodes, 
                    levels, 
                    partition = NULL,
                    counts = NULL, 
                    score = NULL) {
  x <- list(node = node, 
            parentnodes = parentnodes, 
            levels = levels, 
            partition = partition,
            counts = counts,
            score = score)
  return(x)
}

cpt_from_dag <- function(dag, levels, labels = NULL, node) {
  parentnodes <- which(dag[, node] == 1)
  if (!is.null(labels[[node]])) {
    partition <- labels_to_partition(labels[[node]], levels[parentnodes])
  } else {
    partition <- NULL
  }
  new_CPT(j, parentnodes, levels[c(parentnodes, node)], partition)
}




