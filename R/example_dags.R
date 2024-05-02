


example_dags <- function(x) {
  if (x == "minimal_parent") {
     dag <- rbind(L  = c(0, 1, 1, 0, 1),
                  Z1 = c(0, 0, 0, 1, 0),
                  Z2 = c(0, 0, 0, 1, 0),
                  X  = c(0, 0, 0, 0, 1),
                  Y  = rep(0, 5))
     colnames(dag) <- rownames(dag)
     return(dag)
  } else if (x == "minimal_parent_small") {
    dag <- rbind(Z1 = c(0, 0, 1, 0),
                 Z2 = c(0, 0, 1, 1),
                 X  = c(0, 0, 0, 1),
                 Y  = rep(0, 4))
    colnames(dag) <- rownames(dag)
    return(dag)
  }
}

