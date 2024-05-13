


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
  } else if (x == "LDAG_10_nodes") {
    n <- 10
    nlev <- rep(2, n)
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
    
    
    # g <- as(dag, "graphNEL")
    # 
    # indx <- sapply(labels[dag == 1], is.null)
    # prettify_lab <- function(x) {
    #   tmp <- apply(x, 1, function(y) sprintf("(%s)",  paste(y, collapse = ",")))
    #   sprintf("{%s}", paste(gsub("NA", "*", tmp), collapse = ","))
    # }
    # 
    # edgeAttr <- sapply(labels[dag == 1][!indx], prettify_lab)
    # names(edgeAttr) <- Rgraphviz::edgeNames(g)[!indx]
    # Rgraphviz::plot(g, edgeAttrs = list(label = edgeAttr))
    
    # store labels as nested list
    labels <- lapply(seq_len(n), function(i) labels[which(dag[, i] == 1), i])
    return(list(labels, DAG))
  }
}

