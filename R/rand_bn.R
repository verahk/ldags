

#' Title
#'
#' @param dag 
#' @param nlev 
#' @param alpha 
#' @param labels 
#'
#' @return
#' @export
#'
#' @examples
#' dag <- rbind(L  = c(0, 1, 1, 0, 1),
#'              Z1 = c(0, 0, 0, 1, 0),
#'              Z2 = c(0, 0, 0, 1, 0),
#'              X  = c(0, 0, 0, 0, 1),
#'              Y  = rep(0, 5))
#' colnames(dag) <- rownames(dag)
#' nlev <- rep(2, ncol(dag))
#' 
#' # label structure
#' labels <- vector("list", 5)
#' labels[[4]] <- list(0, 0)
#' labels[[5]] <- list(0, NULL)
#' 
#' rand_bn(dag, nlev, alpha = 1, labels = labels)
rand_bn <- function(dag, partitions = NULL, nlev = rep(2, ncol(dag)), alpha = 1) {
  vars <- colnames(dag)
  
  # sample cpts
  cpts <- vector("list", ncol(dag))
  names(cpts) <- vars
  for (j in seq_along(vars)) {
    pa <- which(dag[, j] == 1)
    r <- nlev[j]
    
    if (is.null(partitions[[j]])) {
      q <- prod(nlev[pa])
      p <- replicate(q, bida:::rDirichlet(1, rep(alpha, r), r), simplify = T)
    } else {
      q <- length((partitions[[j]]))
      indx <- get_parts(partitions[[j]])
      p <- replicate(q, bida:::rDirichlet(1, rep(alpha, r), r), simplify = T)[, indx]
    }
    dim(p) <- nlev[c(j, pa)]
    dimnames(p) <- setNames(vector("list", 1+length(pa)), vars[c(j, pa)])
    
    cpts[[j]] <- p
  }

  # return bn.fit object
  g <- bnlearn::empty.graph(vars)
  bnlearn::amat(g) <- dag
  bnlearn::custom.fit(g, cpts)
}
