

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
rand_bn <- function(dag, labels = NULL, nlev = rep(2, ncol(dag)), alpha = 1) {
  vars <- colnames(dag)
  cpts <- lapply(vars, function(v) rand_cpt_array(dag, labels, nlev, v, alpha, TRUE))
  names(cpts) <- vars
  
  # return bn.fit object
  g <- bnlearn::empty.graph(vars)
  bnlearn::amat(g) <- dag
  bnlearn::custom.fit(g, cpts)
}
