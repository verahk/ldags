

#' Sample DAGs with local structure
#' 
#' Apply hybrid sampling scheme provided by the `BiDAG`-package with user-defined
#' scoring functions to sample DAGs with optimized local structure.
#' Use `bnlearn::hc` to definine a CPDAG for the initial search space, limiting 
#' the maximum number of parents.
#'
#' @param data (integer matrix)
#' @param nlev (integer vector)
#' @param algo (character)
#'  which MCMC-scheme to employ for sampling DAGs given a search space. Either
#'  `"partition"` or `"order"`. See [BiDAG::sampleBN].
#' @param ess (numeric)
#' @param edgepf (numeric)
#'  edge penalty
#' @param hardlimit (integer)
#'  maximum number of parents.
#' @param local_struct (character)
#'  greedy optimization routine to optimize the local structure. Either `tree` or 
#'  `ldag`. If `NULL` (default), no local structure is optimized.
#' @param verbose (logical)
#' @param lookup (environment)
#'  a environment in which the optimized local structures is stored. 
#'  Generated e.g. with [rlang::new_environment]. See [score_from_lookup].
#' 
#' @return 
#' @export
#'
#' @examples
#' bn <- example_bn("LDAG10")
#' data <- bida:::sample_data_from_bn(bn, 1000)
#' nlev <- rep(2, 10) # binary data 
#' lookup <- rlang::new_environment() 
#' 
#' smpl_tree <- sample_dags(data, nlev, algo = "order", local_struct = "tree", verbose = T, lookup = lookup)
#' smpl_ldag <- sample_dags(data, nlev, algo = "order", local_struct = "ldag", verbose = T, lookup = lookup)
#' smpl_dag  <- sample_dags(data, nlev, algo = "order", local_struct = NULL, verbose = T, lookup = lookup)    
#' 
#' edgep <- cbind(tree = c(BiDAG::edgep(smpl_tree)), 
#'                ldag = c(BiDAG::edgep(smpl_ldag)),
#'                dag  = c(BiDAG::edgep(smpl_dag)))
#'                
#' pairs(edgep)
sample_dags <- function(data, nlev, algo = "partition", ess = 1, edgepf = 1, hardlimit = 5, local_struct = NULL, verbose = FALSE, lookup = NULL) {
  hardlimit <- min(ncol(data)-1, hardlimit)
  
  # genereate data.frame, removing unobserved levels in data, as required by BiDAG
  df <- data.frame(apply(data, 2, factor, exclude = NULL, simplify = FALSE))
  
  # define scoreparameters 
  if (is.null(local_struct)) {
    scorepar <- BiDAG::scoreparameters("bdecat", 
                                       data = df, 
                                       bdecatpar = list(chi = ess, 
                                                        edgepf = edgepf))
    
  } else {
    
    scorepar <- BiDAG::scoreparameters("usr", 
                                       data = data, 
                                       usrpar = list(pctesttype = "bdecat"))
    
    # additional parameters for user-specified score params 
    scorepar$ess  <- ess 
    scorepar$edgepf <- edgepf
    scorepar$Cvec <- nlev
    scorepar$levels <- lapply(nlev-1, seq.int, from = 0)
    scorepar$local_struct <- local_struct
    
    # define environment in which to store optimized local structures
    if (!is.null(lookup)) {
      tmp <- list()
      assign(local_struct, tmp, lookup)
      scorepar$lookup <- lookup
    }
    
    
    # define scoring function
    usrDAGcorescore <- function(j, parentnodes, n, scorepar) {
      score_from_lookup(scorepar$data, 
                        levels = scorepar$levels, 
                        nlev = scorepar$Cvec, 
                        j, 
                        parentnodes, 
                        ess = scorepar$ess,
                        method = scorepar$local_struct, 
                        lookup = scorepar$lookup) - length(parentnodes)*log(scorepar$edgepf)
    }
    assignInNamespace("usrDAGcorescore", usrDAGcorescore, ns = "BiDAG")
  }
  
  tic <- Sys.time()
  # learn skeleton
  # if (verbose) cat("\nLearn skeleton using pcalg::skeleton")
  # skel <- pcalg::skeleton(suffStat = list(dm = data, adaptDF = FALSE),
  #                         indepTest = pcalg::disCItest, alpha = 0.05,
  #                         labels = colnames(data), verbose = verbose, m.max = hardlimit) 
  # tic <- c(tic, skel = Sys.time())
  if (verbose) cat("\nLearn initial CPDAG using bnlearn::hc\n")
  cpdag <- bnlearn:::hc(df, maxp = hardlimit)
  skel <- bnlearn::amat(cpdag)
  tic <- c(tic, hc = Sys.time())
  
  # define search space 
  if (verbose) cat("\nDefine search space using BiDAG::learnBN\n")
  iterfit <- BiDAG::learnBN(scorepar, "orderIter", hardlimit = hardlimit, scoreout = T, verbose = verbose, startspace = skel)
  tic <- c(tic, learnBN = Sys.time())

  # sample DAGs 
  if (verbose) cat("\nDefine search space using BiDAG::sampleBN\n")
  smpl     <- BiDAG::sampleBN(scorepar, algo, scoretable = BiDAG::getSpace(iterfit), verbose = T)
  tic <- c(tic, sampleBN = Sys.time())
  
  # add time-tracking as an attribute
  attr(smpl, "toc") <- diff(tic)
  attr(smpl, "call") <- match.call()
  
  return(smpl)
}
