

#' Title
#'
#' @param scorepar (object of class [BiDAG::scoreparameters])
#'  See \link{define_scorepar}.
#' @param algo_init (character)
#'  Name of optimization routine for local structure. 
#'  See \link{optimize_partition}.
#' @param algo_sample (character)
#'  Name of MCMC implementation. Either `partition` or `order`. 
#'  See [BiDAG::sampleBN].
#' @param hardlimit (integer)
#'  Maximal number of parents allowed. See [BiDAG::learnBN].
#' @param verbose (logical)
#' @return
#' @export
#'
#' @examples
sample_dags <- function(scorepar, algo_init = "pcskel", algo_sample = "order", hardlimit = 5, verbose = F) {
  
  tic <- Sys.time()
  
  # init search space 
  if (verbose) cat("\nInit search space:")
  startspace <- init_search_space(scorepar, algo_init, hardlimit, verbose = verbose)
  tic <- c(tic, init = Sys.time())
  
  # expand search space
  if (verbose) cat("\nExpand search space using BiDAG::iterativeMCMC:")
  #iterfit <- BiDAG::iterativeMCMC(scorepar, hardlimit = hardlimit, startspace = startspace, scoreout = T, verbose = verbose)
  iterfit <- BiDAG::learnBN(scorepar, "orderIter", hardlimit = 5, startspace = startspace, scoreout = T, verbose = verbose)
  tic <- c(tic, expand = Sys.time())
  
  # run MCMC 
  if (verbose) cat("\nSample DAGs using BiDAG::sampleBN:\n")
  smpl <- BiDAG::sampleBN(scorepar, algo = algo_sample, scoretable = BiDAG::getSpace(iterfit), verbose = verbose)
  tic  <- c(tic, sample = Sys.time())
  
  # add time-tracking as an attribute
  attr(smpl, "toc") <- diff(tic)
  
  return(smpl)
  
}
