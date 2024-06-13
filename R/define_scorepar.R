#' Title
#'
#' @param data (matrix)
#'  each column is assumed to be outcomes of a categorical variable, with 
#'  possible outcomes `0, ..., K-1.`
#' @param nlev (integer vector)
#'  cardinality of each variable.
#' @param ess (numeric constant)
#'  equivalent sample size
#' @param edgepf (numeric constant)
#'  a factor `-log(edgepf)*length(parentnodes)` is added to the marginal likelihood score.
#' @param local_struct, regular 
#'  parameters controlling the local-structure-optimization routine. See [optimize_partition()].
#' @param lookup (environment)
#'  an environment in which the optimized structures are saved, initiated by rlang:::new_environment().
#'  If `NULL` the structures are not  saved.
#'
#' @return an object of class [BiDAG::scoreparameters]
#' @export
#'
#' @examples
define_scorepar <- function(data, nlev, ess = 1, edgepf = 2, local_struct = NULL, regular = FALSE, lookup = NULL) {
  
  if (! is.null(local_struct) && local_struct %in% c("tree", "ldag", "part")) {
    
    scorepar <- BiDAG::scoreparameters("usr", 
                                       data = data, 
                                       usrpar = list(pctesttype = "bdecat"))
    
    # additional parameters for user-specified score params 
    scorepar$ess  <- ess 
    scorepar$edgepf <- edgepf
    scorepar$Cvec <- nlev
    scorepar$levels <- lapply(nlev-1, seq.int, from = 0)
    scorepar$local_struct <- local_struct
    scorepar$regular <- regular
    
    # define environment in which to store optimized local structures
    if (!is.null(lookup)) {
      tmp <- list()
      assign(local_struct, tmp, lookup)
      scorepar$lookup <- lookup
      
      usrDAGcorescore <- function(j, parentnodes, n, scorepar) {
        score_from_lookup(scorepar$data,
                          levels = scorepar$levels,
                          nlev = scorepar$Cvec,
                          j,
                          parentnodes,
                          ess = scorepar$ess,
                          method = scorepar$local_struct,
                          regular = scorepar$regular,
                          lookup = scorepar$lookup) - length(parentnodes)*log(scorepar$edgepf)
      }
    } else {
      usrDAGcorescore <- function(j, parentnodes, n, scorepar) {
        c(compute_local_bdeu_score(data, 
                                   scorepar$levels, 
                                   scorepar$Cvec, 
                                   j, parentnodes, 
                                   ess = scorepar$ess, 
                                   method = scorepar$local_struct, 
                                   regular = scorepar$regular)) - length(parentnodes)*log(scorepar$edgepf)
      }
    }
    
    assignInNamespace("usrDAGcorescore", usrDAGcorescore, ns = "BiDAG")
  } else {
    
    # genereate data.frame, removing unobserved levels in data, as required by BiDAG
    df <- data.frame(apply(data, 2, factor, exclude = NULL, simplify = FALSE))
    
    scorepar <- BiDAG::scoreparameters("bdecat", 
                                       data = df, 
                                       bdecatpar = list(chi = ess, 
                                                        edgepf = edgepf))
  }
  return(scorepar)
}