#' Title
#'  
#' @param data (matrix)
#'  Each column is assumed to be outcomes of a categorical variable, with 
#'  possible outcomes `0, ..., K-1.`
#' @param nlev (integer vector)
#'  Cardinality of each variable.
#' @param ess 
#' @param edgepf 
#' @param local_struct 
#' @param regular 
#' @param lookup 
#'
#' @return
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
        ldags:::score_from_lookup(scorepar$data,
                                  levels = scorepar$levels,
                                  nlev = scorepar$Cvec,
                                  j,
                                  parentnodes,
                                  ess = scorepar$ess,
                                  method = scorepar$local_struct,
                                  lookup = scorepar$lookup) - length(parentnodes)*log(scorepar$edgepf)
      }
    } else {
      usrDAGcorescore <- function(j, parentnodes, n, scorepar) {
        c(ldags:::compute_local_bdeu_score_from_data(data, 
                                                   scorepar$levels, 
                                                   scorepar$Cvec, 
                                                   j, parentnodes, 
                                                   ess = scorepar$ess, 
                                                   struct = scorepar$local_struct, 
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