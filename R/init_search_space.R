

#' Initialize search space
#' 
#' Wrapper around structural learning routines ([bnlearn::hc] and [pcalg::pc]),
#' intended to learn a start-space for the [Bidag::iterativeMCMC] procedure. 
#' Forces the maximal number of parents (`hardlimit`) to be respected, 
#' by re-running the structure learning procedure with stricter add-edge-policies
#' (`maxp` and `alpha`, respectively) until a structure that satifies the `hardlimit`
#' criteria is inferred.
#' 
#' @inheritParams sample_dags
#' @param maxp (integer)
#'  Maximal number of parents in CPDAG learned by [bnlearn::hc]. 
#'  Ignored if not `algo_int=="hc"` or `algo_init == "hcskel"`.
#' @param alpha (numeric)
#'  Significance level for the conditional independencies test in [pcalg::pc].
#'  Ignored if not `algo_int=="pc"` or `algo_init == "pcskel"`. 
#' 
#' @return a `n-by-n` adjacency matrix
#' @export
#'
#' @examples
#' @seealso \link{sample_dags}
init_search_space <- function(scorepar, algo_init, hardlimit, maxp = hardlimit, alpha = .05, verbose) {

  if (algo_init == "hc") {
    df <- data.frame(apply(scorepar$data, 2, factor, exclude = NULL, simplify = FALSE))
    fit  <- bnlearn:::hc(df, maxp = maxp)
    startspace <- bnlearn::amat(fit)
  } else if (algo_init == "hcskel") {
    df <- data.frame(apply(scorepar$data, 2, factor, exclude = NULL, simplify = FALSE))
    fit  <- bnlearn:::hc(df, maxp = maxp)
    startspace <- bnlearn::amat(fit)
    startspace <- startspace + t(startspace)
  } else if (algo_init == "pc") {
    suffStat <- list(dm = scorepar$data, nlev = scorepar$Cvec, adaptDF = FALSE)
    fit <- pcalg::pc(suffStat, pcalg::disCItest, alpha = alpha, 
                     labels = colnames(scorepar$data), verbose = verbose)
    startspace <- as(fit@graph, "matrix")
  } else if (algo_init == "pcskel") {
    suffStat <- list(dm = scorepar$data, nlev = scorepar$Cvec, adaptDF = FALSE)
    fit <- pcalg::skeleton(suffStat, pcalg::disCItest, alpha = alpha, 
                           labels = colnames(scorepar$data), verbose = verbose)
    startspace <- as(fit@graph, "matrix")
  } else if (algo_init == "BiDAG") {
    fit <- try(iterativeMCMC(scorepar, hardlimit = hardlimit, alpha = alpha))
  }

  if (! any(colSums(startspace) > hardlimit)) {
    return(startspace)
  } else if (switch(substr(algo_init, 1, 2), "hc" = maxp > 1, "pc" = alpha > 10**-10)) {
    init_search_space(scorepar, algo_init, hardlimit, maxp = maxp-1, alpha = alpha/2, verbose)
  } else {
    cat("\nCould not find startspace satisfying `hardlimit`\n")
    cat(sprintf("algo_init: %s hardlimit: %s maxp: %s alpha: %s", algo_init, hardlimit, maxp, alpha))
    return(startspace*0)
  }
}

