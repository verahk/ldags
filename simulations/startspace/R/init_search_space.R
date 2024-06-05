

init_search_space <- function(scorepar, algo, hardlimit, maxp = hardlimit, alpha = .05, verbose) {

  if (algo == "hc") {
    df <- data.frame(apply(scorepar$data, 2, factor, exclude = NULL, simplify = FALSE))
    fit  <- bnlearn:::hc(df, maxp = maxp)
    startspace <- bnlearn::amat(fit)
  } else if (algo == "hcskel") {
    df <- data.frame(apply(scorepar$data, 2, factor, exclude = NULL, simplify = FALSE))
    fit  <- bnlearn:::hc(df, maxp = maxp)
    startspace <- bnlearn::amat(fit)
    startspace <- startspace + t(startspace)
  } else if (algo == "pc") {
    suffStat <- list(dm = scorepar$data, nlev = scorepar$Cvec, adaptDF = FALSE)
    fit <- pcalg::pc(suffStat, pcalg::disCItest, alpha = alpha, labels = colnames(scorepar$data), verbose = verbose)
    startspace <- as(fit@graph, "matrix")
  } else if (algo == "pcskel") {
    suffStat <- list(dm = scorepar$data, nlev = scorepar$Cvec, adaptDF = FALSE)
    fit <- pcalg::skeleton(suffStat, pcalg::disCItest, alpha = alpha, labels = colnames(scorepar$data), verbose = verbose)
    startspace <- as(fit@graph, "matrix")
  } else if (algo == "BiDAG") {
    fit <- try(iterativeMCMC(scorepar, hardlimit = hardlimit, alpha = alpha))
  }

  if (! any(colSums(startspace) > hardlimit)) {
    return(startspace)
  } else if (maxp > 1 && alpha > .001) {
    init_search_space(scorepar, algo, hardlimit, maxp = maxp-1, alpha = alpha/2, verbose)
  } else {
    cat("\nCould not find startspace satisfying `hardlimit`\n")
    cat(sprintf("algo: %s hardlimit: %s maxp: %s alpha: %s", algo, hardlimit, maxp, alpha))
    stop()
  }
}

