
check_convergence_edgep <- function(bn, local_struct, ess, edgepf, hardlimit){
  data <- bida:::sample_data_from_bn(bn, 1000)
  nlev <- sapply(bn, function(x) dim(x$prob)[1])
 
  smpls <- list()
  par(mfrow = c(1, 2))
  
  for (algo in c("order", "partition")) {
    tmp <- replicate(2, ldags::sample_dags(data, nlev, algo = algo,
                                             ess = ess, edgepf = edgepf, hardlimit = hardlimit, 
                                             local_struct = local_struct, verbose = T), simplify = F)
    edgeps <- lapply(tmp, BiDAG::edgep, pdag = T)
    BiDAG::plotpcor(edgep, main = paste0(algo))
    mtext(sprintf("local_struct = %s, ess = %s, edgepf = %s, hardlimit = %s", 
                  local_struct, ess, edgepf, hardlimit), 3)
    smpls[[algo]] <- tmp
  }
  
  return(smpls)
}
