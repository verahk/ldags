
# Check convergence of MCMC-chain of DAGs
# - For each local-struct and MCMC-scheme, run 2 MCMC-runs and compare the edge-weights
bn <- readRDS("./data/child.rds")
N <- 1000
ess <- 1
edgepf <- .25
hardlimit <- 4
set.seed(007)
data <- bida:::sample_data_from_bn(bn, N)
nlev <- sapply(bn, function(x) dim(x$prob)[1])

for (local_struct in c("tree", "ldag", "dag")) {
  for (algo in c("order", "partition")) {
    cat("\n Local structure:", local_struct, "Algo:", algo)
    tic <- Sys.time()
    smpls <- replicate(2, ldags::sample_dags(data, nlev, algo = algo,
                                           ess = ess, edgepf = edgepf, hardlimit = hardlimit, 
                                           local_struct = local_struct, verbose = T), simplify = F)
    print(Sys.time()-tic)
    
    trace <- sapply(smpls, "[[", "trace")
    edgeps_dag <- lapply(smpls, BiDAG::edgep, pdag = F)
    edgeps_pdag <- lapply(smpls, BiDAG::edgep, pdag = T)
    
    png(filename = sprintf("./simulations/bn/edgeps_child_%s_%s.png", local_struct, algo))
    par(mar = c(2, 2, 0, 1))
    nf <- layout(matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, byrow = T), heights = c(.6, 3))
    layout.show(nf)
    
    title <- paste0(local_struct, " + ", algo)
    sub <- sprintf("N = %s, ess = %s, edgepf = %s, hardlimit = %s", 
                   N, ess, edgepf, hardlimit)
    plot.new()
    text(.5, .6, title, cex = 1.5, font = 2)
    text(.5, .1, sub, cex = 1.25, font = 2)

    matplot(trace, type = "l", xlab = "iteration")
   
    plot(edgeps_dag[[1]], edgeps_dag[[2]])
    abline(a = 0, b = 1, col = "red")

    plot(edgeps_pdag[[1]], edgeps_pdag[[2]])
    abline(a = 0, b = 1, col = "red")
    
    dev.off()
  }
}
