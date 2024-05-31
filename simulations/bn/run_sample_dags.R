
# Simulation: Sample DAGs with local-structure using BiDAG


library(doSNOW)
library(BiDAG)
library(ldags)

# prep ----
bnname <- "child"
nClusters <- 5
simpar <- expand.grid(list(method = c( "tree", "ldag", "dag"),
                           N = c(300, 1000, 3000, 10000),
                           r = 20:30))
simpar$method <- as.character(simpar$method)
doRun <- T

run <- function(bnname, r, N, m, algo = "order", write_to_file = T, verbose = T)  {
  
  bn <- readRDS(paste0("./data/", bnname, ".rds"))
  n  <- length(bn)
  edgepf <- 4
  
  # sample data
  set.seed(r+N)
  data <- bida:::sample_data_from_bn(bn, N)
  nlev <- sapply(bn, function(x) dim(x$prob)[1])
  
  filename <- sprintf("%s_%s_%s_N%s_epf%s_r%02.0f.rds", bnname, algo, m, N, edgepf, r)
  filepath <- here::here("./simulations/bn/MCMCchains/", filename)
  if (file.exists(filepath) && write_to_file) return(NULL)
  cat(filename)
  
  lookup <- rlang::new_environment()
  smpl <- ldags::sample_dags(data, nlev, algo, 
                             ess = 1, edgepf = edgepf, hardlimit = 4, 
                             local_struct = m, lookup = lookup, verbose = verbose)
  out <- list(smpl = smpl, lookup = lookup)
  if (write_to_file) {
    saveRDS(out, filepath)
  } else {
    return(out)
  }
}


# test ----
if (FALSE) {
  
  res <- run(bnname, 1, 1000, "tree", "order", T, T)
  res <- run(bnname, 1, 100, "ldag", "order", F, T)
  test <- readRDS(here::here("./simulations/bn/MCMCchains/sachs_order_tree_N1000_r02.rds"))

}

# run ----
if (doRun) {

  if (nClusters > 1) {
    cl <- makeCluster(4, outfile = "")
    registerDoSNOW(cl)
    foreach(r = simpar$r, N = simpar$N, m = simpar$method) %dopar% run(bnname, r, N, m, "order", T, F)
   # foreach(r = simpar$r, N = simpar$N, m = simpar$method,
  #          .packages = c("ldags", "BiDAG")) %dopar% run(bn, r, N, m, "partition", T, T)
    stopCluster(cl)
  } else {
    for (i in 1:nrow(simpar)) {
      run(bnname, simpar$r[i], simpar$N[i], simpar$method[i], algo = "order", T, T)
    }
  }
}

