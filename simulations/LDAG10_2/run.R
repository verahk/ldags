



library(ldags)
library(doSNOW)
library(foreach)
library(BiDAG)
library(ldags)

simpar <- expand.grid(list(bnname = c("LDAG10"), 
                           init = c("pcskel"),
                           struct = c("ldag", "dag", "tree"),
                           sample = "partition", 
                           edgepf = c(1, 2, 10**3),
                           regular = c(FALSE),
                           N = c(300, 1000, 3000, 10000),
                           r = 1:30),
                      stringsAsFactors = F)
indx <- simpar$regular == T & (simpar$edgepf > 1 | !simpar$struct == "tree")
simpar <- simpar[!indx, ]

run <- function(bnname, init, struct, sample, edgepf, regular, N, r, write_to_file = F, verbose = F) {
  
  filename <- paste(bnname, init, struct, sample, 
                    sprintf("epf%s_reg%s_N%s_r%02.0f.rds", edgepf, regular*1, N, r),
                    sep = "_")
  filepath <- paste0("./simulations/LDAG10_2/MCMCchains/", filename)
  if (file.exists(filepath) && write_to_file) return(NULL)
  
  if (verbose) cat(filename)
  #bn <- readRDS(paste0("./data/", bnname, ".rds"))
  
  # draw data
  set.seed(N+r)
  bn <- ldags:::example_bn(bnname)
  data <- bida:::sample_data_from_bn(bn, N)
  nlev <- sapply(bn, function(x) dim(x$prob)[1])
  
  # define scorepars
  scorepar <- ldags:::define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = struct)
  
  # run MCMC
  smpl <- ldags:::sample_dags(scorepar, init, sample, hardlimit = 4, verbose = verbose)
  
  if (write_to_file) {
    saveRDS(smpl, filepath)
    return(NULL)
  } else {
    return(smpl)
  }
}


# test 
if (FALSE) {
  bnname <- "LDAG10"
  test <- run(bnname, init = "hcskel", struct = "tree", sample = "partition",  
              edgepf = 1, regular = FALSE, N = 1000, r = 1, write_to_file = F, verbose = T)
}

# run ---- 
for (r in 1:nrow(simpar)) run(simpar$bnname[r], 
                              simpar$init[r], 
                              simpar$struct[r],
                              simpar$sample[r],
                              simpar$edgepf[r],
                              simpar$regular[r],
                              simpar$N[r],
                              simpar$r[r],
                              write_to_file = T)

file.remove("tmp2.out")
cl <- makeCluster(6, type = "SOCK", outfile = "")
clusterExport(cl, c("run", "simpar", "init_search_space", "define_scorepar", "sample_dags"))
registerDoSNOW(cl)

foreach(r = 1:nrow(simpar)) %dopar% run(simpar$bnname[r], 
                                        simpar$init[r], 
                                        simpar$struct[r],
                                        simpar$sample[r],
                                        simpar$edgepf[r],
                                        simpar$regular[r],
                                        simpar$N[r],
                                        simpar$r[r],
                                        write_to_file = T,
                                        verbose = T)
stopCluster(cl)


