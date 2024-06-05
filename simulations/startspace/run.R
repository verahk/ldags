



library(ldags)
library(doSNOW)
library(BiDAG)
library(ldags)

source("./simulations/startspace/R/define_scorepar.R")
source("./simulations/startspace/R/init_search_space.R")
source("./simulations/startspace/R/sample_dags.R")

simpar <- expand.grid(list(bnname = c("LDAG10"), 
                           init = c("hc", "hcskel", "pcskel"),
                           sample = "partition", 
                           edgepf = c(1, 2, 10**4),
                           regular = c(TRUE, FALSE),
                           N = 1000,
                           r = 1:10, 
                           struct = c("dag", "tree", "ldag")),
                      stringsAsFactors = F)

run <- function(bnname, init, struct, sample, edgepf, regular, N, r, write_to_file = F, verbose = T) {
  
  filename <- paste(bnname, init, struct, sample, 
                    sprintf("epf%s_reg%s_N%s_r%02.0f.rds", edgepf, regular*1, N, r),
                    sep = "_")
  filepath <- paste0("./simulations/startspace/MCMCchains/", filename)
  if (file.exists(filepath) && write_to_file) return(NULL)
  
  if (verbose) cat(filename)
  bn <- readRDS(paste0("./data/", bnname, ".rds"))
  
  # draw data
  set.seed(N+r)
  data <- bida:::sample_data_from_bn(bn, N)
  nlev <- sapply(bn, function(x) dim(x$prob)[1])
  
  # define scorepars
  scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = struct)
  
  # run MCMC
  smpl <- sample_dags(scorepar, init, sample, hardlimit = 5, verbose = verbose)
  
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
  test <- run(bnname, init = "hcskel", struct = "tree", sample = "partition",  edgepf = 1, regular = FALSE, N = 1000, r = 1, write_to_file = F)
}

# run ---- 
# for (r in 1:nrow(simpar)) run(simpar$bnname[r], 
#                               simpar$init[r], 
#                               simpar$struct[r],
#                               simpar$sample[r],
#                               simpar$edgepf[r],
#                               simpar$N[r],
#                               simpar$r[r],
#                               write_to_file = T)

file.remove("tmp.out")
cl <- makeCluster(4, type = "SOCK", outfile = "tmp.out")
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
                                        write_to_file = T)



