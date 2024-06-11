



# prep ---- 
library(ldags)
library(doSNOW)
library(foreach)
library(BiDAG)
library(ldags)

source("./simulations/run_sample_dags.R")
run <- run_sample_dags

simname <- "LDAG10_2"
outdir  <- paste0("simulations/", simname, "/MCMCchains/")
if (!dir.exists(outdir)) dir.create(outdir)

logdir <- paste0("simulations/", simname, "/MCMCchains_log/")
if (!dir.exists(logdir)) dir.create(logdir)


# params ----
simpar <- expand.grid(list(bnname = c("LDAG10"), 
                           init = c("hcskel"),
                           struct = c("ldag", "dag", "tree"),
                           sample = "partition", 
                           edgepf = c(1, 2, 10**3),
                           ess = 1,
                           hardlimit = 5, 
                           regular = c(FALSE),
                           N = c(300, 1000, 3000, 10000),
                           r = 1:30),
                      stringsAsFactors = F)
indx <- simpar$regular == T & (simpar$edgepf > 1 | !simpar$struct == "tree")
simpar <- simpar[!indx, ]


# test ----
if (FALSE) {
  bnname <- "LDAG10"
  N <- 1000
  r <- 1
  test <- run(simpar[1, ], outdir = NULL, verbose = T)
  str(test, max.level = 2)
  attributes(test)
}

# run ---- 

tag <- paste(simname, format(Sys.time(), "%Y%m%d_%H%M%S"), sep = "-")
filepath <- paste0(logdir, tag, "-session_info.rds")
saveRDS(list(simpar, run, sessionInfo()), filepath)

filepath <- paste0(logdir, tag, ".out")
cl <- makeCluster(6, type = "SOCK", outfile = filepath)
clusterExport(cl, c("run", "simpar", "tag", "outdir"))

registerDoSNOW(cl)
foreach(r = 1:nrow(simpar)) %dopar% run(simpar[r, ],
                                        outdir = outdir,
                                        simid = tag,
                                        verbose = T)
stopCluster(cl)


