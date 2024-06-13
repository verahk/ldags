

devtools::install_github("verahk/ldags")

# prep ---- 
library(ldags)
library(doSNOW)
library(foreach)
library(BiDAG)
library(ldags)

source("./simulations/run_sample_dags.R")
run <- run_sample_dags

simname <- "bn_2"
outdir  <- paste0("simulations/", simname, "/MCMCchains/")
if (!dir.exists(outdir)) dir.create(outdir)

logdir <- paste0("simulations/", simname, "/MCMCchains_log/")
if (!dir.exists(logdir)) dir.create(logdir)

nClusters <- 6

# params ----
simpar <- expand.grid(list(bnname = c("asia", "sachs"), 
                           init = c("hcskel"),
                           struct = c("dag", "tree", "ldag"),
                           sample = "partition", 
                           edgepf = 10**c(0, 1, 2),
                           ess = 1,
                           hardlimit = 4, 
                           regular = c(FALSE, TRUE),
                           N = c(300, 1000, 3000, 10000),
                           r = 1:30),
                      stringsAsFactors = F)

indx <-  with(simpar, 
              (struct == "dag" & !(edgepf == 1 & regular == TRUE)) | (struct != "dag" & (edgepf == 1 | regular == TRUE)))

simpar <- simpar[!indx, ]


# test ----
if (FALSE) {
  N <- 1000
  r <- 1
  test <- run(simpar[1, ], outdir = NULL, verbose = T)
  
  testpar <- simpar[1, ]
  testpar$bnname <- "sachs"
  
  str(test, max.level = 2)
  attributes(test)
}

# run ---- 

tag <- paste(simname, format(Sys.time(), "%Y%m%d_%H%M%S"), sep = "-")
filepath <- paste0(logdir, tag, "-session_info.rds")
saveRDS(list(simpar, run, sessionInfo()), filepath)

filepath <- paste0(logdir, tag, ".out")
cl <- makeCluster(nClusters, type = "SOCK", outfile = filepath)
clusterExport(cl, c("run", "simpar", "tag", "outdir"))

if (nClusters > 1) {
  registerDoSNOW(cl)
  foreach(r = 1:nrow(simpar)) %dopar% run(simpar[r, ],
                                          outdir = outdir,
                                          simid = tag,
                                          verbose = T)
  stopCluster(cl)
  
}

for (r in 1:nrow(simpar)) run(simpar[r, ], outdir, simid = tag, verbose = T)

