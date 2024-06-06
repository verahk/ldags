
library(ldags)
library(doSNOW)
library(BiDAG)
library(ldags)

source("./simulations/startspace/R/define_scorepar.R")
source("./simulations/startspace/R/init_search_space.R")
source("./simulations/startspace/R/sample_dags.R")

bn <- readRDS(paste0("./data/", bnname, ".rds"))

# params 
N <- 1000
r <- 1
ess <- 1 
edgepf <- 2
hardlimit <- 4 
algo_init <- "pcskel"
algo_sample <- "partition"


# draw data
set.seed(N+r)
bn <- ldags:::example_bn(bnname)
data <- bida:::sample_data_from_bn(bn, N)
nlev <- sapply(bn, function(x) dim(x$prob)[1])


# tree ----
struct <- "tree"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = struct)
profvis::profvis(sample_dags(scorepar, algo_init, algo_sample, hardlimit = hardlimit))

# ldag ----
struct <- "ldag"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = struct)
profvis::profvis(sample_dags(scorepar, algo_init, algo_sample, hardlimit = hardlimit, verbose = T))

# part ----
struct <- "part"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = struct)
profvis::profvis(sample_dags(scorepar, algo_init, algo_sample, hardlimit = hardlimit, verbose = T))

