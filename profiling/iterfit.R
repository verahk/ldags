

library(ldags)
library(doSNOW)
library(BiDAG)
library(ldags)

# params 
N <- 1000
r <- 1
ess <- 1 
edgepf <- 2
hardlimit <- 5 
algo_init <- "hcskel"
algo_sample <- "partition"
verbose <- T

# draw data
set.seed(N+r)
bn <- ldags:::example_bn("LDAG10")
data <- bida:::sample_data_from_bn(bn, N)
nlev <- sapply(bn, function(x) dim(x$prob)[1])


# tree ----
struct <- "tree"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = struct)
startspace <- ldags:::init_search_space(scorepar, "hcskel", hardlimit = hardlimit, verbose = T)
profvis::profvis(BiDAG::learnBN(scorepar, "orderIter", hardlimit = hardlimit, startspace = startspace, scoreout = T, verbose = verbose))


# ldag ----
struct <- "ldag"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = struct)
startspace <- ldags:::init_search_space(scorepar, "hcskel", hardlimit = hardlimit, verbose = T)
profvis::profvis(BiDAG::learnBN(scorepar, "orderIter", hardlimit = hardlimit, startspace = startspace, scoreout = T, verbose = verbose))


# part ----
struct <- "part"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = struct)
startspace <- ldags:::init_search_space(scorepar, "hcskel", hardlimit = hardlimit, verbose = T)
profvis::profvis(BiDAG::learnBN(scorepar, "orderIter", hardlimit = hardlimit, startspace = startspace, scoreout = T, verbose = verbose))

# dag ----
struct <- "dag"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = struct)
startspace <- ldags:::init_search_space(scorepar, "hcskel", hardlimit = hardlimit, verbose = T)
profvis::profvis(BiDAG::learnBN(scorepar, "orderIter", hardlimit = hardlimit, startspace = startspace, scoreout = T, verbose = verbose))


# sachs ----
bnname <- "sachs"
bn <- readRDS(paste0("./data/", bnname, ".rds"))
N <- 1000
r <- 1

# draw data 
set.seed(N+r)
data <- bida:::sample_data_from_bn(bn, N)
nlev <- sapply(bn, function(x) dim(x$prob)[1])

# define scorepars
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = 1000, local_struct = "ldag", regular = T)
startspace <- init_search_space(scorepar, "hcskel", hardlimit = 4, verbose = T)
profvis::profvis(BiDAG::learnBN(scorepar, "orderIter", hardlimit = 4, startspace = startspace, scoreout = T, verbose = T))

