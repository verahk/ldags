
bn <- readRDS(paste0("./data/", bnname, ".rds"))
dag <- bnlearn::amat(bn)
N <- 1000
r <- 1

set.seed(N+r)
data <- bida:::sample_data_from_bn(bn, N)
nlev <- sapply(bn, function(x) dim(x$prob)[1])
levels <- lapply(nlev-1, seq.int, from = 0)

# REGULARITY ----
j <- 1
parentnodes <- 2:5
local_struct <- "tree"
regular  <- T 
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = local_struct, regular = regular)
famscore <- ldags:::compute_local_bdeu_score_from_data(data, levels, nlev, j, parentnodes, ess = 1, regular = regular, local_struct = NULL)
score <- BiDAG:::usrDAGcorescore(j, parentnodes, length(bn), scorepar)
abs(score - (famscore - length(parentnodes)*log(edgepf))) < 10**-10
cat("score:", score, "famscore:", famscore, "famscore+pen:",  (famscore - length(parentnodes)*log(edgepf)))


# EDGE PENALTY ---- 
# - test that score with edgepenality is correctly computed
j <- 7
parentnodes <- which(dag[, j] == 1)
edgepf <- 2 

local_struct <- "tree"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = local_struct)
famscore <- compute_local_bdeu_score_from_data(data, levels, nlev, j, parentnodes, ess = 1, local_struct = NULL)
score <- BiDAG:::usrDAGcorescore(j, parentnodes, length(bn), scorepar)
abs(score - (famscore - length(parentnodes)*log(edgepf))) < 10**-10
cat("score:", score, "famscore:", famscore, "famscore+pen:",  (famscore - length(parentnodes)*log(edgepf)))

local_struct <- "ldag"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = local_struct)
famscore <- compute_local_bdeu_score_from_data(data, levels, nlev, j, parentnodes, ess = 1, local_struct = NULL)
score <- BiDAG:::usrDAGcorescore(j, parentnodes, length(bn), scorepar)
abs(score - (famscore - length(parentnodes)*log(edgepf))) < 10**-10
cat("score:", score, "famscore:", famscore, "famscore+pen:",  (famscore - length(parentnodes)*log(edgepf)))


local_struct <- "dag"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = local_struct)
famscore <- compute_local_bdeu_score_from_data(data, levels, nlev, j, parentnodes, ess = 1, local_struct = NULL)
score <- BiDAG:::DAGcorescore(j, parentnodes, length(bn), scorepar)
abs(score - (famscore - length(parentnodes)*log(edgepf))) < 10**-10
cat("score:", score, "famscore:", famscore, "famscore+pen:",  (famscore - length(parentnodes)*log(edgepf)))



# SCORE OF MAP DAG -----
bnname <- "LDAG10"
bn <- readRDS(paste0("./data/", bnname, ".rds"))
dag <- bnlearn::amat(bn)
N <- 1000
r <- 1

# sample data
set.seed(N+r)
data <- bida:::sample_data_from_bn(bn, N)
nlev <- sapply(bn, function(x) dim(x$prob)[1])
levels <- lapply(nlev-1, seq.int, from = 0)
test <- run(bnname, init = "hcskel", struct = "tree", smpl = "partition", N = 1000, edgepf = 2, write_to_file = F)

scores <- numeric(length(bn))
for (j in seq_along(bn)) {
  parentnodes <- which(test$DAG[, j] == 1)
  famscore <- compute_local_bdeu_score_from_data(data, levels, nlev, j, parentnodes, ess = 1, local_struct = "tree")
  scores[j] <- famscore - length(parentnodes)*log(2)
}
sum(scores) -test$score

# INIT SEARCH SPACE ----
set.seed(N+r)
data <- bida:::sample_data_from_bn(bn, N)
nlev <- sapply(bn, function(x) dim(x$prob)[1])
levels <- lapply(nlev-1, seq.int, from = 0)



