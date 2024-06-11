
library(ggplot2)

bnname <- "LDAG10"
bn <- readRDS(paste0("./data/", bnname, ".rds"))
dag <- bnlearn::amat(bn)
N <- 1000
r <- 1

set.seed(N+r)
data <- bida:::sample_data_from_bn(bn, N)
nlev <- sapply(bn, function(x) dim(x$prob)[1])
levels <- lapply(nlev-1, seq.int, from = 0)

# REGULARITY ----
j <- 6
parentnodes <- 2:4
local_struct <- "tree"

regular  <- T 
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = local_struct, regular = regular)
score <- BiDAG:::usrDAGcorescore(j, parentnodes, length(bn), scorepar)
famscore <- ldags:::compute_local_bdeu_score_from_data(data, levels, nlev, j, parentnodes, ess = 1, struct = local_struct, regular = regular)
abs(score - (famscore - length(parentnodes)*log(edgepf))) < 10**-10
cat("score:", score, "famscore:", famscore, "famscore+pen:",  (famscore - length(parentnodes)*log(edgepf)))

regular  <- F
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = local_struct, regular = regular)
score <- BiDAG:::usrDAGcorescore(j, parentnodes, length(bn), scorepar)
famscore <- ldags:::compute_local_bdeu_score_from_data(data, levels, nlev, j, parentnodes, ess = 1, struct = local_struct, regular = regular)
abs(score - (famscore - length(parentnodes)*log(edgepf))) < 10**-10
cat("score:", score, "famscore:", famscore, "famscore+pen:",  (famscore - length(parentnodes)*log(edgepf)))


# EDGE PENALTY ---- 
# - test that score with edgepenality is correctly computed
j <- 7
parentnodes <- which(dag[, j] == 1)
edgepf <- 2 

local_struct <- "tree"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = local_struct)
score <- BiDAG:::usrDAGcorescore(j, parentnodes, length(bn), scorepar)
famscore <- ldags:::compute_local_bdeu_score_from_data(data, levels, nlev, j, parentnodes, ess = 1, struct = local_struct)
abs(score - (famscore - length(parentnodes)*log(edgepf))) < 10**-10
cat("score:", score, "famscore:", famscore, "famscore+pen:",  (famscore - length(parentnodes)*log(edgepf)))

local_struct <- "ldag"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = local_struct)
score <- BiDAG:::usrDAGcorescore(j, parentnodes, length(bn), scorepar)
famscore <- ldags:::compute_local_bdeu_score_from_data(data, levels, nlev, j, parentnodes, ess = 1, struct = local_struct)
abs(score - (famscore - length(parentnodes)*log(edgepf))) < 10**-10
cat("score:", score, "famscore:", famscore, "famscore+pen:",  (famscore - length(parentnodes)*log(edgepf)))


local_struct <- "dag"
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = local_struct)
score <- BiDAG:::DAGcorescore(j, parentnodes, length(bn), scorepar)
famscore <- ldags:::compute_local_bdeu_score_from_data(data, levels, nlev, j, parentnodes, ess = 1, struct = NULL)
abs(score - (famscore - length(parentnodes)*log(edgepf))) < 10**-10
cat("score:", score, "famscore:", famscore, "famscore+pen:",  (famscore - length(parentnodes)*log(edgepf)))


# SCORE TRUE DAG ---- 
seqn <- seq_len(ncol(dag))
parentnodes <- apply(dag, 2, function(x) seqn[x == 1])
famscore <- function(j, pa, struct, data) {
  ldags:::compute_local_bdeu_score_from_data(data, levels, nlev, j, pa, ess = 1, struct = switch(struct, "dag" = NULL, struct))
}
famscores <- function(struct, data) {
  mapply(famscore,  j = seqn, pa = parentnodes, MoreArgs = list(struct = struct, data = data))
}

# no edge penalization
tmp <- sapply(c("tree", "ldag", "dag"), famscores, data = data)
colSums(tmp)

# with edge penalization
tmp <- tmp + sapply(parentnodes, length)*log(2)
colSums(tmp)

test <- function() {
  data <- bida:::sample_data_from_bn(bn, N)
  tmp  <- sapply(c("tree", "ldag", "dag"), famscores, data = data) 
  cbind(epf1 = colSums(tmp),
        epf2 = colSums(tmp + sapply(parentnodes, length)*log(2)), 
        epf1000 = colSums(tmp - sapply(parentnodes, length)*log(1000)))
}

tmp <- replicate(10, test())
df <- setNames(reshape2::melt(tmp), c("struct", "epf", "iter", "value"))
df %>%   
  tidyr::pivot_wider(names_from = "struct", values_from = "value") %>% 
  tidyr::pivot_longer(c("tree", "ldag"), names_to = "struct") %>% 
ggplot(aes(dag, value, color = struct, group = interaction(epf, iter))) +
  facet_grid(.~epf) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() 

# SCORE OF MAP DAG -----
bnname <- "LDAG10"
bn <- readRDS(paste0("./data/", bnname, ".rds"))
dag <- bnlearn::amat(bn)
N <- 1000
r <- 1

# sample data
set.seed(N+r)
bn <- ldags:::example_bn(bnname)
data <- bida:::sample_data_from_bn(bn, N)
nlev <- sapply(bn, function(x) dim(x$prob)[1])
levels <- lapply(nlev-1, seq.int, from = 0)
run <- run(bnname, init = "hcskel", struct = "ldag", sample = "partition", 
            N = N, r = 1, edgepf = 2, regular = F, write_to_file = F, verbose = T)
scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = 2, "ldag", F)
scores <- numeric(length(bn))
for (j in seq_along(bn)) {
  parentnodes <- which(test$DAG[, j] == 1)
  famscore <- compute_local_bdeu_score_from_data(data, levels, nlev, j, parentnodes, ess = 1, struct = "ldag")
  scores[j] <- famscore - length(parentnodes)*log(2)
  stopifnot(scores[j] == BiDAG:::usrDAGcorescore(j, parentnodes, n, scorepar))
}
sum(scores) -test$score

for (struct in c("ldag", "dag")) {
  set.seed(N+r)
  data <- bida:::sample_data_from_bn(bn, 1000)
  scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = 2, local_struct = struct)
  smpl <- sample_dags(scorepar, "pcskel", "order", 5, T)
  cat("\nLocal structure:", struct, "score: ", smpl$score)
  plot(smpl)
}






