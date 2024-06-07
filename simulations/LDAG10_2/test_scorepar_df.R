


bn <- readRDS("./data/LDAG10.rds")
data <- bida:::sample_data_from_bn(bn, 100)
n <- length(bn)
ess <- 1
edgepf <- 2 
df <- data.frame(apply(data, 2, factor, exclude = NULL, simplify = FALSE))

scorepar1 <- BiDAG::scoreparameters("bdecat", 
                                   data = df, 
                                   bdecatpar = list(chi = ess, 
                                                    edgepf = edgepf))
scorepar2 <- BiDAG::scoreparameters("bdecat", 
                                    data = as.data.frame(df), 
                                    bdecatpar = list(chi = ess, 
                                                     edgepf = edgepf))

for (j in seq_along(bn)) {
  parentnodes <- which(bnlearn::amat(bn)[, j] == 1)
  score1 <- BiDAG:::DAGcorescore(j, parentnodes, n, scorepar1)
  score2 <- BiDAG:::DAGcorescore(j, parentnodes, n, scorepar2)
  stopifnot(score1 == score2)
}
  