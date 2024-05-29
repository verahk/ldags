
# example for op

# example frequency counts for a CPT with 2**3 rows 
counts <- cbind(1000, c(rep(1000, 4), rep(100, 2), 10, 1))
lev  <- rep(list(0:1), 3)
ess <- 1
res <- optimize_partition(counts, lev, "tree", ess = ess, verbose = TRUE)
res$partition
res$tree
cbind(counts, unlist_partition(res$partition))

gr <- unlist_partition(res$partition)
tmp <- rowsum(counts, gr)
scores <- famscore_bdeu_byrow(tmp, ess, r = ncol(counts), q = nrow(counts), s = lengths(res$partition))
all.equal(scores, res$scores)

# profile 
if (FALSE) {
  
  set.seed(007)
  lev    <- rep(list(0:1), 8)
  counts <- cbind(100, sample(prod(lengths(lev)))*10)
  res <- optimize_partition(counts, lev, "tree", ess = 1, TRUE)
  stopifnot(sum(lengths(res$partition)) == nrow(counts))
  cbind(counts, unlist_partition(res$partition))
  #tmp <- do.call("rbind", res$partition)
  #stopifnot(all(nrow(tmp) == nrow(counts)))
  
 
  microbenchmark::microbenchmark(f(tree), f2(tree))
  
  profvis::profvis(optimize_partition(counts, lev, "tree", ess = 1))
}
