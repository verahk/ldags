test_that("sample_dags handles multi-level bn", {
  
  bnname <- "sachs"
  bn <- readRDS(paste0("./data/", bnname, ".rds"))
  nlev <- 
  
  # draw data 
  N <- 1000
  data <- bida:::sample_data_from_bn(bn, N)
  nlev <- sapply(bn, function(x) dim(x$prob)[1])
  
  # define scorepars
  struct <- "ldag"
  regular <- T
  scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = 2, local_struct = struct, regular = regular)
  
  # run MCMC
  init   <- "hcskel"
  sample <- "partition"
  smpl <- sample_dags(scorepar, init, sample, hardlimit = 5, verbose = T)
  
})
