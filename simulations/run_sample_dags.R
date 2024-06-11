

run_sample_dags <- function(simpar, outdir = NULL, simid = "", verbose = F) {
  
  args <- as.list(environment())
  for (arg in names(simpar)) assign(arg, simpar[[arg]])
    
  filename <- paste(bnname, init, struct, sample, 
                    sprintf("epf%s_reg%s_N%s_r%02.0f.rds", edgepf, regular*1, N, r),
                    sep = "_")
  filepath <- paste0(outdir, filename)
  if (!is.null(outdir) && file.exists(filepath)) return(NULL)
  
  if (verbose) cat(filename)

  # load BN
  if (bnname == "LDAG10") {
    set.seed(N+r)
    bn <- ldags:::example_bn(bnname)
  } else {
    bn <- readRDS(paste0("./data/", bnname, ".rds"))
  }
  
  # draw data 
  set.seed(N+r)
  data <- bida:::sample_data_from_bn(bn, N)
  nlev <- sapply(bn, function(x) dim(x$prob)[1])
  
  # define scorepars
  scorepar <- ldags:::define_scorepar(data, nlev, ess = 1, edgepf = edgepf, local_struct = struct)
  
  # run MCMC
  smpl <- ldags:::sample_dags(scorepar, init, sample, hardlimit = 4, verbose = verbose)
  attr(smpl, "args") <- args 
  
  if (!is.null(outdir)) {
    saveRDS(smpl, filepath)
    cat("output saved to ", filepath)
    return(NULL)
  } else {
    return(smpl)
  }
}
