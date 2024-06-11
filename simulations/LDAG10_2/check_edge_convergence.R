library(ldags)
library(doSNOW)
library(BiDAG)
library(ldags)


# prep ----
bnname <- "LDAG10"
bn <- readRDS(paste0("data/", bnname, ".rds"))

set.seed(007)
N <- 1000
data <- bida:::sample_data_from_bn(bn, N)
nlev <- sapply(bn, function(x) dim(x$prob)[1])

simpar <- expand.grid(list(init = c("pcskel"),
                           struct = c("dag", "tree", "ldag"),
                           sample = c("partition", "order"),
                           regular = c(TRUE, FALSE)), 
                      stringsAsFactors = F)

# functions ----
run <- function(data, init, struct, sample, regular) {
  cat("\nStarted run with:", init, struct, sample, regular)
  ess <- 1
  hardlimit <- 5
  edgepf <- 2
  
  if (struct == "dag" && regular == FALSE) return(NULL) 
  scorepar <- ldags:::define_scorepar(data, nlev, ess = ess, edgepf = edgepf, local_struct = struct, regular = regular)
  smpls <- replicate(2, ldags:::sample_dags(scorepar, init, sample, hardlimit = hardlimit), simplify = F)
  
  trace <- sapply(smpls, "[[", "trace")
  edgeps_dag <- lapply(smpls, BiDAG::edgep, pdag = F)
  edgeps_pdag <- lapply(smpls, BiDAG::edgep, pdag = T)
  
  filename = paste(bnname, init, struct, sample, regular, ".png", sep = "_")
  png(here::here("./simulations/LDAG10_2/convergence_plots", filename)) 
  
  par(mar = c(2, 2, 0, 1))
  nf <- layout(matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, byrow = T), heights = c(.6, 3))
  layout.show(nf)
  
  title <- paste(init, struct, regular, sample, sep = "+")
  sub <- sprintf("N = %s, ess = %s, edgepf = %s, hardlimit = %s, regular = %s", 
                 N, ess, edgepf, hardlimit, regular)
  plot.new()
  text(.5, .6, title, cex = 1.5, font = 2)
  text(.5, .1, sub, cex = 1.25, font = 2)
  
  matplot(trace, type = "l", xlab = "iteration")
  
  plot(edgeps_dag[[1]], edgeps_dag[[2]])
  abline(a = 0, b = 1, col = "red")
  
  plot(edgeps_pdag[[1]], edgeps_pdag[[2]])
  abline(a = 0, b = 1, col = "red")
  
  dev.off()
  return(NULL)
}


# test ----
run(data, "hc", "dag", "order", T)

# run ----


file.remove("tmp.out")
cl <- makeCluster(4, type = "SOCK", outfile = "tmp.out")
clusterExport(cl, c("data", "run", "simpar", "init_search_space", "define_scorepar", "sample_dags"))
registerDoSNOW(cl)

foreach(r = 1:nrow(simpar)) %dopar% run(data,
                                        simpar$init[r],
                                        simpar$struct[r],
                                        simpar$sample[r],
                                        simpar$regular[r])

# results ----
# collect graphs in one markdown file 
filepath <- "./simulations/LDAG10_2/convergence_plots/"
filename <- paste0(filepath, "plots.md")
cat("# Convergence diagnostics for simulations: LDAG10_2\n", file = filename, append = F)
files <- list.files(filepath, full.names = F)

for (m in c("dag", "ldag", "tree")) {
  cat("\n##", toupper(m), "\n", file = filename, append = T)
  for (reg in c("TRUE", "FALSE")) {
    indx <- grepl(paste("", m, "", sep = "_"), files) & grepl(paste("", reg, "", sep = "_"), files)
    for (f in files[indx]) cat(sprintf("![%s](%s)\n", f, f), file = filename, append = T)
  }
}
