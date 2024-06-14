
library(doSNOW)
library(foreach)
library(dplyr)
library(ggplot2)

simname <- "bn_2"
bnnames <- c("asia", "sachs", "child")

indir  <- paste0("./simulations/", simname, "/MCMCchains")
outdir <- paste0("./simulations/", simname, "/results/")
names  <- c("network", "init", "struct", "sample", "epf",  "regular", "N", "r")

dir.create(outdir)

eval_files <- function(filepaths, bnname, prmethod, verbose = F) {
  
  indx <- grepl(bnname, filepaths, T)
  if (!all(indx)) stop("All files do not match bnname.")
  
  bn <- readRDS(here::here("./data/", paste0(bnname, ".rds")))
  n  <-  length(bn)
  dag  <- bnlearn::amat(bn)
  dmat <- bida:::descendants(bn)
  
  foreach(f = filepaths, 
          .combine  = "rbind",
          .packages = "Matrix") %dopar% eval_MCMCchain(f, dag, dmat, n, prmethod, seq_len(200), verbose = verbose)

}

eval_MCMCchain <- function(f, dag, dmat, n, prmethod, burninsamples, verbose = T) {
  if (verbose) cat(f)
  MCMCchain  <- readRDS(f)
  if ("smpl" %in% names(MCMCchain))  MCMCchain <- MCMCchain$smpl
 
  dindx <- diag(n) == 1
  dags <- lapply(MCMCchain$traceadd$incidence[-burninsamples], as.matrix)
  
  # list unique DAGs
  u    <- unique(dags)
  support  <- bida:::rowsum_fast(rep(1/length(dags), length(dags)), dags, u)
  dags <- u 
  dmats <- lapply(dags, bida:::descendants)
  
  # compute posterior edge probs
  edgep <- Reduce("+", Map("*", dags, support))[!dindx]
  ancp  <- Reduce("+", Map("*", dmats, support))[!dindx]
  
  res <- matrix(nrow = 2, ncol = 3)
  colnames(res) <- c("TPR", "FPR", "avgppv")
  rownames(res) <- c("edgep", "ancp")
  
  # edge probs
  n1 <- sum(dag)

  res[1, 1] <- sum(edgep[dag[!dindx] == 1])/n1
  res[1, 2] <- sum(edgep[dag[!dindx] == 0])/(n**2-n-n1)
  res[1, 3] <- ldags:::compute_avgppv(edgep, dag[!dindx], method = prmethod)
  
  # arp
  pr <- ldags:::compute_prec_recall(ancp, dmat[!dindx], method = prmethod)
  
  n1 <- sum(dmat)-n
  res[2, 1] <- sum(ancp[dmat[!dindx] == 1])/n1
  res[2, 2] <- sum(ancp[dmat[!dindx] == 0])/(n**2-n-n1)
  #res[2, 3] <- compute_avgppv(ancp, dmat)
  res[2, 3] <- pr$avgppv 
  
  
  # parent set size 
  ps <- bida:::parent_support_from_dags(dags, support)
  wmean <- mapply(function(m, w) ncol(m)-rowSums(is.na(m))%*%w,
                  m = ps$sets, 
                  w = ps$support)
  #stopifnot(abs(wmean[5] - sum(sapply(dags, function(g) sum(g[, 5]))*support)) < 10**-10)
  
  list(edgep = edgep, 
       rates = res, 
       pr = pr$df,
       npar = wmean,
       toc = as.matrix(attr(MCMCchain, "toc")))
  
}

res_to_df <- function(res, filepaths, names = c("network", "init", "struct", "sample", "epf", "regular", "N", 
"r")) {
  
  tmp <- stringr::str_split(filepaths, ".+MCMCchains/|_|.rds", simplify = F)
  par <- data.frame(do.call(rbind, tmp)[, seq_along(names)+1])
  colnames(par) <- names
  char_to_factor <- function(x, sub) {
    u <- unique(x)
    ordr <- order(as.numeric(gsub(sub, "", u)))
    factor(x, u[ordr])
  }
  par$N <- char_to_factor(par$N, "N")
  par$epf <- char_to_factor(par$epf, "epf")

  tmp <- do.call(rbind, res)
  if (nrow(tmp) == length(filepaths)) {
    return(cbind(par, tmp))
  } else {
    nrows <- vapply(res, nrow, integer(1))
    df <- cbind(par[rep.int(seq_along(filepaths), nrows), ], tmp)
    df$name <- rownames(tmp)
    return(df)
  }
}


plot_box_plot <- function(df, x, y, title = "", fill = "struct", facets = "algo+name~network+epf") {
  plot <- ggplot(df, aes(x = !! sym(x), y = !! sym(y), fill = !! sym(fill))) +
    facet_grid(as.formula(facets), scales = "free", space = "free") +
    geom_boxplot() +
    ggtitle(title) +
    theme(legend.position = "bottom")
  
  if (fill == "struct") plot <- plot + scale_fill_manual(values = c(dag = "green", ldag = "red", tree = "blue"))
  return(plot)
}
plot_prec_recall <- function(df, title = "", facets = "algo+struct+epf~network+N") {
  rate <- df$PPV[nrow(df)]
  ggplot(df, aes(TPR, PPV, group = r)) +
    facet_grid(as.formula(facets)) +
    geom_line() +
    geom_hline(yintercept = rate, color = "red") +
    ggtitle(title)
}
plot_prec_recall_step <- function(df, title = "", facets = "algo+struct+epf~network+N") {
  rate <- df$PPV[nrow(df)]
  ggplot(df, aes(TPR, PPV, group = r)) +
    facet_grid(as.formula(facets)) +
    geom_hline(yintercept = rate, color = "red") +
    geom_step(direction = "vh") +
    geom_point(data = dplyr::filter(df, !is.infinite(df$x)),
               color = "blue", size = .5) +
    ggtitle(title)
}



plots_to_file <- function(res, filepaths, names, tag, outdir, outfile,  print = TRUE) {

  
  plotter <- function(plot, filepath, outfile, print, ...) {
    ggsave(filepath, plot, ...)
    cat(sprintf("\n![%s](%s)", filepath, filepath), file = outfile, append = T)
    if (print) print(plot)
  }
  
  facets <- "init+sample+N~network+struct"
  
  # Parent set size -----
  cat("\n\n# PARENT SET SIZE\n", file = outfile, append = T)
  df <- res_to_df(res[, "npar"], filepaths, names)
  df <- tidyr:::pivot_longer(df, names(df)[!names(df) %in% names], values_to = "avg_npar")
  y  <- "avg_npar"
  plot <- plot_box_plot(df, x = "epf", y = y, fill = "regular", facets = facets) 
  plotter(plot, outfile = outfile, print = print, 
          filepath = here::here(paste0(outdir, y, "_", tag, ".png")))
  
  # Plot rates ----
  df <- res_to_df(res[, "rates"], filepaths, names)
  for (y in c("avgppv", "TPR", "FPR")) {
    cat("\n\n#", ifelse(y == "avgppv", "AVERAGE PRECISION", y), "\n",
        file = outfile, append = T)
    plot <- plot_box_plot(df, x = "epf", y = y, fill = "regular", facets = facets) + ylim(0, 1) 
    plotter(plot, outfile = outfile, print = print, 
            filepath = here::here(paste0(outdir, y, "_", tag, ".png")))
  }
  
  # Plot prec-recall curves ----
  cat("\n\n# PRECISION RECALL \n", file = outfile, append = T)
  indx <- grepl("_epf1_|_epf10_|_epf100_", filepaths) 
  df <- res_to_df(res[indx, "pr"], filepaths[indx], names)
  
  for (x in unique(df$init)) {
    for (N in unique(df$N)) {
      facets <- "init+epf+sample~network+N+struct+regular"
      indx <- df$init == x & df$N == N
      plot <- plot_prec_recall(df[indx, ], facets = facets)
      plotter(plot, outfile = outfile, print = print, 
              filepath = here::here(outdir, paste("pr", tag, x, N, ".png", sep = "_")))
    }
  }
  

  
  # Plot running times 
  cat("\n# RUNTIMES \n", file = outfile, append = T)
  df <- res_to_df(res[, "toc"], filepaths, names)
  df_agg <- aggregate(df$tmp, df[, c("network", "init", "struct", "sample", "epf", "regular", "N", "name")], mean)
  
  plot <- ggplot(df_agg, aes(N, y = x, fill = name)) +
    facet_grid(init+epf+sample~network+struct+regular) +
    geom_col()
  plotter(plot, outfile = outfile, print = print, 
          filepath = here::here(outdir, paste("avgruntimes_", tag, ".png", sep = "")))
}


# run ----
# count number of runs for each simulation setting
filepaths <- list.files(indir, full.names = T)
df <- res_to_df(as.list(filepaths), filepaths, names)
head(df)
tidyr::pivot_wider(aggregate(rep(1, nrow(df)), df[, names[!names == "r"]], sum), 
                   names_from = "struct", values_from = "x") %>%  print(n = 100)


# set up cluster for foreach loop 
cl <- makeCluster(5, "SOCK", outfile = "")
export <- as.vector(lsf.str(envir = .GlobalEnv))
clusterExport(cl, export)
registerDoSNOW(cl)


# run plotting routine ----
res <- list()
for (bnname in bnnames) {
  filepaths <- list.files(indir, bnname, full.names = T)
  if (length(filepaths) == 0) next 
  
  tag <- paste0(simname, "_", bnname)
  outfile <- paste0(outdir, "results_", tag, ".md") 
  file.remove(outfile)
  res <-  eval_files(filepaths, bnname, prmethod = "noise", verbose = T)
  plots_to_file(res, filepaths, names, tag, outdir = outdir, outfile = outfile)
}
stopCluster(cl)




