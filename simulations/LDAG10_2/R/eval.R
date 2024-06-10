eval_files <- function(filepaths, bnname) {
  
  indx <- grepl(bnname, filepaths, T)
  if (!all(indx)) stop("All files do not match bnname.")
  
  bn <- readRDS(here::here("./data/", paste0(bnname, ".rds")))
  n  <-  length(bn)
  dindx <- diag(n) == 1
  dag  <- bnlearn::amat(bn)[!dindx]
  dmat <- bida:::descendants(bn)[!dindx]
  
  export <- c("eval_MCMCchain", "compute_prec_recall", "compute_avgppv")
  res <- foreach(f = filepaths, 
                 .combine  = "rbind",
                 .packages = "Matrix",
                 .export = export) %dopar% eval_MCMCchain(f, seq_len(200), dag, dmat, n)
}


eval_MCMCchain <- function(f, burninsamples, dag, dmat, n, verbose = T) {
  # if (verbose) cat(f)
  MCMCchain  <- readRDS(f)
  if ("smpl" %in% names(MCMCchain))  MCMCchain <- MCMCchain$smpl
 
  dindx <- diag(n) == 1
  dags <- lapply(MCMCchain$traceadd$incidence[-burninsamples], as.matrix)
  
  # list unique DAGs
  u    <- unique(dags)
  support  <- bida:::rowsum_fast(rep(1/length(dags), length(dags)), dags, u)
  dags <- u 
  dmats <- lapply(dags, bida:::descendants)
  
  edgep <- Reduce("+", Map("*", dags, support))[!dindx]
  ancp  <- Reduce("+", Map("*", dmats, support))[!dindx]
  
  res <- matrix(nrow = 2, ncol = 3)
  colnames(res) <- c("TPR", "FPR", "avgppv")
  rownames(res) <- c("edgep", "ancp")
  
  # edge probs
  n1 <- sum(dag)

  res[1, 1] <- sum(edgep[dag == 1])/n1
  res[1, 2] <- sum(edgep[dag == 0])/(n**2-n-n1)
  res[1, 3] <- compute_avgppv(c(edgep), c(dag))
  
  # arp
  pr <- compute_prec_recall_noise(ancp, dmat)
  
  n1 <- sum(dmat)
  probs <- Reduce("+", Map("*", dmats, support))[!dindx]
  res[2, 1] <- sum(ancp[dmat == 1])/n1
  res[2, 2] <- sum(ancp[dmat == 0])/(n**2-n-n1)
  #res[2, 3] <- compute_avgppv(ancp, dmat)
  res[2, 3] <- pr$avgppv 
  
  list(edgep = edgep, 
       rates = res, 
       pr = pr$df,
       toc = as.matrix(attr(MCMCchain, "toc")))
  
}

compute_prec_recall_noise <- function(x, y) {
  indx <- order(x+runif(x)/10**5, decreasing = TRUE)
  tp   <- cumsum(y[indx])
  list(df = data.frame(x = x[indx], TPR = tp/sum(y), PPV = tp/seq_along(x)),
       avgppv = mean((tp/seq_along(x))[y[indx]>0]))
}

compute_prec_recall <- function(x, y, expand = F) {
  
  indx <- order(x, decreasing = TRUE)
  tp   <- cumsum(y[indx])
  pp   <- seq_along(x)
  n1   <- tp[length(y)]
  
  dups <- duplicated(x[indx], fromLast = TRUE)
  PPV  <- tp[!dups]/pp[!dups]
  w    <- diff(c(0, tp[!dups]))

  if (expand){
    rate <- PPV[length(PPV)]
    list(df = data.frame(x = c(Inf, x[indx][!dups], -Inf),
                         TPR = c(0, tp[!dups]/n1, 1),
                         PPV = c(1, PPV, rate)),
         avgppv = 1/n1*sum(w*PPV))
  } else {
    list(df = data.frame(x = x[indx][!dups],
                         TPR = tp[!dups]/n1,
                         PPV = PPV),
         avgppv = 1/n1*sum(w*PPV))
  }

}

if (FALSE) {
  x <- c(rep(1, 10), 0, 0, 0, .5)
  y <- c(rep(1, 3), rep(0, 7), c(0, 1, 1, 0))
  tmp <- compute_prec_recall(x, y)
  plot(tmp$df$TPR, tmp$df$PPV, xlim = c(0, 1), ylim = c(0, 1))
  lines(tmp$df$TPR, tmp$df$PPV)
  cat("\navgppv", tmp$avgppv)
  
  res <- numeric(100)
  for (i in seq_along(res)) {
    tmp <- compute_prec_recall_noise(x, y)
    lines(tmp$df$TPR, tmp$df$PPV, col = "red")
    res[i] <- tmp$avgppv
  }
  cat("\navgppv", mean(res))
}

compute_avgppv <- function(x, y) {
  
  indx <- order(x, decreasing = TRUE)
  tp   <- cumsum(y[indx])
  pp   <- seq_along(x)
  
  dups <- duplicated(x[indx], fromLast = TRUE)
  #cbind(dups, x[indx], c(y[indx]), tp, pp)
  
  w <- diff(c(0, tp[!dups]))/tp[length(x)]
  #cbind(tp = tp[!dups], pp = pp[!dups], tpr = tp[!dups]/tp[length(x)], ppv = tp[!dups]/pp[!dups], w)
  
  sum(w*tp[!dups]/pp[!dups])
}

compute_avgppv_noise <- function(x, y) {
  indx <- order(x+runif(length(x))/1000, decreasing = TRUE)
  tp <- cumsum(y[indx])
  pp <- seq_along(x)
  mean((tp/pp)[y[indx] == 1])
}
res_to_df <- function(res, filepaths, names = c("network", "algo", "struct", "N", "epf", "r")) {
  
  tmp <- stringr::str_split(filepaths, ".+MCMCchains/|_|.rds", simplify = F)
  par <- data.frame(do.call(rbind, tmp)[, seq_along(names)+1])
  colnames(par) <- names
  par$N <- factor(par$N, c("N100", "N300", "N1000", "N3000", "N10000"))
  par$epf <- factor(par$epf, paste0("epf", c(1, 2, 4, 10**3, 10**4, 10**5)))
  
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


