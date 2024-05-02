


new_bida_pair_cat <- function(x, y, params, support, zerosupp, nlev, ess) {
  structure(
    list(x = x,
         y = y, 
         params = params, 
         support = support, 
         zerosupp = zerosupp,
         nlev = nlev,
         ess = ess),
    class = c("bida_pair", "bida_pair_cat"))
}


bida_pairs_local_csi <- function(data, nlev, x, yvars, sets, support, P, ess, kappa) {
 
  counts <- matrix(list(), nrow = length(support), ncol = length(yvars))
  seqy <- seq_along(yvars)
  
  for (r in seq_along(support)) {
    z <- sets[r, ]
    z <- z[!is.na(z)]
    
    if (length(z) == 0) {
      
      # map observations to joint outcomes
      xy    <- data[, x] + nlev[x]*data[, yvars, drop = FALSE] + 1
      nbins <- nlev[x]*nlev[yvars]
      for (yy in seqy) {
        nyx <- matrix(tabulate(xy[, yy], nbins = nbins[yy]), 
                      ncol = nlev[x], byrow = T)
        counts[[r, yy]] <- nyx + ess/length(nyx)
      }
    } else if (length(z) == 1) {
      
      zxy   <- data[, z] + nlev[z]*data[, x] + nlev[x]*nlev[z]*data[, yvars, drop = FALSE] + 1
      nbins <- nlev[z]*nlev[x]*nlev[yvars]
      for (yy in seqy) {
        y <- yvars[yy]
        if (!z == y) {
          nzxy <- tabulate(zxy[, yy], nbins = nbins[yy])
          dim(nzxy) <- nlev[c(z, x, y)]
          counts[[r, yy]] <- aperm(nzxy, c(3, 2, 1)) + ess/length(nzxy)
        }
      }
    } else {
      
      if (is.list(P)) {
        part <- P[[r]]
        npart <- length(part)
        p <- observed_partition(part, data[, z], nlev[z]) 
      } else if (P == "opt") {
        part   <- optimize_CSI_structure(data, nlev, x, z, ess, kappa)$partition
        npart <- length(part)
        p <- observed_partition(part, data[, z], nlev[z])
      } else if (P == "none") {
          npart <- prod(nlev[z])
          part <- as.list(seq_len(npart)-1)
          p <- c(data[, z]%*%c(1, cumprod(nlev[z[-length(z)]]))) +1
      }
      
      
      #paconf <- data[, z]%*%c(1, cumprod(nlev[z[-length(z)]]))
      #table(paconf, p)
     
      # compute marginal counts
      pxy <- p + npart*data[, x] + nlev[x]*npart*data[, yvars, drop = FALSE]
      nbins <- npart*nlev[x]*nlev[yvars]
      for (yy in seqy) {
        y <- yvars[yy]
        if (!any(z == y)) {
          npxy <- tabulate(pxy[, yy], nbins = nbins[yy])
          dim(npxy) <- c(npart, nlev[x], nlev[y])
          counts[[r, yy]] <- aperm(npxy, c(3, 2, 1)) + ess*lengths(part)/prod(nlev[c(x, y, z)])
          
          # check
          # df <- as.data.frame(data) 
          # tab <- table(df[c(z, x, y)])
          # dim(tab) <- c(prod(nlev[z]), nlev[x], nlev[y])
          # for (i in seq_along(part)) {
          #   stopifnot(all(colSums(tab[part[[i]]+1,,, drop = FALSE]) == npxy[i,,]))
          # } 
        }
      }
    }
  }
  
  # return bida pairs
  out <- list()
  zero_effects <- lengths(counts) == 0  # 
  for (yy in seqy) {
    
    # param for zero-effects
    indx <- zero_effects[, yy]
    if (any(indx)) {
      # add counts for non-zero effects
      y <- yvars[yy]
      ny <- tabulate(data[, y]+1, nlev[y]) + ess/nlev[y]
      zerosupp <- sum(support[indx])
 
      out[[yy]] <- new_bida_pair_cat(x, y, c(list(ny), counts[!indx, yy]),
                                     support = c(zerosupp, support[!indx]),
                                     zerosupp = zerosupp,
                                     nlev = nlev[c(x, y)],
                                     ess = 0)
    } else {
      out[[yy]] <- new_bida_pair_cat(x, yvars[yy], counts[, yy], support, zerosupp = 0, nlev = nlev[c(x, y)], ess = 0)
    }
  }
  
  return(out)
}

# test ----
if (FALSE) {
  bn <- readRDS("./data/asia.rds")
  n <- length(bn)
  nlev <- sapply(bn, function(x) dim(x$prob)[1])
  data <- bida:::sample_data_from_bn(bn, 100)
  
  x <- 6
  sets  <- rbind(c(NA, rep(NA, 2)), 
                 c(1,  rep(NA, 2)),
                 c(2, 4, rep(NA, 1)),
                 c(1, 2, 4))
  support <- rep(1, nrow(sets))/nrow(sets)
 
  pairs <- bida_pairs_local_csi(data, nlev, x, yvars = 1:5, sets, support, "opt", 1, 1)
  pairs[[1]]
  pairs[[2]]
  pairs[[3]]
  
  bida:::posterior_mean.backdoor_params_cat(pairs[[1]]$params[[1]], pairs[[1]]$ess, pairs[[1]]$nlev[1])
  bida:::posterior_mean(pairs[[1]])
  means <- lapply(pairs, bida:::posterior_mean)
}


