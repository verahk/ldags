






optimize_local_csi <- function(counts, nlev, ess, r, q) {
  q   <- nrow(counts)
  r   <- ncol(counts)
  
  # init
  S  <- uS <- letters[seq_len(q)] # seq_len(q)   # partition and list with regions
  nS <- rep(1, q)
  scores <- 
    
  n <- length(nlev)
  labels <- vector("list", n)
  
  conf <- bida:::expand_grid_fast(lapply(nlev-1, seq.int, from = 0), nlev)
  
  q <- nrow(conf)
  rowindx <- seq_len(q)
  
  part <- list()
  part$value  <- part$unique <- seq_len(q)
  part$size   <- rep(1, q)
  part$counts <- counts
  part$score  <- famscore_bdeu_byrow(counts, ess, r, q, s = part_n)
  
  while (TRUE) {
    max_diff <- kappa
    for (i in seq_len(n)) {
      context <- conf[, -i, drop = FALSE]%*%stride[-i] +1
      
      # find the context that give the highest rise in the score
      # if added to current label
      for (x in context[conf[, i] == 0]) {
        if (x in labels[[i]]) next
        
        collapse <- match(part[context == x], part_unique)
        tmp_counts <- matrix(rowSums(part_counts[collapse, , drop = FALSE]), nrow = 1)
        tmp_score  <- famscore_bdeu_byrow(tmp_counts, ess, r, q, sum(part_size[collapse]))
        diff <- tmp_score - sum(part_score[collapse])
        if (diff > max_diff) {
          max_diff <- diff
          upd <- list(i = i, x = x, collapse, counts = tmp_counts, score = tmp_score)
        }
      }
      
      if (diff > max_diff) {
        
        # add context to label
        labels[[upd$i]] <- c(labels[[upd$i]], upd$x)
        
        # update counts, score and number of cases in each region
        old_regions <- part$unique[upd$collapse]
        keep <- old_regions[1]
        drop <- old_regions[2]
        
        part$score  <- replace(part$score[-drop], keep, upd$score)
        part$size   <- replace(part$size[-drop],  keep, upd$size)
        part$score  <- replace(part$score[-drop], keep, upd$score)
        
        part$counts[keep, ] <- upd$counts
        part$counts <- part$counts[-drop, ]
        
        # update partition
        indx <- part$value %in% old_regions
        part$value[indx] <- keep
      }
    }
  }
}

optimize_local_csi <- function(counts, nlev, ess) {
  q   <- nrow(counts)
  r   <- ncol(counts)
  
  # init
  S  <- uS <- letters[seq_len(q)] # seq_len(q)   # partition and list with regions
  nS <- rep(1, q)
  scores <- 
  
  n <- length(nlev)
  labels <- vector("list", n)

  conf <- bida:::expand_grid_fast(lapply(nlev-1, seq.int, from = 0), nlev)
  
  q <- nrow(conf)
  rowindx <- seq_len(q)
  
  part <- list()
  part$value  <- part$unique <- seq_len(q)
  part$size   <- rep(1, q)
  part$counts <- counts
  part$score  <- famscore_bdeu_byrow(counts, ess, r, q, s = part_n)
  
  while (TRUE) {
    max_diff <- kappa
    for (i in seq_len(n)) {
      context <- conf[, -i, drop = FALSE]%*%stride[-i] +1
      
      # find the context that give the highest rise in the score
      # if added to current label
      for (x in context[conf[, i] == 0]) {
        if (x in labels[[i]]) next
        
        collapse <- match(part[context == x], part_unique)
        tmp_counts <- matrix(rowSums(part_counts[collapse, , drop = FALSE]), nrow = 1)
        tmp_score  <- famscore_bdeu_byrow(tmp_counts, ess, r, q, sum(part_size[collapse]))
        diff <- tmp_score - sum(part_score[collapse])
        if (diff > max_diff) {
          max_diff <- diff
          upd <- list(i = i, x = x, collapse, counts = tmp_counts, score = tmp_score)
        }
      }
      
      if (diff > max_diff) {
        
        # add context to label
        labels[[upd$i]] <- c(labels[[upd$i]], upd$x)

        # update counts, score and number of cases in each region
        old_regions <- part$unique[upd$collapse]
        keep <- old_regions[1]
        drop <- old_regions[2]

        part$score  <- replace(part$score[-drop], keep, upd$score)
        part$size   <- replace(part$size[-drop],  keep, upd$size)
        part$score  <- replace(part$score[-drop], keep, upd$score)
        
        part$counts[keep, ] <- upd$counts
        part$counts <- part$counts[-drop, ]
        
        # update partition
        indx <- part$value %in% old_regions
        part$value[indx] <- keep
      }
    }
  }
}
