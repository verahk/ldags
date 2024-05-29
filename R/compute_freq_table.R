

compute_freq_table <- function(data, nlev, j, parentnodes, lookup = NULL) {
  
  npar <- length(parentnodes)
  r <- nlev[j]
  if (npar == 0) {
    matrix(tabulate(data[, j] +1, r), nrow = 1)
  } else  if (is.null(lookup)) {
    stride <- c(1, cumprod(nlev[parentnodes]))
    r <- nlev[j]
    q <- stride[length(stride)] 
    joint  <- data[, c(parentnodes, j), drop = FALSE]%*%stride + 1 # joint observed outcomes
    matrix(tabulate(joint, q*r), nrow = q, ncol = r)
  } else {
    parentnodes <- sort(parentnodes)
    id <- paste0(c(j, parentnodes), collapse = "|")
    counts <- lookup$counts[[id]]
    if (is.null(counts)) {
      counts <- compute_freq_table(data, nlev, j, parentnodes, lookup = NULL)
      lookup$counts[[id]] <- counts
    }
    return(counts)
  }
}

if (FALSE) {
  # compare search for counts from lookup vs computing directly
  data <- replicate(10, sample(0:2, 100, T))
  nlev <- rep(3, 10)
  tab <- rlang::new_environment(list(counts = list()))
  res <- list()
  for (n in 1:9) {
    par <- 1 + seq.int(n)
    # pre-compute
    compute_freq_table(data, nlev, j = 1, par, lookup = tab)
    tmp <- microbenchmark::microbenchmark(lookup = compute_freq_table(data, nlev, j = 1, par, lookup = tab), 
                                          nolookup = compute_freq_table(data, nlev, j = 1, par, lookup = NULL))
    res[[n]] <- summary(tmp)[, "mean"]
  }

}

compute_freq_table_from_cpt <- function(cpt, data) {
  n <- length(cpt$levels)
  j <- cpt$node
  r <- length(cpt$levels[[n]])
  if (n == 1) {
    matrix(tabulate(data[, j] +1, r), nrow = 1)
  } else if (is.null(cpt$partition)) {
    nlev <- lengths(cpt$levels)
    stride <- c(1, cumprod(nlev[-n]))
    q <- stride[n] 
    joint  <- data[, c(cpt$parentnodes, j), drop = FALSE]%*%stride  # joint observed outcomes
    matrix(tabulate(joint+1, q*r), nrow = q, ncol = r)
  } else {
    q <- length(cpt$partition)
    parts <- unlist_partition(cpt$partition) -1    
    joint <- data[, j] + r*parts                 # joint observed outcomes
    matrix(tabulate(joint+1, q*r), nrow = q, ncol = r)
  }
}