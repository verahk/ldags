

unique_integer_fast <- function(x, nbins = as.integer(max(1L, x, na.rm = TRUE)), seqn = seq_len(nbins)) {
  seqn[tabulate(x, nbins)>0]
}

unique_integer_fast <- function(x, nbins = as.integer(max(1L, x, na.rm = TRUE)), seqn = seq_len(nbins)) {
  seqn[tabulate(x, nbins)>0]
}
if (FALSE) {
  dims <- list(n = c(2, 3, 4, 5), 
               nbins = c(5, 10, 20, 50, 100),
               method = c("unique", "fast", "fast_w_args"))
  
  res <- array(NA, lengths(dims), dims)
  for (n in dims$n) {
    for (nbins in dims$nbins) {
      x <- sample(nbins, len, T)
      seqn <- seq_len(nbins)
      
      stopifnot(all(sort(unique(x)) == unique_integer_fast(x)))
      stopifnot(all(sort(unique(x)) == unique_integer_fast(x, nbins, seqn)))
      
      tmp <- microbenchmark::microbenchmark(unique(x),
                                            unique_integer_fast(x),
                                            unique_integer_fast(x, nbins, seqn))
      
      res[paste(n), paste(nbins), ] <- summary(tmp)[, "mean"] 
    }
  }
  
  library(ggplot2)
  df <- reshape2::melt(res)
  ggplot(df, aes(n, value, color = method)) +
    facet_grid(.~nbins) +
    geom_line()
}

