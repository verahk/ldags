
unlist_partition <- function(x, lens = lengths(x)) {
  rep.int(seq_along(x), lens)[order(unlist(x))]
}

observed_partition <- function(x, data, nlev) {
  stride <- c(1, cumprod(nlev[-length(nlev)]))
  unlist_partition(x)[data%*%stride+1]
}
#table(data[, 1:2]%*%stride+1, observed_partition(res$partition, data[, 1:2], c(2, 2)))