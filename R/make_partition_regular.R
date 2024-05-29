


#' Make a partition regular
#' 
#' If the partition implies a conditional independence, return a finer partitioning
#' that does not. Specifically, if all outcomes of one or more variables are included 
#' in all parts of the partition, divide each part by the outcomes of these variables.
#'
#' @param P 
#' @param nlev 
#'
#' @return
#' @export
#'
#' @examples
#' nlev <- c(2, 2)
#' P <- list(c(0, 2), c(1, 3))
#' make_regular_partition(P, nlev)
make_partition_regular <- function(P, nlev) {

  if (all(lengths(P) >= min(nlev))) {
    stride <- c(1, cumprod(nlev[-length(nlev)]))

    for (i in seq_along(nlev)) {
      j <- 1
      while(j <= length(P)) {
        vals <- unique((P[[j]]%/%stride[i])%%nlev[i])
        if (length(vals) < nlev[i]) {
          break
        } else if (j == length(P)) {
          PP <- lapply(P, function(p) split(p, (p%/%stride[i])%%nlev[i]))
          P <- unlist(PP, recursive = F, use.names = FALSE)
        }
        j <- j+1
      }
    }
  }

  return(P)
}
