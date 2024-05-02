

#' Title
#'
#' @param label 
#' @param nlev 
#'
#' @return
#' @export
#'
#' @examples
#' nlev  <- rep(2, 3)
#' label <- vector("list", length(nlev))
#' label[[1]] <- rbind(c(NA, 0), c(0, 1))
#' part <- label_to_partition(label, nlev)
#' conf <- expand.grid(lapply(nlev-1, seq.int, from = 0))
#' cbind(conf, part = part)
labels_to_partition <- function(labels, nlev) {
  
  n      <- length(nlev)
  seqn   <- seq_len(n)
  cump   <- cumprod(nlev)
  part   <- seq_len(cump[n])      # enumerate rows in full CPT

  stride <- c(1, cump[-n])
  lev <- lapply(nlev-1, seq.int, from = 0)
  
  for (i in seq_along(labels)[lengths(labels) > 0]) {
    for (r in 1:NROW(labels[[i]])) {
      context <- label[[i]][r, ]
      isna <-  is.na(context)
      rows0 <- c(context[!isna]%*%stride[-i][!isna])
      
      vars  <- c(i, seqn[-i][isna])
      rows  <- rows0 + bida:::expand_grid_fast(lev[vars])%*%stride[vars]
      
      part[rows+1] <- part[rows0+1]
    }
  }
  
  unname(split(seq_along(part)-1, part))
}
