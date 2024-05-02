

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
#' label[[1]] <- rbind(c(0, 0), c(0, 1))
#' part <- label_to_partition(label, nlev)
#' conf <- expand.grid(lapply(nlev-1, seq.int, from = 0))
#' cbind(conf, part = part)
label_to_partition <- function(label, nlev = rep(2, length(label))) {
  
  n      <- length(nlev)
  cump   <- cumprod(nlev)
  stride <- c(1, cump[-n])
  part   <- seq_len(cump[n])      # enumerate rows in full CPT
  
  for (i in seq_along(labels)[lengths(label) > 0]) {

    contexts <- label[[i]]%*%stride[-i] 
    rows <- rep(contexts, nlev[i]) + rep(seq.int(0, nlev[i]-1), each = length(contexts)) +1
    dim(rows) <- c(2, 2)
    
    for (r in 1:nrow(rows)) {
      
      # rows in full CPT corresp to context
      indx    <- rows[r, ]
      
      # regions in partition corresp to context
      regions <- part[indx]            # regions that are collapsed by x
      region  <- min(regions)          # region regions is collapsed into
      
      # collapse regions in partition
      part[part %in% regions] <- region
    }
  }
  
  return(part)
}
