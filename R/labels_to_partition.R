

#' Title
#'
#' @param label (list)
#'  list with labels on edge from each parent to their common child. 
#'  Labels specified as a matrix with each row the joint outcomes of co-parents, 
#'  or as a vector that enumerates these joint outcomes (the rowindex of the CPT
#'  table starting from 0). See example.
#' @param levels (list)
#'  list with levels of each parent
#' @param nlev (integer vector)
#' @param type (character)
#'  how to read labels. Either `"outcome_vectors"` (default) or `"enumerated"`.
#'
#' @return
#' @export
#'
#' @examples
#' levels <- list(0:1, 0:2, 0:3)
#' label <- vector("list", length(levels))
#' label[[1]] <- rbind(c(NA, 0), c(0, 1))
#' part1 <- labels_to_partition(label, levels, type = "outcome_vectors")
#' 
#' label[[1]] <- c(0, 2, 4, 6)
#' label[[2]] <- c(0)
#' part2 <- labels_to_partition(label, levels, type = "enumerated")
#' 
#' cbind(expand.grid(levels),  get_partitions(part1), get_partitions(part2))
#'
labels_to_partition <- function(labels, levels, nlev = lengths(levels), type = "outcome_vectors") {
  
  n      <- length(nlev)
  cump   <- cumprod(nlev)
  stride <- c(1, cump[-n])
  parts   <- seq_len(cump[n])      # enumerate rows in full CPT
    
  if (type == "outcome_vectors") {
    for (i in seq_along(labels)[lengths(labels) > 0]) {
      # from set of joint outcomes to row in CPT
      lab <- labels[[i]]
      for (r in 1:nrow(lab))  {
        indx <- is.na(lab[r, ])
        if (any(indx)) {
          tmp <- bida:::expand_grid_fast(levels[-i][indx])
          mat <- matrix(NA, nrow = nrow(tmp), ncol = n-1)
          mat[, indx] <- tmp
          mat[, !indx] <- rep(lab[r, !indx], each = nrow(tmp))
          tmp <-mat%*%stride[-i]
        } else {
          tmp <- lab[r, ]%*%stride[-i]
        }
        rows <- outer(c(tmp), levels[[i]]*stride[i] +1, "+")
        regions_to_be_collapsed <- parts[rows]
        parts[parts %in% regions_to_be_collapsed] <- min(regions_to_be_collapsed)
      }
    }

  } else if (type == "enumerated") {
    for (i in seq_along(labels)[lengths(labels) > 0]) {
      rows <- outer(labels[[i]], levels[[i]]*stride[i] +1, "+")
      for (r in 1:nrow(rows)) {
        regions_to_be_collapsed <- parts[rows[r, ]]
        parts[parts %in% regions_to_be_collapsed] <- min(regions_to_be_collapsed)
      }
    }
  }
   
  unname(split(seq_along(parts)-1, parts))
}


