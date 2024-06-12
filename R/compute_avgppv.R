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