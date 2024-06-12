

#' Title
#'
#' @param x 
#' @param y 
#' @param method 
#'
#' @return
#' @keywords internal
#' @examples
#'
#' x <- c(0, 1, 1)
#' y <- c(0, 0, 1)
#' compare_prec_recall(x, y, niter = 20)
#' 
#' x <- c(1, 1, 3, 3, 3)
#' y <- c(0, 1, 1, 1, 0)
#' compare_prec_recall(x, y, niter = 20)
#' 
#' x <- c(rep(0, 20), 1:80)
#' y <- c(runif(20) < .2, runif(80) < c(1:80)/100)
#' compare_prec_recall(x, y, niter = 20)
#' 
#' 
#' 


compute_prec_recall <- function(x, y, method) {
  switch(method, 
         "noise" = compute_prec_recall_noise(x, y),
         "dups"  = compute_prec_recall_dups(x, y),
         "interp" = compute_prec_recall_interp(x, y))
}

compute_prec_recall_noise <- function(x, y) {
  indx <- order(x+runif(x)/10**5, decreasing = TRUE)
  tp   <- cumsum(y[indx])
  list(df = data.frame(x = x[indx], TPR = tp/sum(y), PPV = tp/seq_along(x)),
       avgppv = mean((tp/seq_along(x))[y[indx]>0]))
}

compute_prec_recall_dups <- function(x, y) {
  
  indx <- order(x, decreasing = TRUE)
  tp   <- cumsum(y[indx])
  pp   <- seq_along(x)
  n1   <- tp[length(y)]
  
  dups <- duplicated(x[indx], fromLast = TRUE)
  #cbind(x = x[indx], dups, y[indx], tp, pp, ppv = tp/pp)
  
  PPV  <- tp[!dups]/pp[!dups]
  w    <- diff(c(0, tp[!dups]))
  list(df = data.frame(x = x[indx][!dups],
                       TPR = tp[!dups]/n1,
                       PPV = PPV),
       avgppv = 1/n1*sum(w*PPV))
}

compute_prec_recall_interp <- function(x, y) {
  
  indx <- order(x, decreasing = TRUE)
  tp   <- cumsum(y[indx])
  pp   <- seq_along(x)
  fp   <- pp-tp
  n1   <- tp[length(y)]
  
  dups <- duplicated(x[indx], fromLast = TRUE)
  cbind(x = x[indx], y = y[indx], dups, tp, fp = pp-tp, pp, ppv = tp/pp)
  
  # interpolate false positives 
  tmp <- approx(c(0, tp[!dups]), 
                c(0, fp[!dups]),
                xout = seq_len(n1),
                ties = min)
  # plot(tp[!dups], fp[!dups])
  # lines(tmp$x, tmp$y)
  
  TPR <- tmp$x/n1
  PPV <- tmp$x/(tmp$x+tmp$y)
  
  list(df = data.frame(TPR = TPR, PPV = PPV),
       avgppv = 1/n1*sum(PPV))
}

compare_prec_recall <- function(x, y, niter = 10) {
  avgppvs <- numeric(niter)
  
  res <- compute_prec_recall(x, y, method = "dups")
  plot(res$df$TPR, res$df$PPV, 
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "recall", ylab = "precision")
  avgppvs[1] <- res$avgppv
  
  res <- compute_prec_recall(x, y, method = "interp")
  points(res$df$TPR, res$df$PPV, col = "blue")
  lines(res$df$TPR, res$df$PPV, col = "blue")
  avgppvs[2] <- res$avgppv
  
  
  for (i in seq_len(niter)) {
    res <- compute_prec_recall(x, y, method = "noise")
    lines(res$df$TPR, res$df$PPV, col = "red")
    avgppvs[i+2] <- res$avgppv
    cat("\navgppv, noise, ", i, ":", res$avgppv)
    
  }
  
  cat("\navgppv, dups:", avgppvs[1])
  cat("\navgppv, interp:", avgppvs[2])
  cat("\navgppv, noise, mean:",   mean(avgppvs[-c(1:2)]))
}

