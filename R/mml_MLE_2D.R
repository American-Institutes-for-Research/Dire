# updates the i,j element (and j,i element) of Sigma
# using a simple one-dimensional optimization of the likelihood
mml_MLE_2D <- function(values, i, j, verbose=0, tol=1e-3) {
  lnlf <- MLE_lnl_factory(xbi=values$Xb[,i], xbj=values$Xb[,j],
                          rr1i=values$rr1_list[[i]], rr1j=values$rr1_list[[j]],
                          w = values$w,
                          si2=values$Sigma[i,i], sj2=values$Sigma[j,j],
                          nodes=values$nodes, verbose=verbose)
  opt <- stats::optimize(f=lnlf, interval=c(-3,3), maximum=TRUE, tol=tol)
  values$Sigma[i,j] <- values$Sigma[j,i] <- tanh(opt$maximum) * sqrt(values$Sigma[i,i] * values$Sigma[j,j])
  return(values)
}

MLE_lnl_factory <- function(xbi, xbj, rr1i, rr1j, w, si2, sj2, nodes, verbose=-1) {
  # subset Xb
  # just cases in both
  ijss <- !is.na(xbi) & !is.na(xbj)
  # subset xs, rr1, weights to valid data
  xbi <- xbi[ijss]
  xbj <- xbj[ijss]
  w <-w[ijss]
  rr1i <- rr1i[ , ijss]
  rr1j <- rr1j[ , ijss]
  # this is the matrix of nodes that conforms to rr1
  Xb0 <- t(matrix(nodes, nrow=length(nodes), ncol=length(xbi)))
  Q <- length(nodes)
  # calculate q - Xb_i 
  nodes.minus.XBi <- t(Xb0 - xbi)
  nodes.minus.XBj <- t(Xb0 - xbj)
  Sigma <- matrix(0, nrow=2, ncol=2)
  Sigma[1,1] <- si2
  Sigma[2,2] <- sj2
  dTheta <- nodes[2] - nodes[1]
  function(rho, vec=FALSE) {
    Sigma[1,2] <- Sigma[2,1] <- tanh(rho) * sqrt(si2 * sj2)
    # temporary vectors to calculate integrals. One element per person in the sample
    pr_sum <- rep(0, length(xbi))
    Sigma_inv <- solve(Sigma)
    det_sigma_12 <- 1/(2*pi) * 1/sqrt(Sigma[1,1]*Sigma[2,2] - Sigma[1,2] * Sigma[2,1])
    # for all nodes in the i direction
    for(i in 1:Q) {
      resid_i <- nodes[i] - xbi
      for(j in 1:Q) { # this is 5.2 generalized to covariance
        resid_j <- nodes[j] - xbj
        mvn0 <- det_sigma_12 * exp(dmvn_log_exp(resid_i, resid_j, Sigma_inv[1,1], Sigma_inv[2,2], Sigma_inv[1,2])) 
        # this is pr(theta|responses) * phi term
        pr_part <- rr1i[i,] * rr1j[j,] * mvn0
        # for denominator
        pr_sum <- pr_sum + pr_part
      }
    }
    res <- sum(w*(log(pr_sum) + 2*log(dTheta)))
    if(verbose > 0) {
      cat(paste0("rho=", tanh(rho), "\tlnl=",res,"\n"))
    }
    if(vec) {
      return(pr_sum * dTheta^2)
    }
    return(res)
  }
}

# helper that finds the log of the multivariate normal distribution exponential part.
dmvn_log_exp <- function(x,y,Sinv11, Sinv22, Sinv12) {
  -0.5* (x*x*Sinv11 + y*y*Sinv22 + 2*x*y*Sinv12)
}
