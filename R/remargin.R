# formula 3.3
getM <- function(rr2, nodes, w) {
  # get E(theta), this is intended to be unweighted
  cs <- colSums(rr2*nodes)
  return(as.vector(cs))
}

# formula 3.6
getSoneD <- function(rr2, nodes.minus.XB, w) {
  res <- colSums( (nodes.minus.XB)^2 * rr2)
  return(sum(w*res)/sum(w))
}

updateGamma <- function(rr1, nodes, X, w, tol=10*sqrt(.Machine$double.eps), maxiter=10000, verbose=0) {
  # initalize B to all zeros, s2 to 1 (which is way too large)
  B <- rep(0, ncol(X))
  s2 <- 1
  # oldB is for finding the change in B
  oldB <- B
  dB <- B + 1
  # count iterations
  iter <- 1
  # find cases with any variance in rr1, can fit with all cases, but it just slows things down
  used <- ((apply(rr1, 2, sd)/pmax(.Machine$double.xmin , apply(rr1, 2, mean))) > .Machine$double.eps)
  if(verbose > 0 & sum(!used) > 0) {
    cat(paste0("", sum(!used), " cases not used (of ", length(used)," total cases; ",round(100*mean(used),1),"%) because of no variance.\n"))
  }
  X <- X[used, , drop=FALSE]
  w <- w[used]
  rr1 <- rr1[,used]
  # make W matrix
  W <- Matrix::Diagonal(x=w)
  Xb0 <- t(matrix(nodes, nrow=length(nodes), ncol=nrow(X)))
  while(max(abs(dB)/pmax(1, abs(B))) > tol & iter < maxiter) {
    XB <- as.vector(X %*% B)
    nodes.minus.XB <- t(Xb0 - XB) 
    rr2 <- rr1 * ((1/(sqrt(2 * pi * s2))) * exp(-((nodes.minus.XB)^2/(2 * s2))))
    w2 <- 1/colSums(rr2)
    rr2 <- t(t(rr2) * w2)
    M <- getM(rr2, nodes, w=w)
    s2 <- getSoneD(rr2=rr2, nodes.minus.XB=nodes.minus.XB, w=w)
    B2 <- qr.coef(qr(sqrt(W) %*% X), sqrt(W) %*% M)
    B <- B2
    dB <- (B2 - oldB)
    oldB <- B
    ff <- format(c(max(abs(dB)/pmax(1, abs(B))),tol), scientific=FALSE, digits=2, zero.print=TRUE)
    if(verbose > 0) {
      cat(paste0("iter ", iter, "\tgrad: ", ff[1], "\ttarget grad: ", ff[2], "\n"))
    }
    iter <- iter + 1
  }
  if(iter < maxiter) {
    convergence <- "converged"
  } else {
    convergence <- "over iteration max"
  }
  res <- list(par=c(unname(B),sqrt(s2)),
              iterations=iter-1,
              convergence=convergence)
  return(res)
}
