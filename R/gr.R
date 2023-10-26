fn_reg_Hessian <- function(X_, i=NULL, wv=NULL, rr1, stuDat, nodes, inside=TRUE, fast=TRUE) {
  if(length(wv) > 1) {
    stop("Must specify NULL weight or exactly one.")
  }
  # any lets the NULL case pass
  if(any(! wv %in% colnames(stuDat))) {
    stop("weight wv not in X_ wv=", dQuote(wv), " column names of X_", pasteItems(dQuote(colnames(X_))))
  }
  K <- ncol(X_)
  # fix rr1 to just regard this data
  if(!is.null(i)) {
    rr1p <- rr1[,i,drop=FALSE]
  } else{
    rr1p <- rr1
    i <- 1:nrow(X_)
  }
  # fix weights
  if(!is.null(wv)) {
    w <- stuDat[i, wv]
  } else {
    w <- rep(1, nrow(X_))
  }
  dTheta <- nodes[2] - nodes[1]
  # for numerical stability, premultiply by rr1m, then remove it when calculating sum(w * log(cs))
  rr1m <- (1/(.Machine$double.eps)^0.5)/apply(rr1p, 2, max) 
  rr1q <- t(t(rr1p) * rr1m) 
  insd <- inside
  function(par) {
    # break up par into beta and residual components
    B <- par[1:K]
    if(insd) {
      s2 <- ifelse(par[(K+1)] < 1, exp(par[(K+1)] - 1), par[(K+1)]^2)
    } else {
      s2 <- par[(K+1)]^2
    }
    s2 <- max(s2, 1e-6)
    s <- sqrt(s2)
    # form prediction
    XB <- as.vector(X_ %*% B)
    # residual, per indivudal (outer sapply), per node (inner sapply)
    nodes.minus.XB <- t(t(matrix(nodes, nrow=length(nodes), ncol=length(XB))) - XB)
    # likelihood of normal distribution for residuals nodes.minus.XB
    rr2 <- rr1p * ((1/(sqrt(2 * pi * s2))) * exp(-((nodes.minus.XB)^2/(2 * s2))))
    gr_res <- matrix(0, ncol=length(par),nrow=nrow(X_))
    trr2mxb <- t(rr2 * nodes.minus.XB)
    denom <- colSums(rr2)
    for(i in 1:K) {
      gr_res[,i] <- -2 * rowSums( trr2mxb * X_[,i] / s2)/denom
    }
    gr_res[,K+1] <- -2 * colSums( rr2 * (-1/s + 1*(nodes.minus.XB)^2/s^3))/denom
    if(par[(K+1)] < 1) {
      gr_res[,K+1] <- gr_res[,K+1] * 0.5 * exp(0.5*(par[(K+1)]-1)) 
    }
    gr_res[is.na(gr_res)] <- 0
    hess <- matrix(0, nrow=length(par), ncol=length(par))
    gr_res <- gr_res
    for(i in 1:nrow(gr_res)) {
      hess <- hess + w[i] * gr_res[i,] %o% gr_res[i,]
    }
    return(list(hess=hess,gr=apply(gr_res, 2, sum),gri = gr_res))
  }
}
