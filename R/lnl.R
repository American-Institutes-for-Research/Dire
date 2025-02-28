### This is the marginal likelihood function to be maximized
# it is weighted if weight var is not null
#par - the score vector
#getIndLnL - logical, should the likelihood be calculated only for a subset
#X_ - the matrix of X values, defaults to all X , used for calculating group likelihoods
#i - list of indexes for individual likelihood (used to subset weights and r22)
#repweightVar  - the replicate weight for replicate variance estimation
#inside - TRUE if inside mix, otherwise linear and inside=FALSE
fn.regression <- function(X_, i=NULL, wv=NULL, rr1, stuDat, nodes, inside=TRUE, fast=TRUE, verbose=FALSE) {
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
  rm(stuDat) # this can save a lot of memory
  dTheta <- nodes[2] - nodes[1]
  # for numerical stability, premultiply by rr1m, then remove it when calculating sum(w * log(cs))
  rr1m <- (1/(.Machine$double.eps)^0.5)/apply(rr1p, 2, max)
  rr1q <- t(t(rr1p) * rr1m) 
  insd <- inside
  Xnodes <- t(matrix(nodes, nrow=length(nodes), ncol=nrow(X_)))

  prev_gr <- rep(Inf, ncol(X_))
  max_x_gr <- ncol(X_) * 1e-5
  function(par, returnPosterior=FALSE, gr=FALSE, hess=FALSE) {
    # break up par into beta and residual components
    B <- par[1:K]
    if(insd) {
      s2 <- ifelse(par[(K+1)] < 1, exp(par[(K+1)] - 1), par[(K+1)]^2)
    } else {
      s2 <- par[(K+1)]^2
    }
    s2 <- max(s2, 1e-6)
    s <- sqrt(s2)
    if(hess) {
      stop("Calculating the Hessian directly has been deprecated. Please use the list of functions which contains an alternative Hessian.")
    }
    if(gr) {
      XB <- as.vector(X_ %*% B)
      # residual, per indivudal (outer sapply), per node (inner sapply)
      nodes.minus.XB <- t(Xnodes - XB)
      # likelihood of normal distribution for residuals nodes.minus.XB
      rr2 <- rr1p * ((1/(sqrt(2 * pi * s2))) * exp(-((nodes.minus.XB)^2/(2 * s2))))
      gr_res <- rep(0, length(par))
      trr2mxb <- t(rr2 * nodes.minus.XB)
      denom <- colSums(rr2)
      rs <- w * rowSums( trr2mxb / s2)/denom
      for(i in 1:K){
        gr_res[i] <- -2 * sum( X_[,i] * rs)
        # works like: gr_res[i] <- -2 * sum(w * rowSums( trr2mxb * X_[,i] / s2)/denom)
      }
      gr_res[K+1] <- -2* sum(w * colSums( rr2 * (-1/s + 1*(nodes.minus.XB)^2/s^3))/denom)
      if(par[(K+1)] < 1) {
        gr_res[K+1] <- gr_res[K+1] * 0.5 * exp(0.5*(par[(K+1)]-1)) 
      }
      gr_res[is.na(gr_res)] <- 0
      prev_gr <<- gr_res
      return(gr_res) 
    }
    # form prediction
    XB <- X_ %*% B 
    # residual, per indivudal (outer sapply), per node (inner sapply)
    nodes.minus.XB <- t(t(matrix(nodes, nrow=length(nodes), ncol=nrow(XB))) - as.vector(XB)) 
    # likelihood of normal distribution for residuals nodes.minus.XB
    rr2 <- rr1q * ((1/(sqrt(2 * pi * s2))) * exp(-((nodes.minus.XB)^2/(2 * s2))))
    # aggregate likelihood, weight, multiply by -2 to make it a deviance
    cs <- pmax(.Machine$double.eps, colSums(rr2))
    if(returnPosterior) {
      sd <- mu <- rep(0, ncol(rr1q)) 
      for(i in 1:ncol(rr1q)) {
        fi <- rr2[,i]
        mu[i] <- sum(nodes*fi)/sum(fi)
        sd[i] <- sqrt(sum((nodes - mu[i])^2 * fi)/sum(fi))
      }
      return(data.frame(id=rownames(X_), mu=mu, sd=sd, stringsAsFactors=FALSE))
    }
    res <- -2*( sum(w * (log(cs) - log(rr1m)) ) + sum(w) * log(dTheta) )
    return( res )
  }
}

fn.regression_lnl <- function(X_, w, rr1, nodes) {
  K <- ncol(X_)
  rr1p <- rr1
  dTheta <- nodes[2] - nodes[1]
  # for numerical stability, premultiply by rr1m, then remove it when calculating sum(w * log(cs))
  rr1m <- (1/(.Machine$double.eps)^0.5)/apply(rr1p, 2, max)
  rr1q <- t(t(rr1p) * rr1m) 
  Xnodes <- t(matrix(nodes, nrow=length(nodes), ncol=nrow(X_)))
  prev_gr <- rep(Inf, ncol(X_))
  max_x_gr <- ncol(X_) * 1e-5
  Xb0 <- t(matrix(nodes, nrow=length(nodes), ncol=nrow(X_)))
  function(par, returnCS=FALSE) {
    # break up par into beta and residual components
    B <- par[1:K]
    s1 <- par[(K+1)]
    s2 <- s1^2
    s <- sqrt(s2)
    # form prediction
    XB <- X_ %*% B 
    # residual, per indivudal (outer sapply), per node (inner sapply)
    nodes.minus.XB <- t(Xb0 - as.vector(XB)) 
    # likelihood of normal distribution for residuals nodes.minus.XB
    rr2 <- rr1q * ((1/(sqrt(2 * pi * s2))) * exp(-((nodes.minus.XB)^2/(2 * s2))))
    # aggregate likelihood, weight, multiply by -2 to make it a deviance
    cs <- pmax(.Machine$double.eps, colSums(rr2))
    if(returnCS) {
      return(log(cs))
    }
    res <- -2*( sum(w * (log(cs) - log(rr1m)) ) + sum(w) * log(dTheta) )
    return( res )
  }
}

fn.regression_gr_lnl <- function(X_, w, rr1, nodes) {
  K <- ncol(X_)
  dTheta <- nodes[2] - nodes[1]
  # for numerical stability, premultiply by rr1m, then remove it when calculating sum(w * log(cs))
  Xnodes <- t(matrix(nodes, nrow=length(nodes), ncol=nrow(X_)))
  prev_gr <- rep(Inf, ncol(X_))
  max_x_gr <- ncol(X_) * 1e-5
  Xb0 <- t(matrix(nodes, nrow=length(nodes), ncol=nrow(X_)))
  function(par, sqrtW=FALSE) {
    # break up par into beta and residual components
    B <- par[1:K]
    s1 <- par[(K+1)]
    s2 <- s1^2
    s <- sqrt(s2)
    # form prediction
    XB <- X_ %*% B 
    # residual, per indivudal (outer sapply), per node (inner sapply)
    nodes.minus.XB <- t(Xb0 - as.vector(XB)) 
    # terms that are not a function of the node nor b cancel out
    rr20 <- rr1 * exp(-((nodes.minus.XB)^2/(2 * s2)))
    rr2p <- rr1 * exp(-((nodes.minus.XB)^2/(2 * s2))) * ((nodes.minus.XB)/(s2))
    # aggregate likelihood, weight, multiply by -2 to make it a deviance
    cs0 <- pmax(.Machine$double.eps, colSums(rr20))
    csp <- colSums(rr2p)
    if(sqrtW) {
      res <- sqrt(w)*csp/cs0
    } else {
      res <- w*csp/cs0
    }
    return( res )
  }
}

fn.regression_hess_lnl <- function(X_, w, rr1, nodes) {
  K <- ncol(X_)
  dTheta <- nodes[2] - nodes[1]
  # for numerical stability, premultiply by rr1m, then remove it when calculating sum(w * log(cs))
  Xnodes <- t(matrix(nodes, nrow=length(nodes), ncol=nrow(X_)))
  prev_gr <- rep(Inf, ncol(X_))
  max_x_gr <- ncol(X_) * 1e-5
  Xb0 <- t(matrix(nodes, nrow=length(nodes), ncol=nrow(X_)))
  function(par) {
    # break up par into beta and residual components
    B <- par[1:K]
    s1 <- par[(K+1)]
    s2 <- s1^2
    s <- sqrt(s2)
    # form prediction
    XB <- as.vector(X_ %*% B) 
    # residual, per indivudal (outer sapply), per node (inner sapply)
    nodes.minus.XB <- t(Xnodes - XB)
    # likelihood of normal distribution for residuals nodes.minus.XB
    rr2 <- rr1 * ((1/(sqrt(2 * pi * s2))) * exp(-((nodes.minus.XB)^2/(2 * s2))))
    trr2mxb <- t(rr2 * nodes.minus.XB)
    colSums_rr2 <- colSums(rr2)
    # Hess loop in rcpp
    H <- slowHess_no_s(K, rr2, trr2mxb, X_, nodes.minus.XB, w, s2)
    return(H)
  }
}

slowHess_no_s <- function(K, rr2, trr2mxb, X_, nodes.minus.XB, w, s2) {
  H <- matrix(0, nrow=K, ncol=K)
  denom <- colSums(rr2)
  denom2 <- denom^2
  denomPrimeList <- list()
  for(i in 1:K) {
    fi <- t( trr2mxb * X_[,i]/s2)
    num <- colSums(fi)
    for(j in i:K) {
      fj <- t( trr2mxb * X_[,j]/s2)
      numPrime <- colSums(t(t(rr2) * X_[,i] * X_[,j] / s2) - t( t(fj * nodes.minus.XB) * X_[,i])/s2)
      denomPrime <- colSums(fj)
      H[j,i] <- H[i,j] <- 2*sum(w*(numPrime * denom + num * denomPrime) / denom2)
    }
  }
  return(H)
}
