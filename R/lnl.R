### This is the marginal likelihood function to be maximized
# it is weighted if weight var is not null
#par - the score vector
#getIndLnL - logical, should the likelihood be calculated only for a subset
#X_ - the matrix of X values, defaults to all X , used for calculating group likelihoods
#i - list of indexes for individual likelihood (used to subset weights and r22)
#repweightVar  - the replicate weight for replicate variance estimation
fn.regression <- function(X_, i=NULL, wv=NULL, rr1, stuDat, nodes) {
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
  function(par, returnPosterior=FALSE, gr=FALSE, hess=FALSE) {
    if(hess) {
      # break up par into beta and residual components
      B <- par[1:K]
      s2 <- par[(K+1)]^2
      s <- sqrt(s2)
      # form prediction
      XB <- X_ %*% B
      # residual, per indivudal (outer sapply), per node (inner sapply)
      nodes.minus.XB <- t(t(matrix(nodes, nrow=length(nodes), ncol=nrow(XB))) - as.vector(XB))
      # likelihood of normal distribution for residuals nodes.minus.XB
      rr2 <- rr1p * ((1/(sqrt(2 * pi * s2))) * exp(-((nodes.minus.XB)^2/(2 * s2))))
      H <- matrix(0,nrow=K+1,ncol=K+1)
      for(i in 1:K) {
        for(j in i:K) {
          fi <- t(   t((rr2 * nodes.minus.XB)) * X_[,i]/s2)
          fj <- t(   t((rr2 * nodes.minus.XB)) * X_[,j]/s2)
          num <- apply(fi, 2, sum)
          numPrime <- apply(t(t(rr2) * X_[,i]*X_[,j]/s2) - t( t(fj * nodes.minus.XB) * X_[,i])/s2, 2, sum)
          denom <- colSums(rr2)
          denomPrime <- apply(fj, 2, sum)
          H[j,i] <- H[i,j] <- 2*sum(w*(numPrime * denom + num * denomPrime) / denom^2) # not quite
        }
      }
      for(i in 1:K) {
        gr0 <- -2 * sum(w*apply( t(   t((rr2 * nodes.minus.XB)) * X_[,i]/s2), 2, sum)/colSums(rr2))
        # break up par into beta and residual components
        s_ <- s + 1e-6
        # likelihood of normal distribution for residuals nodes.minus.XB
        rr2_ <- rr1p * ((1/(sqrt(2 * pi * s_^2))) * exp(-((nodes.minus.XB)^2/(2 * s_^2))))
        gr_ <- -2 * sum(w*apply( t(   t((rr2_ * nodes.minus.XB)) * X_[,i]/s_^2), 2, sum)/colSums(rr2_))
        H[K+1,i] <- H[i,K+1] <- (gr_ - gr0)/1e-6
      }
      gr0 <- -2*sum(w*apply( rr2 * (-1/s + 1*(nodes.minus.XB)^2/s^3), 2, sum)/colSums(rr2))
      # break up par into beta and residual components
      s_ <- s + 1e-6
      # likelihood of normal distribution for residuals nodes.minus.XB
      rr2_ <- rr1p * ((1/(sqrt(2 * pi * s_^2))) * exp(-((nodes.minus.XB)^2/(2 * s_^2))))
      gr_ <- -2*sum(w*apply( rr2_ * (-1/s_ + 1*(nodes.minus.XB)^2/s_^3), 2, sum)/colSums(rr2_))
      H[K+1,K+1] <- (gr_ - gr0)/1e-6
      return(H)
    }
    if(gr) {
      # break up par into beta and residual components
      B <- par[1:K]
      s2 <- par[(K+1)]^2
      s <- sqrt(s2)
      # form prediction
      XB <- X_ %*% B
      # residual, per indivudal (outer sapply), per node (inner sapply)
      nodes.minus.XB <- t(t(matrix(nodes, nrow=length(nodes), ncol=nrow(XB))) - as.vector(XB))
      # likelihood of normal distribution for residuals nodes.minus.XB
      rr2 <- rr1p * ((1/(sqrt(2 * pi * s2))) * exp(-((nodes.minus.XB)^2/(2 * s2))))
      gr_res <- rep(0,length(par))
      for(i in 1:K) {
        gr_res[i] <- -2 * sum(w*apply( t(   t((rr2 * nodes.minus.XB)) * X_[,i]/s2), 2, sum)/colSums(rr2))
      }
      gr_res[K+1] <- -2*sum(w*apply( rr2 * (-1/s + 1*(nodes.minus.XB)^2/s^3), 2, sum)/colSums(rr2))
      return(gr_res) 
    }

    # break up par into beta and residual components
    B <- par[1:K]
    s2 <- par[(K+1)]^2
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
    return( -2*( sum(w * (log(cs) - log(rr1m)) ) + sum(w) * log(dTheta) ) )
  }
}


interleaveMatrixRows <- function(A,B) {
  stopifnot(nrow(A) == nrow(B))
  stopifnot(ncol(A) == ncol(B))
  res <- matrix(NA, nrow=nrow(A)*2, ncol=ncol(A))
  for(i in 1:nrow(A)) {
    res[1+(i-1)*2,] <- A[i,]
    res[2+(i-1)*2,] <- B[i,]
  }
  return(res)
}


#' @importFrom Matrix Diagonal
#' @importFrom stats spline predict
fnCor <- function(Xb1, Xb2, s1, s2, w, rr1, rr2, nodes, fine=FALSE, fast = TRUE) {
  # distance between nodes
  dTheta <- nodes[2] - nodes[1]
  # get residuals
  numNodes <- length(nodes) 
  numUnits <- length(Xb1)
  nodes.minus.Xb1 <- t(Xb1 - t(matrix(nodes, nrow=numNodes, ncol=numUnits)))
  nodes.minus.Xb2 <- t(Xb2 - t(matrix(nodes, nrow=numNodes, ncol=numUnits)))
  function(par, returnPosterior=FALSE, plot=FALSE, newXb1=NULL, newXb2=NULL) {
    if(!is.null(newXb1)) {
      if(length(newXb1) != numUnits) {
        stop("The dimensions of Xb1 must remain the same.")
      }
      nodes.minus.Xb1 <- t(newXb1 - t(matrix(nodes, nrow=numNodes, ncol=numUnits)))
    }
    if(!is.null(newXb2)) {
      if(length(newXb2) != numUnits) {
        stop("The dimensions of Xb2 must remain the same.")
      }
      nodes.minus.Xb2 <- t(newXb2 - t(matrix(nodes, nrow=numNodes, ncol=numUnits)))
    }
    # par is propose correlation in Fisher-Z space, r is in correlation space
    r <- tanh(par)
    cov <- s1 * s2 * r
    Sigma <- matrix(c(s1^2, cov, cov, s2^2), ncol=2)
    detSigma <- s1^2*s2^2 - cov^2
    # invert Sigma
    SigmaInv <- 1/detSigma * matrix(c(s2^2, -1*cov, -1*cov, s1^2), ncol=2)
    # calculates vector of (x-mu)T %*% SigmaInv %*% (x-mu)
    # where x-mu is for an individual member is nodes.minux.Xb1/2, a vector of length 2
    rr_ <- vector(length=numUnits)
    # this is a 2D integral
    # for nodes in dimension 1
    mat <- matrix(NA, nrow=numNodes, ncol=numNodes)
    if(returnPosterior) {
      # set constants to zero
      mui2 <- muj2 <- mui <- muj <- m2 <- rep(0, numUnits)
    }
    for(i in 1:numNodes) {
      # node differences on subscale i
      nodesi <- nodes.minus.Xb1[i, ]
      # student probability of the nodes on subscale i
      rr1i <- rr1[i, ]
      # for nodes in dimension 2
      for(j in 1:numNodes) {
        # node differences on subscale j
        nodesj <- nodes.minus.Xb2[j, ]
        # student probability of the node on subscale j
        rr2j <- rr2[j, ]
        # for each 
        epsilon <- rbind(nodesi, nodesj)
        # the next line works, these are equal: mvnResid[4] (next line) and
        #                                       t(epsilon[,4]) %*% SigmaInv %*% epsilon[,4]
        mvnResid <- apply(epsilon * (SigmaInv %*% epsilon), 2, sum)
        # check a single entry with this code:
        # a <- c(nodes.minus.Xb1[i,k], nodes.minus.Xb2[j,k])
        # t(a) %*% SigmaInv %*% a - mvnResid[k] # = 0, they are the same
        # rr_:
        # put 1/sqrt(detSigma) into exp because it can be incredibly large while mvnResid is
        # incredibly small, making the entire expression resolve to zero.
        rr_ij <- rr1i * rr2j * (2/pi) * exp(-1/2 * (log(detSigma) + mvnResid))
        if(returnPosterior) {
          mui <- mui + rr_ij * nodes[i]
          muj <- muj + rr_ij * nodes[j]
          mui2 <- mui2 + rr_ij  * nodes[i]^2
          muj2 <- muj2 + rr_ij  * nodes[j]^2
          m2  <- m2  + rr_ij * nodes[i] * nodes[j]
        }
        rr_ <- rr_ + rr_ij
        mat[i,j] <- sum(rr_ij)
      }
    }
    if(plot) {
      plot(mat, main=paste0("r=",round(r,4)))
    }
    cs <- pmax(.Machine$double.eps, rr_)

    finef <- min(floor( (abs(par))^1.5 ), 10) # cap at 10 before factor

    if(returnPosterior){ # & finef <= 1) {
      if(finef <= 1) {
        mu1 <- mui/rr_ # mean of var1
        mu2 <- muj/rr_ # mean of var2
        Cov <- (m2/rr_ - mu1*mu2) # covariance
        mu12 <- mui2/rr_ # mean of mu1^2
        mu22 <- muj2/rr_ # mean of mu2^2
        rho <- Cov/sqrt( (mu1^2 - mu12) * (mu2^2 - mu22) )
        return(rho)
      } else {
        # reset sums to zero
        mui2 <- muj2 <- mui <- muj <- m2 <- m2 * 0
      }
    }
    if(!fine | finef <=1) {
      return(-2 * (sum(w * log(cs)) + 2 * sum(w) * log(dTheta)))
    }
    rr1f <- matrix(NA, nrow=nrow(rr1)*finef, ncol=ncol(rr1))
    nodesf <- nodes
    for(i in 2:finef) {
      nodesf <- c(nodesf, nodes + dTheta * (i-1)/finef)
    }
    nodesf <- sort(nodesf) - dTheta/2
    for(i in 1:ncol(rr1)) {
      rr1f[,i] <- spline(x=nodes, y=rr1[,i], n=length(nodesf), xmin=min(nodesf), xmax=max(nodesf))$y
    }
    rr2f <- matrix(NA, nrow=nrow(rr2)*finef, ncol=ncol(rr2))
    for(i in 1:ncol(rr1)) {
      rr2f[,i] <- spline(x=nodes, y=rr2[,i], n=length(nodesf), xmin=min(nodesf), xmax=max(nodesf))$y
    }
    nodes.minus.Xb1f <- t(Xb1 - t(matrix(nodesf, nrow=numNodes*finef, ncol=numUnits)))
    nodes.minus.Xb2f <- t(Xb2 - t(matrix(nodesf, nrow=numNodes*finef, ncol=numUnits)))
    rr2_ <- vector(length=numUnits)
    if(plot) {
      matf <- matrix(0, nrow=numNodes * finef, ncol=numNodes * finef)
    }
    for(i in 2:(numNodes-1)) {
      nodesfi <- nodes.minus.Xb1f[(i-1)*finef + 1:finef, ]
      rr1fi <- rr1f[(i-1)*finef + 1:finef, ]
      for(j in 2:(numNodes-1)) {
        if(any( mat[(i-1):(i+1),(j-1):(j+1)] * 1e4 >= max(mat)) | (mat[i,j] * 1e8 >= max(mat))) {
          # close to an important value, 
          rr2fj <- rr2f[(j-1)*finef + 1:finef, ]
          nodesfj <- nodes.minus.Xb2f[(j-1)*finef + 1:finef, ]
          for(fine_i in 1:finef) {
            for(fine_j in 1:finef) {
              epsilon <- rbind(nodesfi[fine_i, ], nodesfj[fine_j, ])
              mvnResid <- apply(epsilon * (SigmaInv %*% epsilon), 2, sum)
              if(fast) {
                rr_ij <- calcRrij(fine_i-1, fine_j-1, rr1fi, rr2fj, detSigma, mvnResid)
              } else {
                rr_ij <- rr1fi[fine_i,] * rr2fj[fine_j,] * (2/pi) * exp(-1/2 * (log(detSigma) + mvnResid))
              }
              if(returnPosterior) {
                ii <- (i-1)*finef+fine_i
                jj <- (j-1)*finef+fine_j
                mui <- mui + rr_ij * nodesf[ii]
                muj <- muj + rr_ij * nodesf[jj]
                mui2 <- mui2 + rr_ij  * nodesf[ii]^2
                muj2 <- muj2 + rr_ij  * nodesf[jj]^2
                m2  <- m2  + rr_ij * nodesf[ii] * nodesf[jj]
              }
              rr2_ <- rr2_ + rr_ij
              if(plot) {
                matf[(i-1)*finef+fine_i,(j-1)*finef+fine_j] <- sum(rr_ij)
              }
            }
          }
        } else {
          # approximately zero
        }
      }
    }
    if(plot) {
      plot(matf, main=paste0("r=",round(r,4)," fine=",finef))
    }
    # log requires a probability above zero, use the smallest double
    # as a floor
    cs2 <- pmax(.Machine$double.eps, rr2_)
    if(returnPosterior) {
      mu1 <- mui/rr2_ # mean of var1
      mu2 <- muj/rr2_ # mean of var2
      Cov <- (m2/rr2_ - mu1*mu2) # covariance
      mu12 <- mui2/rr2_ # mean of mu1^2
      mu22 <- muj2/rr2_ # mean of mu2^2
      rho <- Cov/sqrt( (mu1^2 - mu12) * (mu2^2 - mu22) )
#TODO: how should we handle zero division in calculating posterior. We get all three cases: 0/0, 0 /int, int/0
      rho[!is.finite(rho)] <- 0 # replace division by zero, we might want to make this more strict
      return(rho)
    }
    # likelihood of normal distribution for residuals nodes.minus.XB
    # aggregate likelihood, weight, multiply by -2 to make it a deviance
    # dTheta squared because there are two dimensions to the integral
    # simpler version of formula
    # return(-2*sum(w * (log(cs*dTheta^2))))
    # rearange to:
    return(-2 * (sum(w * log(cs2)) + 2 * sum(w) * log(dTheta/finef)))
  }
}

# graded response model
grm <- function (theta, d, score, a, D) {
  maxD <- length(d)   
  if(score == 0) {
    pr <- 1/(1 + exp(D*a*(theta - d[(score+1)])))
  }
  else if(score == maxD) {
    pr <- 1/(1 + exp(-D*a*(theta - d[score])))
  } else {
    pr <- 1/(1 + exp(D*a*(theta - d[(score+1)]))) - 1/(1 + exp(D*a*(theta - d[score])))
  }
  pr
}

# log of density of binomial where size=1, accepts x in [0,1] for partial credit
ldbinom3 <- function(x,pr) {
  return(x*log(pr) + (1-x)*log(1-pr))
}

# graded partial credit model
gpcm <- function (theta, d, score, a, D) {
  if(is.na(score)){
    return(NA)
  } 
  if (score > length(d)) {
    stop (paste0("Score of ", score," higher than maximum (", length(d),")"))
  }
  if(score <= 0) {
    stop (paste0("Score of ", score," lower than minimum (1)"))
  }
  Da <- D * a
  exp(sum(Da * (theta - d[1:score]))) / sum(exp(cumsum(Da * (theta - d))))
}

# helper, get the graded partial credit model likelihood
gpcmLikelihood <- function (theta, d, score, a, D=1.7) {
  exp(sum(D * (theta - d[1:score]))) / sum(exp(cumsum(D * (theta - d))))
}

# helper, get the graded response model likelihood
grmLikelihood <- function (theta, d, score, a, D=1.7) {
  maxD <- length(d)   
  if(score == 0) {
    pr <- 1/(1 + exp(D*a*(theta - d[(score+1)])))
  }
  else if(score == maxD) {
    pr <- 1/(1 + exp(-D*a*(theta - d[score])))
  } else {
    pr <- 1/(1 + exp(D*a*(theta - d[(score+1)]))) - 1/(1 + exp(D*a*(theta - d[score])))
  }
  pr
}
