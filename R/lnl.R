### This is the marginal likelihood function to be maximized
# it is weighted if weight var is not null
#par - the score vector
#getIndLnL - logical, should the likelihood be calculated only for a subset
#X_ - the matrix of X values, defaults to all X , used for calculating group likelihoods
#i - list of indexes for individual likelihood (used to subset weights and r22)
#repweightVar  - the replicate weight for replicate variance estimation
#inside - TRUE if inside mix, otherwise linear and inside=FALSE
fn.regression <- function(X_, i=NULL, wv=NULL, rr1, stuDat, nodes, inside=TRUE, fast=TRUE) {
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
      # form prediction
      XB <- as.vector(X_ %*% B)
      # residual, per indivudal (outer sapply), per node (inner sapply)
      nodes.minus.XB <- t(t(matrix(nodes, nrow=length(nodes), ncol=length(XB))) - XB)
      # likelihood of normal distribution for residuals nodes.minus.XB
      rr2 <- rr1p * ((1/(sqrt(2 * pi * s2))) * exp(-((nodes.minus.XB)^2/(2 * s2))))
      # final parameter done numerically
      par_ <- par[(K+1)] + 1e-6
      # apply non-linear transform
      s2_ <- ifelse(par_ < 1, exp(par_ - 1), par_^2)
      s_ <- sqrt(s2_)
      # likelihood of normal distribution for residuals nodes.minus.XB
      rr2_ <- rr1p * ((1/(sqrt(2 * pi * s_^2))) * exp(-((nodes.minus.XB)^2/(2 * s_^2))))
      trr2mxb <- t(rr2 * nodes.minus.XB)
	    colSums_rr2 <- colSums(rr2)
      
      # Hess loop in rcpp
      if(fast) {
        H <- calcHess(K, rr2, rr2_, trr2mxb, X_, nodes.minus.XB, w, s2, s_)
        if(any(H %in% c(NA, NaN, Inf, -Inf))) {
          H <- slowHess(K, rr2, rr2_, trr2mxb, X_, nodes.minus.XB, w, s2, s_)
        }
      } else {
        H <- slowHess(K, rr2, rr2_, trr2mxb, X_, nodes.minus.XB, w, s2, s_)
      }
      # now variance of variance term
      gr0 <- -2*sum(w * colSums( rr2 * (-1/s + 1*(nodes.minus.XB)^2/s^3))/colSums_rr2)
      # break up par into beta and residual components
      par_ <- par[(K+1)] + 1e-6
      if(insd) {
        s2_ <- ifelse(par_ < 1, exp(par_ - 1), par_^2)
      } else {
        s2_ <- par_
      }
      s_ <- sqrt(s2_)
      # likelihood of normal distribution for residuals nodes.minus.XB
      rr2_ <- rr1p * ((1/(sqrt(2 * pi * s_^2))) * exp(-((nodes.minus.XB)^2/(2 * s_^2))))
      gr_ <- -2*sum(w*colSums( rr2_ * (-1/s_ + 1*(nodes.minus.XB)^2/s_^3))/colSums(rr2_))
      H[K+1,K+1] <- (gr_ - gr0)/1e-6
      return(H)
    }
    if(gr) {
      # form prediction
      XB <- as.vector(X_ %*% B)
      # residual, per indivudal (outer sapply), per node (inner sapply)
      nodes.minus.XB <- t(t(matrix(nodes, nrow=length(nodes), ncol=length(XB))) - XB)
      # likelihood of normal distribution for residuals nodes.minus.XB
      rr2 <- rr1p * ((1/(sqrt(2 * pi * s2))) * exp(-((nodes.minus.XB)^2/(2 * s2))))
      gr_res <- rep(0, length(par))
      trr2mxb <- t(rr2 * nodes.minus.XB)
      denom <- colSums(rr2)
      for(i in 1:K) {
        if(fast) {
          gr_res[i] <-  grSum2(w, trr2mxb, X_, i, s2, denom)
          if(gr_res[i] %in% c(NA, Inf, -Inf, NaN)) {
            gr_res[i] <- -2 * sum(w * rowSums( trr2mxb * X_[,i] / s2)/denom)
          }
        } else {
          gr_res[i] <- -2 * sum(w * rowSums( trr2mxb * X_[,i] / s2)/denom)
        }
      }
      gr_res[K+1] <- -2* sum(w * colSums( rr2 * (-1/s + 1*(nodes.minus.XB)^2/s^3))/denom)
      if(par[(K+1)] < 1) {
        gr_res[K+1] <- gr_res[K+1] * 0.5 * exp(0.5*(par[(K+1)]-1)) 
      }
      gr_res[is.na(gr_res)] <- 0
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
  function(par, returnPosterior=FALSE, newXb1=NULL, newXb2=NULL, gr=FALSE, hess=FALSE, superFine=0) {
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
    if(returnPosterior & gr) {
      stop("cannot return posterior gradient")
    }
    if(hess & gr) {
      stop("can only return one of Hessian or gradient")
    }
    # will we go fine
    finef <- min(floor( (abs(par))^1.5 ), 10) # cap at 10 before factor
    if(superFine > 0) {
      finef <- superFine
    }
    dofine <- finef > 1
    # par is propose correlation in Fisher-Z space, r is in correlation space
    r <- tanh(par)
    cov <- s1 * s2 * r
    Sigma <- matrix(c(s1^2, cov, cov, s2^2), ncol=2)
    detSigma <- s1^2*s2^2 - cov^2
    if(gr) {
      dcov <- s1 * s2 * 1/(cosh(par)^2)
      dDetSigma <- -2 * cov * s1 * s2 * (1/ (cosh(par)^2) )
    }
    if(hess) {
      dcov <- s1 * s2 * 1/(cosh(par)^2)
      ddcov <- -2 * s1 * s2 * sinh(par)/(cosh(par)^3)
      dDetSigma <- -2 * cov * s1 * s2 * (1/ (cosh(par)^2) )
      ddDetSigma <- -2 * dcov * s1 * s2 * (1/ (cosh(par)^2) ) + 4 * cov * s1 * s2 * sinh(par)/(cosh(par)^3)
    }
    #
    # calculates vector of (x-mu)T %*% SigmaInv %*% (x-mu)
    # where x-mu is for an individual member is nodes.minux.Xb1/2, a vector of length 2
    rr_ <- vector(length=numUnits)
    if( (gr | hess) & !dofine) {
      drr_ <- rr_ # another container for derivatives
      if(hess) {
        ddrr_ <- rr_ # another container for derivatives
      }
    }
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
        # the next line works, these are equal: mvnResid[4] (next line) and
        #                                       t(epsilon[,4]) %*% SigmaInv %*% epsilon[,4]
        # also equivalent to:  apply(epsilon * (SigmaInv %*% epsilon), 2, sum)
        # these use this sigma inverse, though solve could be used as well
        #         SigmaInv <- 1/detSigma * matrix(c(s2^2, -1*cov, -1*cov, s1^2), ncol=2)
        # and this defintion of epsilon
        #         epsilon <- rbind(nodesi, nodesj)
        mvnResid <- (1/detSigma) * (s1^2 * nodesi^2 - 2 * cov * nodesi * nodesj + s2^2 * nodesj^2)
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
        if(gr & !dofine){
          # mvnresid gradient
          mvnResidNum <- (s1^2 * nodesi^2 - 2*cov * nodesi*nodesj + s2^2 * nodesj^2)      
          dMvnResidNum <- - 2* dcov * nodesi * nodesj
          dMvnResid <- (1/detSigma^2) * (detSigma *  dMvnResidNum  - mvnResidNum * dDetSigma)
          # rr_ij gradient
          drr_ij <- -(1/2) * rr_ij * (dDetSigma/detSigma + dMvnResid)
          drr_ <- drr_ + drr_ij
        }
        if(hess & !dofine) {
          mvnResidNum <- (s1^2 * nodesi^2 - 2*cov * nodesi*nodesj + s2^2 * nodesj^2)      
          dMvnResidNum <- - 2* dcov * nodesi * nodesj
          ddMvnResidNum <- - 2* ddcov * nodesi * nodesj
          dMvnResid <- (1/detSigma^2) * (detSigma *  dMvnResidNum  - mvnResidNum * dDetSigma)
          num <- (detSigma *  dMvnResidNum  - mvnResidNum * dDetSigma)
          dnum <- (dDetSigma *  dMvnResidNum + detSigma *  ddMvnResidNum  - dMvnResidNum * dDetSigma - mvnResidNum * ddDetSigma)
          ddMvnResid <- (1/detSigma^4) * (detSigma^2 * dnum - num * 2 * detSigma * dDetSigma)
          # rr_ij gradient
          drr_ij <- -(1/2) * rr_ij * (dDetSigma/detSigma + dMvnResid)
          drr_ <- drr_ + drr_ij
          ddrr_ij <- -(1/2) * rr_ij * ((detSigma*ddDetSigma-dDetSigma^2)/detSigma^2 + ddMvnResid) -(1/2) * drr_ij * (dDetSigma/detSigma + dMvnResid)
          ddrr_ <- ddrr_ + ddrr_ij
        }
        mat[i,j] <- sum(rr_ij)
      }
    } # end for(i in 1:numNodes)
    cs <- pmax(.Machine$double.eps, rr_)
    if( (gr | hess)  & !dofine) {
      dcs <- ifelse(rr_ < .Machine$double.eps, 0, drr_)
      if(hess) {
        ddcs <- ifelse(rr_ < .Machine$double.eps, 0, ddrr_)
      }
    }

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
      if(hess) {
        return(-2 * sum(w * (cs*ddcs - dcs^2)/(cs^2)))
      }
      if(gr) {
        return(-2 * sum(w * dcs/cs))
      }
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
    if(gr | hess) {
      drr2_ <- rr2_
      if(hess) {
        ddrr2_ <- rr2_ 
      }
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
            fnodesi <- nodesfi[fine_i, ]
            for(fine_j in 1:finef) {
              fnodesj <- nodesfj[fine_j, ]
              mvnResid <- (1/detSigma) * (s1^2 * fnodesi^2 - 2*cov * fnodesi*fnodesj + s2^2 * fnodesj^2)
              if(fast) {
                rr_ij <- calcRrij(fine_i-1, fine_j-1, rr1fi, rr2fj, detSigma, mvnResid)
              } else {
                rr_ij <- rr1fi[fine_i,] * rr2fj[fine_j,] * (2/pi) * exp(-1/2 * (log(detSigma) + mvnResid))
              }
              if(gr) {
                # mvnresid gradient
                mvnResidNum <- (s1^2 * fnodesi^2 - 2*cov * fnodesi*fnodesj + s2^2 * fnodesj^2)      
                dMvnResidNum <- - 2* dcov * fnodesi * fnodesj
                dMvnResid <- (1/detSigma^2) * (detSigma *  dMvnResidNum  - mvnResidNum * dDetSigma)
                # rr_ij gradient
                drr2_ij <- -(1/2) * rr_ij * (dDetSigma/detSigma + dMvnResid)
                drr2_ <- drr2_ + drr2_ij
              }
              if(hess) {
                mvnResidNum <- (s1^2 * fnodesi^2 - 2*cov * fnodesi*fnodesj + s2^2 * fnodesj^2)      
                dMvnResidNum <- - 2* dcov * fnodesi * fnodesj
                ddMvnResidNum <- - 2* ddcov * fnodesi * fnodesj
                dMvnResid <- (1/detSigma^2) * (detSigma *  dMvnResidNum  - mvnResidNum * dDetSigma)
                num <- (detSigma *  dMvnResidNum  - mvnResidNum * dDetSigma)
                dnum <- (dDetSigma *  dMvnResidNum + detSigma *  ddMvnResidNum  - dMvnResidNum * dDetSigma - mvnResidNum * ddDetSigma)
                ddMvnResid <- (1/detSigma^4) * (detSigma^2 * dnum - num * 2 * detSigma * dDetSigma)
                drr_ij <- -(1/2) * rr_ij * (dDetSigma/detSigma + dMvnResid)
                drr2_ <- drr2_ + drr_ij
                ddrr_ij <- -(1/2) * rr_ij * ((detSigma*ddDetSigma-dDetSigma^2)/detSigma^2 + ddMvnResid) -(1/2) * drr_ij * (dDetSigma/detSigma + dMvnResid)
                ddrr2_ <- ddrr2_ + ddrr_ij
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
            }
          }
        } else {
          # approximately zero
        }
      }
    }
    # log requires a probability above zero, use the smallest double
    # as a floor
    cs2 <- pmax(.Machine$double.eps, rr2_)
    if(gr) {
      dcs2 <- ifelse(rr2_ <= .Machine$double.eps, 0, drr2_)
      return(-2 * sum(w * dcs2/cs2))
    }
    if(hess) {
      dcs2 <- ifelse(rr2_ <= .Machine$double.eps, 0, drr2_)
      ddcs2 <- ifelse(rr2_ < .Machine$double.eps, 0, ddrr2_)
      return(-2 * sum(w * (cs2*ddcs2-dcs2^2)/(cs2^2)))
    }
    if(returnPosterior) {
      mu1 <- mui/rr2_ # mean of var1
      mu2 <- muj/rr2_ # mean of var2
      Cov <- (m2/rr2_ - mu1*mu2) # covariance
      mu12 <- mui2/rr2_ # mean of mu1^2
      mu22 <- muj2/rr2_ # mean of mu2^2
      rho <- Cov/sqrt( (mu1^2 - mu12) * (mu2^2 - mu22) )
      rho[!is.finite(rho)] <- 0 # replace division by zero
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

slowHess <- function(K, rr2, rr2_, trr2mxb, X_, nodes.minus.XB, w, s2, s_) {
  colSums_rr2_ <- colSums(rr2_)
  H <- matrix(0, nrow=K+1, ncol=K+1)
  denom <- colSums(rr2)
  denom2 <- denom^2
  for(i in 1:K) {
    fi <- t( trr2mxb * X_[,i]/s2)
    num <- colSums(fi)
    for(j in i:K) {
      fj <- t( trr2mxb * X_[,j]/s2)
      numPrime <- colSums(t(t(rr2) * X_[,i] * X_[,j] / s2) - t( t(fj * nodes.minus.XB) * X_[,i])/s2)
      denomPrime <- colSums(fj)
      H[j,i] <- H[i,j] <- 2*sum(w*(numPrime * denom + num * denomPrime) / denom2)
    }
    # do numerical cross partial of variance term, taking care of the possibility of exp(x-1) transform
    gr0 <- -2 * sum(w * colSums(fi)/denom)
    gr_ <- -2 * sum(w * colSums( t( t((rr2_ * nodes.minus.XB)) * X_[,i]/s_^2))/colSums_rr2_)
    H[K+1,i] <- H[i, K+1] <- (gr_ - gr0)/1e-6
  }
  return(H)
}
