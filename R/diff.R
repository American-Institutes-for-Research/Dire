# This is a helper function to calculate the jacobian of the function with respect to all pars,
# note we have not fixed s
gradgradT <- function(location, ii, X_subset, weightVar, rr1=rr1, stuDat=stuDat, nodes=nodes){
  if(!inherits(rr1, "matrix")) {
    gri <- c()
    k <- length(location)/length(rr1)
    for(i in 1:length(rr1)) {
      fn2 <- fn.regression(X_=X_subset, i=ii, wv=weightVar, rr1=rr1[[i]], stuDat=stuDat, nodes=nodes, inside=FALSE)
      loc <- location[1:k]
      gri <- c(gri, -1/2*getGrad(fn2, loc))
      location <- location[-(1:k)]
    }
    return(gri %*% t(gri))
  }
  fn2 <- fn.regression(X_=X_subset, i=ii, wv=weightVar, rr1=rr1, stuDat=stuDat, nodes=nodes, inside=FALSE)
  # deviance is 2*lnl, so -1/2*grad of deviance is grad of lnl
  gr <- -1/2*getGrad(fn2, location)
  return(gr %*% t(gr))
}

gradInd <- function(location, ii, X_subset, weightVar, rr1=rr1, stuDat=stuDat, nodes=nodes){
  if(!inherits(rr1, "matrix")) {
    gri <- c()
    k <- length(location)/length(rr1)
    for(i in 1:length(rr1)) {
      fn2 <- fn.regression(X_=X_subset, i=ii, wv=weightVar, rr1=rr1[[i]], stuDat=stuDat, nodes=nodes, inside=FALSE)
      loc <- location[1:k]
      gri <- c(gri, -1/2*getGrad(fn2, loc))
      location <- location[-(1:k)]
    }
    return(gri)
  }
  fn2 <- fn.regression(X_=X_subset, i=ii, wv=weightVar, rr1=rr1, stuDat=stuDat, nodes=nodes, inside=FALSE)
  # deviance is 2*lnl, so -1/2*grad of deviance is grad of lnl
  gr <- -1/2*getGrad(fn2, location)
  return(gr)
}


# from WeMix
# author: Paul Bailey
getHessian <- function(func, x, inputs=1:length(x), f0=NULL, h = (abs(x) + 1) * (.Machine$double.eps)^0.25) {
  # use a wide net
  k <- length(inputs)
  if(is.null(f0)){
    f0 <- func(x)
  }
  if(length(f0) > 1) {
    return(lapply(1:length(f0),
                  function(i) {
                    f2 <- function(x) { 
                      func(x)[i]
                    }
                    getHessian(f2, x, inputs, f0=f0[i])
                  }))
  }
  # do not let h get smaller than the fourth root of the machine precision
  # derivatives at a larger x (fp, for f plus h) and smaller x (fm for f minus h)
  hess <- matrix(NA, nrow=k, ncol=k)
  for(i in 1:k) {
    for(j in i:k) {
      xpp <- xmm <- xmp <- xpm <- xm <- xp <- x
      if(i!=j) {
        # e.g. xpm is x with (x[1], x[2], ..., x[i] + h[i], ..., x[j] - h[j], ...)
        # so the first index (p, in the example) indicates x[i] is incrimented
        # by h[i] while the second index (m in the example) indicates x[j]
        # is decrimented by h[j].
        # xpp = (x[1], ..., x[i] + h[i], ..., x[j] + h[j], ...)
        xpp[inputs[i]] <- x[inputs[i]] + h[inputs[i]]
        xpp[inputs[j]] <- x[inputs[j]] + h[inputs[j]]
        # xmm = (x[1], ..., x[i] - h[i], ..., x[j] - h[j], ...)
        xmm[inputs[i]] <- x[inputs[i]] - h[inputs[i]]
        xmm[inputs[j]] <- x[inputs[j]] - h[inputs[j]]
        # xmp = (x[1], ..., x[i] + h[i], ..., x[j] - h[j], ...)
        xpm[inputs[i]] <- x[inputs[i]] + h[inputs[i]]
        xpm[inputs[j]] <- x[inputs[j]] - h[inputs[j]]
        # xmp = (x[1], ..., x[i] - h[i], ..., x[j] + h[j], ...)
        xmp[inputs[i]] <- x[inputs[i]] - h[inputs[i]]
        xmp[inputs[j]] <- x[inputs[j]] + h[inputs[j]]
        # find \partial^2 f(x) / \partial x[i] \partial x[j]
        # Start with u(x) = \partial f(x) / \partial x[i] = ( f(xp_i) - f(xm_i)) / (2h_i)
        # and apply that same formula to hess[i,j] = \partial u(x) / \partial x[j]
        hess[j,i] <- hess[i,j] <- (func(xpp) - func(xpm) - func(xmp) + func(xmm)) / (4* h[inputs[i]]*h[inputs[j]])
      } else {
        # xp = (x[1], ..., x[i] + h[i], ...)
        xp[inputs[i]] <- x[inputs[i]] + h[inputs[i]]
        # xm = (x[1], ..., x[i] - h[i], ...)
        xm[inputs[i]] <- x[inputs[i]] - h[inputs[i]]
        # xpp = (x[1], ..., x[i] + 2*h[i], ...)
        xpp[inputs[i]] <- x[inputs[i]] + 2*h[inputs[i]]
        # xmm = (x[1], ..., x[i] - 2*h[i], ...)
        xmm[inputs[i]] <- x[inputs[i]] - 2*h[inputs[i]]
        # the above formula would suggest
        # (func(xpp) - 2*f0 + func(xmm)) / (4*h[inputs[i]]^2)
        # but this is more o(h^4) as opposed to o(h^2) above,
        # and the diagonal is critical, so use it
        hess[i,i] <- (-1*func(xmm) + 16*func(xp) - 30*f0 + 16*func(xm)-1*func(xpp)) / (12*h[inputs[i]]^2)
      }
    }
  }
  return(hess)
}

# This function returns the first derivative of func with
# respect to each element of x
# from WeMix
# Author: Paul Bailey
getGrad <- function(func, x, inputs=1:length(x), highAccuracy=TRUE) {
  # use a wide net
  k <- length(inputs)
  # calculate the gradient
  grad <- vector(mode="numeric", length=k)
  h <- (abs(x) + 1) * sqrt(.Machine$double.eps)
  for(i in 1:length(inputs)) {
    xpp <- xmm <- xm <- xp <- x
    xp[inputs[i]] <- x[inputs[i]] + h[inputs[i]]
    xpp[inputs[i]] <- x[inputs[i]] + 2*h[inputs[i]]
    xm[inputs[i]] <- x[inputs[i]] - h[inputs[i]]
    xmm[inputs[i]] <- x[inputs[i]] - 2*h[inputs[i]]
    # forward difference, does not work near max (imagine if h takes you past the max)
    # (func(xp) - f0)/h[i]
    # central differences work better
    if(highAccuracy) {
      # this is the o(h^4) five point stencil method (the middle points has weight 0)
      ff <- ( (fxmm <- func(xmm)) - 8*func(xm) + 8*func(xp) - func(xpp))
      if( abs(80*fxmm) * .Machine$double.eps > abs(ff) ) {
        # too close, set to zero
        grad[i] <- 0
      } else {
        # not too close to zero, use normal equation
        grad[i] <- ff / (12*h[inputs[i]])
      }
    } else {
      ff <- (func(xp) - (fxm <- func(xm)))
      if( abs(10*fxm) * .Machine$double.eps > abs(ff) ) {
        grad[i] <- 0
      }
      # this is the o(h^2) central first difference (again, no middle point needed)
      grad[i] <- ff / (2*h[inputs[i]])
    }
  }
  return(grad)
}
