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
