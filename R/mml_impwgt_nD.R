# Cholesky that does not throw an error but instead uses nearPD to modify non-PD matricies to be PD
# does nothing to a 2x2 correlation matrix with 0.5 correlation between the two variables
# test:
# b <- chol_robust(matrix(c(1,0.5,0.5,1),ncol=2)); print(t(b) %*% b)
# changes a non-PD correlation matrix to have correlation of 0.9998
# test:
# b <- chol_robust(m <- matrix(c(1,1.01,1.01,1),ncol=2)); print(t(b) %*% b, digits=16)
# returns Cholesky of m
chol_robust <- function(m) {
  tryCatch(res <- chol(m), error=function(e) {
    res <<- chol(Matrix::nearPD(m, keepDiag=TRUE, base.matrix=TRUE, doSym=TRUE, eig.tol=1e-4, do2eigen=TRUE, posd.tol=1e-4)$mat)
  })
  return(res)
}

#'@importFrom cli cli_progress_along cli_progress_done
#'@importFrom mvnfast rmvt rmvn dmvn dmvt
impwgt_nD_add_ints <- function(values, nint=100, updateGr=TRUE) {
  n <- values$n
  values$nint <- nint
  Sigma <- values$Sigma
  cSigma <- chol_robust(Sigma)
  Sinv <- solve(Sigma)
  dim <- ncol(values$xb)
  k <- values$k
  values$ints <- list()
  if(updateGr) {
    values$gr <- matrix(0, nrow=n, ncol=k)
    Xb <- as.matrix(values$Xb)
  }
  cli_args <- list(x=1:n,
                   format="{cli::pb_spin} Evaluating gradients {cli::pb_bar} ETA:{cli::pb_eta} Elapsed: {cli::pb_elapsed}")
  res <- lapply(do.call(cli_progress_along, cli_args), function(i) {
    sigma <- values$sigma_list[[i]]
    E <- values$E[i,]
    ints <- mvnfast::rmvt(nint, mu=E, sigma=sigma, df=4, isChol=TRUE)
    #wgts <- 1/mnvfast::dmvt(nint, mu=E, sigma=sigma, df=4, isChol=TRUE)
    ints <- pmin(ints, pmax(ints, min(values$nodes))) # make sure none are off the basis for the splines
    values$ints[[i]] <<- ints
    prvec <- 1
    # multiply in all of the likelihoods
    for(ii in 1:k) {
      # this pmax makes sure that the spline never goes non-positive. Clearly, a negative number would be incorrect
      # and a zero value would prevent anything in this dimension from being used
      prvec <- prvec * pmax(.Machine$double.eps, values$rr1_splines[[i]][[ii]](ints[,ii]))
    }
    prvec <- prvec * mvnfast::dmvn(ints, mu=values$Xb[i,], sigma=cSigma, isChol=TRUE) /
                     mvnfast::dmvt(ints, mu=E, sigma=sigma, df=4, isChol=TRUE)
    # make gradient
    if(updateGr) {
      resids <- t(t(ints) - Xb[i,])
      for(ii in 1:k) {
        Xii <- rep(0, k)
        Xii[ii] <- 1
        values$gr[i,] <<- apply(prvec * resids %*% Sinv,2,mean)/mean(prvec)
      }
    }
    values$prvec[[i]] <<- prvec
  })
  cli_progress_done()
  return(values)
}

# uses: n, rrlist
#' @importFrom stats splinefun
add_rr1_splines <- function(values) {
  k <- length(values$rr1_list)
  values$rr1_splines <- lapply(1:values$n, function(i) {
    spline <- list()
    for(ii in 1:k) {
      spline[[ii]] <- splinefun(x=values$nodes, y=values$rr1_list[[ii]][,i])
    }
    return(spline)
  })
  return(values)
}

# uses: Sigma, xb, var_vec, cov_vec
# updates: E, var_vec, cov_vec
e_step_impwgt_nd <- function(values) {
  Sigma <- values$Sigma
  Sinv <- solve(Sigma)
  k <- values$k
  n <- values$n
  Xb <- as.matrix(values$Xb)
  res <- lapply(1:n, function(i) {
    sigma <- values$sigma_list[[i]]
    Ei <- values$E[i,]
    ints <- values$ints[[i]]
    prvec <- values$prvec[[i]]
    mm <- matmap_factory(k)
    new_sigma <- sigma
    for(ii in 1:k) {
      values$E[i,ii] <<- eii <- sum(ints[,ii] * prvec) / sum(prvec) 
      values$cov[i,ii] <<- varii <- sum((ints[,ii] - eii)^2 * prvec) / sum(prvec) 
      new_sigma[ii,ii] <- varii
      if(ii > 1) {
        for(jj in 1:(ii-1)) {
          cov_ij <- sum((ints[,ii] - eii)*(ints[,jj] - values$E[i,jj]) * prvec) / sum(prvec) 
          values$cov[i,mm(ii,jj)] <<- cov_ij
          new_sigma[ii,jj] <- new_sigma[jj,ii] <- cov_ij
        }
      }
    }
    new_sigma <- chol_robust(new_sigma)
    values$sigma_list[[i]] <<- new_sigma
  })
  return(values)
}
