#' Draw plausible values (PVs) from an mml fit
#' @param x a fit from a call to \code{\link{mml}}
#' @param npv integer indicating the number of plausible values to draw
#' @param pvVariableNameSuffix suffix added to new PV variables after construct name and before the plausible value ID. For example, if there is a construct \code{math} and the suffix is the default \code{_dire}, then the fourth plausible value would have a column name, \code{math_dire4}.
#' @param stochasticBeta deprecated.
#' @param \dots additional parameters passed on to methods.
#' 
#' @details
#' 
#' When the argument passed to \code{stocasticBeta} is a data frame then each column is an element that will be used as a
#' regression coefficient for that index of the coefficients vector. The row index used for the nth PV will be the nth row.
#' 
#' @return the \code{stuDat} attribute of the \code{\link{mml}} fit object, \code{x}, with columns for the new plausible values.
#' 
#' @example man/examples/PV.R
#' @author Paul Bailey, Sun-joo Lee, and Eric Buehler 
#' @export
drawPVs <- function(x, npv, pvVariableNameSuffix="_dire", stochasticBeta=deprecated(), ...) {
  UseMethod("drawPVs")
}

#' @rdname drawPVs
#' @method drawPVs summary.mmlMeans
#' @importFrom lifecycle is_present deprecate_soft deprecated
#' @export
drawPVs.summary.mmlMeans <- function(x, npv=5L, pvVariableNameSuffix="_dire", stochasticBeta=deprecated(), ...) {
  if(is_present(stochasticBeta)) {
    deprecate_soft("3.0.2", "drawPVs(stochasticBeta)", details="the latent regression coefficients are always treated as stochastic now.")
  }
  npv <- round(npv)
  stopifnot(npv >= 1)
  oneD <- inherits(x, "summary.mmlMeans")
  p <- ncol(x$values$X)
  coefvec <- x$rawCoef[1:p]
  BigSigma <- x$VCtheta[1:p, 1:p, drop=FALSE]
  # generate one beta per PV
  cs <- chol(as.matrix(nearPD(BigSigma, ensureSymmetry=TRUE, doSym=TRUE)$mat))
  betas <- mvnfast::rmvn(npv, mu=coefvec, sigma=cs, isChol=TRUE)
  
  ts <- x$testScale
  k <- x$values$k
  for(i in 1:npv) {
    # if stochastic, update E
    vlsi <- x$values
    beta <- betas[i,]
    for(j in 1:k) {
      vlsi$Xb[,j] <- vlsi$X %*% beta[(j-1)*p+1:p]
      vlsi$E <- matrix(getOneD_E_sd(vlsi$Xb[,1], vlsi$rr1_list, x$rawCoef[p+1], vlsi$nodes)$E, ncol=1)
    }
    # for each unit, gather E, sigma, and draw PVs
    pvsi <- as.data.frame(drawOnePV(vlsi, intercept=ts$location, scale=ts$scale))
    colnames(pvsi) <- "score"
    pvsi <- cbind(id=x$values$id_vec, pv=i, pvsi)
    if(i == 1) {
      pvs <- pvsi
    } else {
      pvs <- rbind(pvs,pvsi)
    }
  }
  pvs <- as.data.frame(pivot_wider(pvs, id_cols=id, names_from=pv, 
                                          values_from="score", names_prefix=ts$test, names_sep="_dire_"))
  colnames(pvs)[colnames(pvs)!="id"] <- paste0(x$outcome, pvVariableNameSuffix, 1:npv)
  colnames(pvs)[colnames(pvs)=="id"] <- x$idVar
  stuDat <- merge(x$stuDat, pvs, by=x$idVar, all.x=TRUE, all.y=FALSE, order="x")
  return(stuDat)
}

#' @rdname drawPVs
#' @method drawPVs mmlMeans
#' @importFrom stats model.frame
#' @export
drawPVs.mmlMeans <- function(x, npv=5L,
                             pvVariableNameSuffix="_dire", stochasticBeta=deprecated(), ...) {
  stop("This method not implemented. Instead call summary on the mmlMeans object and call drawPVs with that.")
}

#' @rdname drawPVs
#' @method drawPVs summary.mmlCompositeMeans
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @export
drawPVs.summary.mmlCompositeMeans <- function(x, npv=20L, pvVariableNameSuffix="_dire", stochasticBeta=deprecated(), ...) {
  if(is_present(stochasticBeta)) {
    deprecate_soft("3.0.2", "drawPVs(stochasticBeta)", details="the latent regression coefficients are always treated as stochastic now.")
  }
  npv <- round(npv)
  stopifnot(npv >= 1)
  coefvec <- as.vector(t(x$rawCoef))
  BigSigma <- x$VCfull_theta
  # this could be nearly PD but not PD, so make sure that it is
  BigSigma <- Matrix::nearPD(BigSigma)$mat
  # generate one beta per PV
  betas <- mvnfast::rmvn(npv, mu=coefvec, sigma=BigSigma, isChol=FALSE)

  ts <- x$testScale

  cli_args <- list(total=npv,
                   format="{cli::pb_spin} Drawing plausible values {cli::pb_bar} ETA:{cli::pb_eta} Elapsed: {cli::pb_elapsed}")
  res <- do.call(cli_progress_bar, cli_args)
  for(i in 1:npv) {
    cli_progress_update()
    # if stochastic, update E
    k <- x$values$k
    vlsi <- x$values
    beta <- betas[i,]
    p <- ncol(vlsi$X)
    for(j in 1:k) {
      vlsi$Xb[,j] <- vlsi$X %*% beta[(j-1)*p+1:p]
      if(j > 1) {
        gc()
        # for each dimension, update E, this updates 1 too
        vlsi <- e_step_rect_grid(vlsi, j, 1, e_only=TRUE, addGr=FALSE)
      }
    }
    # for each unit, gather E, sigma, and draw PVs
    pvsi <- as.data.frame(drawOnePV(vlsi, intercept=ts$location, scale=ts$scale))
    colnames(pvsi) <- ts$subtest
    pvsi[,ts$test[1]] <- apply(t(pvsi) * ts$subtestWeight, 2, sum)
    pvsi <- cbind(id=x$ids,pv=i,pvsi)
    if(i == 1) {
      pvs <- pvsi
    } else {
      pvs <- rbind(pvs,pvsi)
    }
  }
  cli_progress_done()
  pvs <- as.data.frame(pivot_wider(pvs, id_cols=id, names_from=pv, 
                                          values_from=c(ts$subtest,"composite"), names_sep=pvVariableNameSuffix))
  colnames(pvs)[colnames(pvs)=="id"] <- x$idVar
  stuDat <- merge(x$stuDat, pvs, by=x$idVar, all.x=TRUE, all.y=FALSE, order="x")
  return(stuDat)
}

#' @rdname drawPVs
#' @importFrom MASS mvrnorm
#' @method drawPVs mmlCompositeMeans
#' @export
drawPVs.mmlCompositeMeans <- function(x, npv=20L, pvVariableNameSuffix="_dire", stochasticBeta=deprecated(), ...) {
  stop("This method not implemented. Instead call summary on the mmlCompositeMeans object and call drawPVs with that.")
}

#' @importFrom Matrix nearPD
nearPD2 <- function(X, tol=400, warn="") {
  eig <- eigen(X)
  if(min(eig$values) <= 0 || max(eig$values)/min(eig$values) >= 1/((.Machine$double.eps)^0.25)) {
    X <- nearPD(X, posd.tol=tol*sqrt(.Machine$double.eps))$mat
  }
  return(X)
}


drawOnePV <- function(values, intercept, scale, weights) {
  n <- values$n
  k <- values$k
  pvs <- matrix(0, nrow=n, ncol=k)
  for(i in 1:n) {
    Ei <- values$E[i,]
    sigma <- values$sigma_list[[i]]
    pvs[i,] <- mvnfast::rmvn(1, mu=Ei, sigma=sigma, isChol=TRUE)
  }
  pvs <- t(t(pvs) * scale + intercept)
  return(pvs)
}
