#' @method print mmlMeans
#' @export
print.mmlMeans <- function(x, ...){
  co <- coef(x)
  print(co, ...)
}

#' @method predict mmlMeans
#' @export
predict.mmlMeans <- function(object, newData=NULL, ...) {
  if(!is.null(newData)) {
    trms <- delete.response(terms(object$formula))
    m <- model.frame(trms, data=newData, drop.unused.levels=TRUE)
    X <- model.matrix(trms, m, contrasts.arg = object$contrasts, xlev = object$xlevels)
  } else {
    X <- object$getX(object$stuDat)
  }
  b <- object$coef
  pred <- X %*% b[-length(b)]
  names(pred) <- rownames(X)
  return(pred)
}

#' @method predict mmlCompositeMeans
#' @export
predict.mmlCompositeMeans <- function(object, newData=NULL, ...) {
  if(is.null(newData)) {
    newData <- object$stuDat
  }
  res <- data.frame(id=rownames(newData), stringsAsFactors=FALSE)
  for(i in 1:length(object$resl)) {
    newObj <- list(formula = object$formula,
                   contrasts = object$contrasts[[i]],
                   xlevels = object$xlevels[[i]],
                   coef = object$coefficients[i, , drop=TRUE]
                   )
    res[ , paste0("pred", i)] <- predict.mmlMeans(newObj, newData=newData)
  }
  return(res)
}

#' @method print mmlCompositeMeans
#' @export
print.mmlCompositeMeans <- function(x, ...){
  co <- coef(x)
  print(co, ...)
}

mml.remap.coef <- function(coefficients, location, scale, noloc=FALSE) {
  coefficients <- coefficients * scale
  if(!noloc) {
    coefficients[grepl("(Intercept)", names(coefficients))] <- coefficients[names(coefficients) == "(Intercept)"] + location
  }
  return(coefficients)
}


#' @method summary mmlMeans
#' @export
summary.mmlMeans <- function(object, gradientHessian=TRUE,
                             varType=c("consistent", "robust", "cluster", "Taylor"),
                             clusterVar=NULL, jkSumMultiplier=1, # cluster
                             repWeight=NULL, # replicate
                             strataVar=NULL, PSUVar=NULL, singletonFix=c("drop", "use mean"),# Taylor
                             ...){
  if(!missing(gradientHessian)) {
    warning("gradientHessian has been deprecated. It is always TRUE now, meaning that the outer product of Hessians are used to estimate the Hessian matrix.")
  }
  object$rawCoef <- object$coef
  # check/fix varType argument
  # for EdSurvey 2.7.1 compatibility
  varType0 <- c("consistent", "robust", "cluster", "replicate", "Taylor")
  if(length(varType) == length(varType0) && all(varType == varType0)) {
    varType <- "Taylor"
  }
  varType <- match.arg(varType)
  sumCall <- match.call()
  # keep these
  latentCoef <- object$coefficients
  if(is.null(strataVar)) {
    strataVar <- object$strataVar
  }
  if(is.null(PSUVar)) {
    PSUVar <- object$PSUVar
  }
  H_B_prime <- getIHessian.mmlMeans(object)
  singletonFix <- match.arg(singletonFix)
  if(varType=="consistent") {
    VC <- getVarConsistent(object, H_B_prime)
  }
  if(varType=="robust") {
    VC <- getVarRobust(object, H_B_prime)
  }
  if(varType=="cluster") {
    stuDat <- object$stuDat
    if(is.null(clusterVar)) {
      stop("You must define a valid clusterVar to use cluster variance estimation.")
    }
    if(length(clusterVar) != 1) {
      if("ClusterVar__" %in% colnames(stuDat)) {
        stop("Please rename the variable ", dQuote("ClusterVar__"), " on the ", dQuote("stuDat"), " argument.")
      }
      # paste together variables with colons
      # first, remove colons from existing variables, if necessary
      for(i in 1:length(clusterVar)) {
        if(inherits(stuDat[clusterVar], "character") && sd(nchar(stuDat[clusterVar])) > 0) {
          stuDat[clusterVar] <- gsub(":", "-", stuDat[clusterVar], fixed=TRUE)
        }
        if(inherits(stuDat[clusterVar], "factor")) {
          stuDat[clusterVar] <- gsub(":", "-", as.character(stuDat[clusterVar]), fixed=TRUE)
        }
      }
      stuDat$ClusterVar__ <- apply(stuDat[clusterVar], 1, function(x) { paste(x, collapse=":") } )
      clusterVar <- "ClusterVar__"
    }
    if(!clusterVar %in% colnames(stuDat)) {
      stop(paste0("Could not find clusterVar column named ", dQuote(clusterVar), " on ", dQuote("stuDat"), " data"))
    }
    if(length(unique(stuDat[,clusterVar])) <= 1) {
      stop("There must be more than one cluster for cluster variance estimation.")
    }
    VC <- getVarCluster(object, H_B_prime, clusterVar)
  }
  if(varType=="Taylor") {
    stuDat <- object$stuDat
    if(is.null(strataVar) | is.null(PSUVar)) {
      stop(paste0("the arguments ", dQuote("strataVar"), " and ", dQuote("PSUVar")," must be defined for varType ", dQuote("Taylor"), "."))
    }
    if(!strataVar %in% colnames(stuDat)) {
      stop(paste0("Could not find strataVar column named ", dQuote(strataVar), " on ", dQuote("stuDat"), " data."))
    }
    if(!PSUVar %in% colnames(stuDat)) {
      stop(paste0("Could not find strataVar column named ", dQuote(PSUVar), " on ", dQuote("stuDat"), " data."))
    }
    Tres <- getVarTaylor(object, H_B_prime, strataVar, PSUVar, singletonFix)
    VC <- Tres$VC
    dof <- Tres$dof
  }
  VCtheta <- VC
  VC <- VC * object$scale^2
  se <- sqrt(diag(VC))
  tval <- as.vector(coef(object)/se)

  TAB <- data.frame(Estimate = coef(object),
                    StdErr = se,
                    t.value = tval)
  if(exists("dof")) {
    TAB$dof <- dof
    TAB$"Pr(>|t|)" <- 2*(1-pt(abs(TAB$t.value), df=dof))
  }
  row.names(TAB) <- names(object$coefficients) 
  # get weighted obs, if relevant
  if(is.null(object$weightVar)) {
    wo <- NA
  } else {
    wo <- object$weightedObs
  }

  object$coefficients <- as.matrix(TAB)
  object <- c(object, list("summaryCall" = sumCall,
                           "latentCoef" = latentCoef,
                           "LL" = object$LogLik,
                           "VC" = VC,
                           "iHessian" = H_B_prime,
                           "weightedObs" = wo,
                           "VCtheta" = VCtheta))
  class(object) <- "summary.mmlMeans"
  return(object)
}

#' @importFrom stats printCoefmat
#' @method print summary.mmlMeans
#' @export
print.summary.mmlMeans <- function(x, ...){
  cat(paste0("Call:\n"))
  print(x$call)
  cat(paste0("Summary Call:\n"))
  print(x$summaryCall)
  cat("\n")
  cat("Summary:\n")
  cof <- x$coefficients
  cof1 <- cof[1:(nrow(cof)-1),,drop=FALSE]
  cof2 <- cof[nrow(cof),1:2,drop=FALSE]
  printCoefmat(cof1)
  cat("\n")
  cat("Residual Variance Estimate:\n")
  print(cof2)
  cat("\n")
  if(!is.na(x$Convergence)) {
    cat(paste0("Convergence = ", x$Convergence, "\n"))
  } 
  if(!is.na(x$iterations)) {
   cat(paste0("Iterations = ", x$iterations, "\n"))
  }
  if(!is.na(x$LogLik)) {
    cat(paste0("LogLike = ", round(x$LogLik,2), "\n"))
  }
  cat(paste0("Observations = ", x$obs, "\n"))
  if(!is.na(x$weightedObs)) {
    cat(paste0("Weighted observations = ", round(x$weightedObs,2), "\n"))
  }
  cat(paste0("location = ", x$location, " scale = ", x$scale, "\n"))
  if("SubscaleVC" %in% names(x$Summary)) {
    cat("\nEstimated subscale correlations:")
    print(round(x$Summary$SubscaleVC, 2))
  }
} 


#' @method summary mmlCompositeMeans
#' @export
#' @importFrom stats pt
summary.mmlCompositeMeans <- function(object, gradientHessian=TRUE,
                                      varType=c("Taylor", "Partial Taylor"),
                                      strataVar=NULL, PSUVar=NULL, singletonFix=c("use mean", "drop"),# Taylor
                                      verbose=0,
                                      nint=10000,
                                      JohnsonRust=FALSE,
                                      ...){
  sumCall <- match.call()
  if(!missing(gradientHessian)) {
    warning("gradientHessian has been deprecated. It is always TRUE now, meaning that the outer product of Hessians are used to estimate the Hessian matrix.")
  }
  # EdSurvey 2.7.1 compatibility
  varType0 <- c("consistent", "robust", "cluster", "replicate", "Taylor")
  if(length(varType) == length(varType0) && all(varType == varType0)) {
    varType <- "Taylor"
  }
  varType <- match.arg(varType)
  # get varType and singletonFix cleaned up
  singletonFix <- match.arg(singletonFix)
  if(is.null(strataVar)) {
    if(is.null(object$strataVar) & varType %in% c("Taylor", "Partial Taylor")) {
      stop(paste0("argument ", dQuote("strataVar"), " must be included in the ", dQuote("mml"), " or ", dQuote("summary"), " call."))
    }
    strataVar <- object$strataVar
  }
  if(is.null(PSUVar)) {
    if(is.null(object$PSUVar) & varType %in% c("Taylor", "Partial Taylor")) {
      stop(paste0("argument ", dQuote("PSUVar"), " must be included in the ", dQuote("mml"), " or ", dQuote("summary"), " call."))
    }
    PSUVar <- object$PSUVar
  }
  # first
  k <- nrow(object$coefficients) 
  p <- ncol(object$coefficients) - 1 # drop SD
  rawCoef <- object$coefficients[,1:p,drop=FALSE]
  # the H_B_prime0 and (later) VCfull are block diagonal matrixes
  # "block" is a list where each element is every index for that block
  block <- list()
  for(i in 1:k) {
    block <- c(block, list(1:p + p*(i-1)))
  }
  # build H_B_prime
  H_B_prime0 <- matrix(0, nrow=p*k, ncol=p*k)
  X <- object$resl[[1]]$getX(object$stuDat)
  for(i in 1:k) {
    if(verbose>0) {
      message(paste0("Finding Hessian block ",i," of ", k))
    }
    if(gradientHessian) {
      gr <- object$resl[[i]]$funcs$gr(object$coefficients[i,], sqrtW=TRUE)
      H_B_prime0[block[[i]], block[[i]]] <- accumulator_outer_fast(X*gr)
    } else {
      H_B_prime0[block[[i]], block[[i]]] <- object$resl[[i]]$funcs$hess(object$coefficients[i,])
    }
  }
  
  # get weighted obs, if relevant
  if(is.null(object$weightVar)) {
    wo <- NA
    w <- rep(1, nrow(object$stuDat))
  } else {
    wo <- object$weightedObs
    w <- object$stuDat[,object$weightVar]
  }
  if(varType=="Partial Taylor") {
    stop("Partial Taylor not implemented.")
  } # end if(varType=="Partial Taylor")
  if(varType=="Taylor") {
    dof <- Inf
    if(verbose > 0) {
      message("start Taylor")
    }
    # there are many steps to getting the Taylor series. First,
    # when k > 2 (and so to be general, always) we use the gradient from impwgtnD_add_ints
    # that requires numerical quadrature at arbitrary points, so will add splines for rr1
    if(verbose > 0) {
      cat("interpolating response probabilities\n")
    }
    vls <- add_rr1_splines(object$values) 
    if(verbose > 0) {
      cat("finding expected values and covariances\n")
    }
    # we also need first guess expected values, variances, and covariances so the distribution
    # is correctly centered for important weighting. These turn out to be very acruate despite being
    # calculated based only on 2 dimensions at once
    vls <- e_step_nD_by2D(vls) # find E, var, cov
    if(verbose > 0) {
      cat("forming unit-level covariance matricies\n")
    }
    # this is a convenience function that adds a per-unit sigma matrix
    vls <- add_sigma_list(vls) # construct list of individual-level covariance matrixes
    if(verbose > 0) {
      cat("calculating numerical approximation to unit-level gradient\n")
    }
    # calculate importance weighting sample
    # this allows the n/k-dimensional gradient to be calculated
    # we sill use the gradient calculated in this step for the Taylor series
    vls <- impwgt_nD_add_ints(vls, nint=nint, updateGr=TRUE)
    # this updates sigma (per individual) in an n-dimensional
    # way that leads to better condition numbers
    object$values <- e_step_impwgt_nd(vls)
    if(verbose > 0) {
      cat("calculating Taylor Covariance matrix\n")
    }
    # get Taylor VC matrix
    VC <- getVarTaylor2(object=object, H_B_prime = H_B_prime0,
                        strataVar = strataVar, PSUVar = PSUVar,
                        singletonFix = singletonFix, w=w, JohnsonRust=JohnsonRust)
    Vi <- VC$VC
    td <- object$testScale
    # VCfull is the scale weighted full VC matrix
    VCfull_theta <- Vi
    VCfull <- Vi
    for(i in 1:k) {
      for(j in 1:k) {
        if(i >= j) {
          # block i, j (potentially block i, i)
          VCfull[block[[i]], block[[j]]] <- td$scale[i] * td$scale[j] * VCfull[block[[i]], block[[j]]]
          if(i > j) {
            # symetric matrix, so also set block j, i
            VCfull[block[[j]], block[[i]]] <- t(VCfull[block[[i]], block[[j]]])
          }
        }      
      }
    }
  } # end if(varType="Taylor")
  
  if(is.null(colnames(VCfull))) {
    # renames
    colnames(H_B_prime0) <- rownames(H_B_prime0) <- colnames(VCfull) <- rownames(VCfull) <- 1:nrow(H_B_prime0)
    # first, make sure these are unique
    new_names <- c()
    for(i in 1:k) {
      new_names <- c(new_names, paste0(colnames(object$coefficients), "_", rownames(object$coefficients)[i]))
    }
    if(length(unique(new_names)) == length(new_names)) {
      for(i in 1:k) {
        rownames(H_B_prime0)[1:p + p*(i-1)] <- 
          colnames(H_B_prime0)[1:p + p*(i-1)] <-
          rownames(VCfull)[1:p + p*(i-1)] <-
          colnames(VCfull)[1:p + p*(i-1)] <- paste0(colnames(rawCoef), "_", rownames(rawCoef)[i])
      }
    }
  }
  # make a containder for scaled weights
  scaledCoef <- rawCoef
  # scale them
  for(i in 1:k) {
    scaledCoef[i, ] <- mml.remap.coef(rawCoef[i, ], td$location[i], td$scale[i])
  }
  # remap to weighted levels
  compositeCoef <- coef.mmlCompositeMeans(object)
  names(compositeCoef) <- colnames(object$coefficients)
  VCsmall <- matrix(0, nrow=p, ncol=p)
  for(i in 1:p) {
    vi <- rep(0, p*k)
    vi[i+((1:k)-1)*p] <- td$subtestWeight
    for(j in 1:p) {
      vj <- rep(0, p*k)
      vj[j+((1:k)-1)*p] <- td$subtestWeight
      VCsmall[i,j] <- t(vj) %*% VCfull %*% vi
    }
  }
  rownames(VCsmall) <- colnames(VCsmall) <- colnames(rawCoef)
  se <- sqrt(diag(VCsmall))
  se <- c(se, NA)
  tval <- as.vector(compositeCoef/se)
  TAB <- data.frame(Estimate = compositeCoef,
                    StdErr = se,
                    t.value = tval,
                    dof = rep(VC$simple_dof, p+1))
  if("dof" %in% colnames(TAB)) {
    TAB$"Pr(>|t|)" <- 2*(1-pt(abs(TAB$t.value), df=dof))
  }
  row.names(TAB) <- names(compositeCoef) 
  res <- object
  res$summaryCall <- sumCall
  res$coefficients <- as.matrix(TAB)
  res$VC <- VC
  res$iHessian <- -1*H_B_prime0
  res$weightedObs <- wo
  res$VCfull <- VCfull
  res$rawCoef <- rawCoef
  res$VCfull_theta <- VCfull_theta
  class(res) <- "summary.mmlCompositeMeans"
  return(res)
}



#' @method print summary.mmlCompositeMeans
#' @export
print.summary.mmlCompositeMeans <- function(x, ...){
  cat(paste0("Call:\n"))
  print(x$call)
  cat(paste0("Summary Call:\n"))
  print(x$summaryCall)
  cat("\n")
  cat("Summary:\n")
  cof <- x$coefficients
  cof1 <- cof[1:(nrow(cof)-1),,drop=FALSE]
  cof2 <- cof[nrow(cof),1:2,drop=FALSE]
  printCoefmat(cof1)
  cat("\n")
  cat("Residual Variance Estimate:\n")
  print(cof2)
  cat("\n")
  cat(paste0("Convergence = ", pasteItems(unique(x$Convergence)), "\n"))
  cat(paste0("Iterations = ", sum(x$iterations), "\n"))
  cat(paste0("observations = ", pasteItems(x$obs), "\n"))
  if(!all(is.na(x$weightedObs))) {
    cat(paste0("Weighted observations = ", pasteItems(round(x$weightedObs,2)), "\n"))
  }
} 


getIHessian.mmlMeans <- function(object) {
  lnlf <- object$funcs$lnl
  n <- nrow(object$stuDat)
  x0 <- object$coefficients
  K <- length(x0)
  gr_mat <- matrix(0, nrow=n, ncol=K)
  y0 <- lnlf(x0)
  h <- (abs(x0) + 1) * sqrt(.Machine$double.eps)
  if(is.null(object$weightVar)) {
    weightVar <- "one"
    object$stuDat[[ki]][["one"]] <- 1
  }
  w <- object$stuDat[[object$weightVar]]
  for(ki in 1:K) {
    xm2 <- xp2 <- xm1 <- xp1 <- x0
    xp1[ki] <- x0[ki] + h[ki]
    xp2[ki] <- x0[ki] + 2*h[ki]
    xm1[ki] <- x0[ki] - h[ki]
    xm2[ki] <- x0[ki] - 2*h[ki]
    yp1 <- lnlf(xp1, returnCS=TRUE)
    yp2 <- lnlf(xp2, returnCS=TRUE)
    ym1 <- lnlf(xm1, returnCS=TRUE)
    ym2 <- lnlf(xm2, returnCS=TRUE)
    # sqrt(w) because it gets squared in the outher and should only occur w times
    # -1 here 
    gr_mat[ , ki] <- -1 * sqrt(w) * (ym2 - 8*ym1 + 8*yp1 - yp2) / (12*h[ki])
  }
  true_outer <- function(x) {
    x %o% x
  }

  H <- accumulator_outer_cpp(gr_mat)
  # previous line is the same as this:
  #H <- Reduce("+", gg <- apply(gr_mat, 1, true_outer, simplify=FALSE))
  Wsum <- sum(w)
  H <- -1 * Wsum/(Wsum-1) * solve(nearPD2(H))
  return(H)
}


# estimate covariance matrix (Standard error)
# -1/2 turns a deviance into a likelihood
getVarConsistent <- function(object, H_B_prime) {
  #consistent  estimate
  return(-1 * H_B_prime )
}

getVarRobust <- function(object, H_B_prime) {
  X <- object$getX(object$stuDat)
  x_ind <- split(X, as.factor(1:nrow(X)), drop=FALSE)
  # dl/dBeta * dl/dBeta ' evaluated for each individual 
  object$stuDat$one <- 1
  vars <- lapply(1:nrow(X), FUN=function(i){gradgradT(location=object$coefficients,
                                                      ii=i,
                                                      X_subset=matrix(x_ind[[i]], nrow=1),
                                                      weightVar="one",
                                                      rr1=object$rr1,
                                                      stuDat=object$stuDat,
                                                      nodes=object$nodes)})
  V <- Reduce("+", vars)
  return(H_B_prime %*% V %*% H_B_prime)
}

getVarCluster <- function(object, H_B_prime, clusterVar) {
  stuDat <- object$stuDat
  X <- object$getX(object$stuDat)
  # dl/dBeta * dl/dBeta ' evaluated for each group  before being multipled and then summed
  # first get list of each group index 
  #this is important to ensure group_index and x_group+s are in the same order
  stuDat[[clusterVar]] <- factor(stuDat[[clusterVar]], levels=unique(stuDat[[clusterVar]])) 
  group_index <- lapply(levels(stuDat[[clusterVar]]), FUN=function(x) {
    which(stuDat[[clusterVar]]==x)
  })
  x_groups <- lapply(split(X, stuDat[[clusterVar]]), matrix, ncol=ncol(X))
  vars <- lapply(c(1:length(group_index)), FUN=function(group){
    gradgradT(location=object$coefficients,
              ii=group_index[[group]],
              X_subset = x_groups[[group]],
              weightVar=object$weightVar,
              rr1=object$rr1,
              stuDat=object$stuDat,
              nodes=object$nodes)
  })
  V <-  Reduce("+",vars)
  return(H_B_prime %*% V %*% H_B_prime)
} 

getVarReplicate <- function(object, H_B_prime, repWeight, jkSumMultiplier=1, returnVecs=FALSE) {
  X <- object$getX(object$stuDat)
  B0 <- object$coefficients
 
  B_j <- lapply(repWeight, FUN=function(x){ #restimate with each weight 
    fn2B <- fn.regression(X_=X, wv=x, rr1=object$rr1, stuDat=object$stuDat, nodes=object$nodes, inside=TRUE)
    return(robustOptim(fn2B, par0=object$coefficients, X=X)$par)
  })
  # for composite
  if(returnVecs) {
    return(B_j)
  }
  rep <- lapply(B_j, function(x){(x-B0) %*% t(x-B0)})
  return(jkSumMultiplier * Reduce("+", rep))
} 

getVarTaylor <- function(object, H_B_prime, strataVar, PSUVar,
                         singletonFix=c("drop", "use mean"),
                         returnVecs=FALSE, progress_bar=NULL) {
  #find PSU per strata, and warn about dropping strata with one PSU
  singletonFix <- match.arg(singletonFix)
  stuDat <- object$stuDat
  snames <- sort(unique(stuDat[[strataVar]]))
  strata <- lapply(snames, FUN=function(x) {
    list(strat=x,
         psu=sort(unique(stuDat[stuDat[,strataVar]==x, PSUVar])))
  })
  n_psu <- lapply(strata, function(x) { length(x$psu)}) 
  simple_dof <- sum(unlist(n_psu)) - sum(length(strata))
  if (any(n_psu==1)){
    if(singletonFix == "drop") {
      warning(paste0("Of the ", length(n_psu)," strata, ", sum(n_psu<2), " strata have only one PSU. All strata with only one PSU are excluded from variance estimation. See the ", dQuote("singletonFix"), " argument for other options."))
      strata <- strata[n_psu>1] #variance estimation can only happen for Strata with more than one PSU 
    }
    if(singletonFix == "use mean") {
      warning(paste0("Of the ", length(n_psu), " strata, ", sum(n_psu<2), " strata have only one PSU. All strata with only one PSU have their value compared to the mean. See the ", dQuote("singletonFix"), " argument for more details and other options."))
    }
  }
  #only keep snames in strata 
  snames <- lapply(strata, function(x){x$strat}) 
  #split X based on psu for access later
  if("X" %in% names(object)) {
    X <- object$X
  } else {
    X <- object$getX(object$stuDat) # this object is created in the caller with X in it
  }
  # loop through strata 
  str <- lapply(strata, FUN=function(st) {
    # number of PSUs in this stratum
    n_a <- length(st$psu)
    # PSU index for units in this stratum and PSU
    group_index <- lapply(st$psu, FUN=function(x) {
      which(stuDat[[PSUVar]]==x & stuDat[[strataVar]]==st$strat)
    })
    # extract data for this stratum
    X_strata <- X[stuDat[[strataVar]] %in% st$strat, , drop=FALSE]
    stuDat_strata <- stuDat[stuDat[[strataVar]] %in% st$strat, , drop=FALSE]
    #split up data by PSU, lapply after the split, making each one a matrix
    x_groups <- lapply(split(X_strata, stuDat_strata[[PSUVar]]), matrix, ncol=ncol(X_strata))
    #vector of scores per psu
    s_p <- lapply(c(1:length(group_index)), FUN=function(k) {
      gradInd(location=object$coefficients,
              ii=group_index[[k]],
              X_subset=x_groups[[k]],
              weightVar=object$weightVar,
              rr1=object$rr1,
              stuDat=object$stuDat,
              nodes=object$nodes)
    })
    st$s_p <- s_p 
    
    if(n_a > 1) {
      #average score across as psu in strata
      s_a_bar <- Reduce("+",s_p)/n_a
      #(s_p - s_a_bar)*(s_p - s_a_bar)' 
      if(returnVecs) {
        s <- lapply(s_p, FUN=function(s){
          (s - s_a_bar) 
        })
        st$s_p <- NULL
        names(s) <- st$psu
        st$s <- s
      } else {
        v_a <- lapply(s_p, FUN=function(s){
          (s - s_a_bar) %*% t(s - s_a_bar)
        })
        st$V_a <- (n_a/(n_a-1))*Reduce("+",v_a)
      }
    }
    return(st)
  }) # end lapply(strata, FUN=function(st) {
  # find zeros, zeroV_a
  zeroV_a <- NULL
  zeros <- NULL
  i <- 1
  while(is.null(zeroV_a)) {
    if(returnVecs) {
      # if this is not present and V_a is, then we don't need zeros, it's okay for it to stay NULL
      if(!is.null(str[[i]]$s)) {
        zeros <- 0 * str[[i]]$s[[1]]
      }
    }
    else {
      if(!is.null(str[[i]]$V_a)) {
        zeroV_a <- 0 * str[[i]]$V_a
      }
    }
    i <- i + 1
  }

  if(singletonFix %in% c("use mean")) {
    # get all strata's PSUs to find mean
    s_p_all <- list()
    for(stri in str) {
      s_p_all <- c(s_p_all, stri$s_p)
    }
    # this is the overall average, across strata, should be zero
    n_all <- length(s_p_all)
    s_a_barOverall <- Reduce("+", s_p_all)/n_all
  }
  if(returnVecs) {
    s <- lapply(str, FUN=function(st) {
      npsu <- length(st$psu)
      if(npsu > 1) {
        return(st$s)
      }
      if(singletonFix %in% c("use mean")) {
        s <- st$s_p[[1]] # there is only one element
        s <- 2*(s - s_a_barOverall) 
        return(s)
      }
      if(singletonFix=="drop") {
        return(zeros)
      }
    }) #end aggregate V_a lapply(str, FUN=function(st)
    names(s) <- snames
    return(list(vec=s, simple_dof=simple_dof))
  }
  # aggregate V_a
  V <- zeroV_a
  Vi <- zeroV_a
  num2 <- num <- rep(0, ncol(H_B_prime))
  numStrata <- 0
  for(i in seq_along(str)) {
    st <- str[[i]]
    npsu <- length(st$psu)
    if(npsu > 1) {
      numStrata <- numStrata + 1
      numi <- diag(Vii <- H_B_prime %*% st$V_a %*% H_B_prime)
      num <- num + numi 
      num2 <- num2 + numi^2
      Vi <- Vi + Vii
    }
    if(singletonFix %in% c("use mean")) {
      s <- st$s_p[[1]] # there is only one element
      # use v_a formed by comparing to overall mean
      v_a <- 1*(s - s_a_barOverall) %*% t(s - s_a_barOverall)
      numi <- diag(Vii <- H_B_prime %*% v_a %*% H_B_prime)
      num <- num + numi 
      num2 <- num2 + numi^2
      Vi <- Vi + Vii
      # strata not incrimented
    }
    if(singletonFix=="drop") {
      # do nothing, V does not change, strata not incrimented
    }
  }
  #sum variance across over all strata
  res <- Vi
  m <- length(str)
  dof <- (3.16 - 2.77/sqrt(m)) * num^2/num2
  return(list(VC=res, dof=dof, simple_dof=simple_dof))
}

getVarTaylor2 <- function(object, H_B_prime, strataVar, PSUVar,
                          singletonFix=c("drop", "use mean"), w, JohnsonRust=FALSE) {
  #find PSU per strata, and warn about dropping strata with one PSU
  singletonFix <- match.arg(singletonFix)
  stuDat <- object$stuDat[,c(strataVar, PSUVar)]
  snames <- sort(unique(stuDat[[strataVar]]))
  strata <- lapply(snames, FUN=function(x) {
    list(strat=x,
         psu=sort(unique(stuDat[stuDat[,strataVar]==x, PSUVar])))
  })
  n_psu <- lapply(strata, function(x) { length(x$psu)}) 
  simple_dof <- sum(unlist(n_psu)) - sum(length(strata))
  if (any(n_psu==1)){
    if(singletonFix == "drop") {
      warning(paste0("Of the ", length(n_psu)," strata, ", sum(n_psu<2), " strata have only one PSU. All strata with only one PSU are excluded from variance estimation. See the ", dQuote("singletonFix"), " argument for other options."))
      strata <- strata[n_psu>1] #variance estimation can only happen for Strata with more than one PSU 
    }
    if(singletonFix == "use mean") {
      warning(paste0("Of the ", length(n_psu), " strata, ", sum(n_psu<2), " strata have only one PSU. All strata with only one PSU have their value compared to the mean. See the ", dQuote("singletonFix"), " argument for more details and other options."))
    }
  }
  #only keep snames in strata 
  snames <- lapply(strata, function(x){x$strat}) 
  #split X based on psu for access later
  if("getX" %in% names(object)) {
    X <- object$getX(object$stuDat) # this object is created in the caller with X in it
  } else if("X" %in% names(object)) {
    X <- object$X[[1]]
  } else {
    X <- object$values$X
  }
  # loop through strata 
  if(!"gr" %in% names(object$values)) {
    stop("the values object must have an object named gr in it before you can calculate the Taylor series variance.")
  }
  grmat <- object$values$gr
  if(nrow(grmat) != nrow(stuDat)) {
    stop("the student data and grmat must agree on the number of rows to calculate the Taylor series variance.")
  }
  
  str <- lapply(strata, FUN=function(st) {
    # number of PSUs in this stratum
    n_a <- length(st$psu)
    # PSU index for units in this stratum and PSU
    group_index <- lapply(st$psu, FUN=function(x) {
      which(stuDat[[PSUVar]]==x & stuDat[[strataVar]]==st$strat)
    })
    s_p <- lapply(seq_along(group_index),  FUN=function(k) {
      inds <- group_index[[k]]
      # this is the gradient for the individuals in stratum st, group k
      grmatk <- grmat[inds,,drop=FALSE]
      Xk <- X[inds,,drop=FALSE]
      wk <- w[inds] # the weights for this group
      vecs <- lapply(1:ncol(Xk), function(i) {
        # for the ith coefficient, you multiply the ith column by grmatk to get the gradient for this column
        # it is summed across students to get the total gradient in the group
        apply(wk * grmatk * Xk[,i], 2, sum)
      })
      mat <- do.call(rbind, vecs)
      # as.vector turns a matrix into a vector "by column"
      # and so this returns the stacked s_p vector for the composite
      return(as.vector(mat))
    })
    st$s_p <- s_p 
    if(n_a > 1) {
      #average score across as psu in strata
      s_a_bar <- Reduce("+",s_p)/n_a
      #(s_p - s_a_bar)*(s_p - s_a_bar)' 
      v_a <- lapply(s_p, FUN=function(s){
        (s - s_a_bar) %*% t(s - s_a_bar)
      })
      st$V_a <- (n_a/(n_a-1))*Reduce("+",v_a)
    } else {
      if(singletonFix %in% "use mean") {
        # the mean should be zero, so use that
        v_a <- lapply(s_p, FUN=function(s){
          s %*% s
        })
        st$V_a <- 2*Reduce("+", v_a)
      }
      if(singletonFix %in% "drop") {
        k <- ncol(H_B_prime)
        st$V_a <- matrix(0, nrow=k, ncol=k)
      }
    }
    return(st)
  }) # end lapply(strata, FUN=function(st) {
  # when we reduce these, we will need a zero matrix
  # so grab a valid V_a and multiply it by 0
  zeroV_a <- NULL
  i <- 1
  while(is.null(zeroV_a)) {
    if(!is.null(str[[i]]$V_a)) {
      zeroV_a <- 0 * str[[i]]$V_a
    }
    i <- i + 1
  }
  # aggregate V_a
  V <- zeroV_a
  Vi <- zeroV_a
  num2 <- num <- rep(0, ncol(H_B_prime))
  numStrata <- 0
  for(i in seq_along(str)) {
    st <- str[[i]]
    npsu <- length(st$psu)
    if(npsu > 1) {
      numStrata <- numStrata + 1
    }
    if(npsu > 1 || singletonFix != "drop") {
      numi <- diag(Vii <- solve(H_B_prime) %*% st$V_a %*% solve(H_B_prime))
    }
    num <- num + numi 
    num2 <- num2 + numi^2
    Vi <- Vi + Vii
  }
  #sum variance across over all strata
  m <- length(str)
  dof <- num^2/num2
  if(JohnsonRust) {
    dof <- (3.16 - 2.77/sqrt(m)) * dof
  }
  return(list(VC=Vi, dof=dof, simple_dof=simple_dof))
}

# from EdSurvey
# author: Paul Bailey
pasteItems <- function(vector, final="and") {
  # no need to do anything if there is one or fewer elements
  if(length(vector) <= 1) {
    return(vector)
  }
  if(length(vector) == 2) {
    return(paste0(vector[1], " ", final, " ", vector[2]))
  }
  v <- vector[-length(vector)]
  f <- vector[length(vector)]
  return(paste0(paste(v, collapse=", "), ", ", final, " ", f))
}

#' @importFrom stats vcov
#' @method vcov mmlMeans
#' @export
vcov.mmlMeans <- function(object, ...){
  vcov(summary(object, ...))
}

#' @method vcov summary.mmlMeans
#' @export
vcov.summary.mmlMeans <- function(object, ...){
  object$VC
}

#' @method vcov summary.mmlCompositeMeans
#' @export
vcov.summary.mmlCompositeMeans <- function(object, ...){
  object$VC
}

#' @method vcov mmlCompositeMeans
#' @export
vcov.mmlCompositeMeans <- function(object, ...){
  vcov(summary(object, ...))
}

#' @method subset scoredTest
#' @export
subset.scoredTest <- function(x, subset, select, drop=FALSE, ...) {
  r <- eval(substitute(subset),x$stuDat)
  x$stuDat <- x$stuDat[r, , drop=FALSE]
  idVar <- x$idVar
  uid <- x$stuDat[ , idVar]
  x$stuItems <- x$stuItems[x$stuItems[ , idVar] %in% uid, ]
  x$likelihood_list <- lapply(x$likelihood_list, function(z) {
    # just those IDs in the new subset
    z[z[ , idVar] %in% uid, ]
  })
  return(x)
}

#' @method subset reducedInformationScoredTest
#' @export
subset.reducedInformationScoredTest <- function(x, subset, select, drop=FALSE, ...) {
  xi <- x$scored_test
  r <- eval(substitute(subset),xi$stuDat)
  xi$stuDat <- xi$stuDat[r, , drop=FALSE]
  idVar <- xi$idVar
  uid <- xi$stuDat[ , idVar]
  xi$stuItems <- xi$stuItems[x$istuItems[ , idVar] %in% uid, ]
  xi$likelihood_list <- lapply(xi$likelihood_list, function(z) {
    # just those IDs in the new subset
    z[z[ , idVar] %in% uid, ]
  })
  x$scored_test <- xi
  return(x)
}

#' @method coef mmlCompositeMeans
#' @export
coef.mmlCompositeMeans <- function(object, ...) {
  M <- nrow(object$coefficients)
  k <- ncol(object$coefficients)
  rawCoef <- object$coefficients
  rawCoef <- rawCoef[order(rownames(rawCoef)),]
  td <- object$testScale
  td <- td[order(td$subtest),]
  tdWeight <- c()
  tdScale <- c()
  scaledCoef <- rawCoef
  for(i in 1:M) {
    tdi <- td[i, ]
    scaledCoef[i, ] <- mml.remap.coef(rawCoef[i, ], tdi$location, tdi$scale)
    tdWeight <- c(tdWeight, tdi$subtestWeight)
    tdScale <- c(tdScale, tdi$scale)
  }
  # remap to weighted levels
  compositeCoef <- t(scaledCoef * tdWeight)
  compositeCoef <- rowSums(compositeCoef)
  # compute SD
  # v is the map applied to the covariance matrix
  v <- object$testScale$scale * object$testScale$subtestWeight
  # find the quadratic form of v and subscaleVC, which has the variances on the diagonal.
  names(compositeCoef) <- colnames(object$coefficients)
  return(compositeCoef)
}

#' @method coef mmlMeans
#' @export
coef.mmlMeans <- function(object, ...){
  co <- mml.remap.coef(object$coefficients, object$location, object$scale)
  co
}

getBareEnv <- function() {
  blank_parent <- parent.env(.GlobalEnv)
  default_env_names <- c("package:stats", "package:graphics", "package:grDevices", "package:utils", "package:datasets", "package:methods", "Autoloads", "package:base")
  while(!attr(blank_parent,"name") %in% default_env_names ) {
    blank_parent <- parent.env(blank_parent)
  }
  blank_env <-new.env(hash = TRUE, parent = blank_parent, size = 29L) 
}
reEnvFormula <- function(formula) {
  formula_char <- as.character(formula)
  blank_env <- getBareEnv()
  if(length(formula_char)==2) {
    return(formula(paste0(" ~ ",formula_char[2]),env=blank_env))
  } else {
    return(formula(paste0(formula_char[2]," ~ ",formula_char[3]),env=blank_env))
  }
}
