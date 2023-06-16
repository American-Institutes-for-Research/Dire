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
    X <- object$X
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
    newData <- as.data.frame(object$Xfull)
  }
  res <- data.frame(id=rownames(newData), stringsAsFactors=FALSE)
  for(i in 1:length(object$X)) {
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

#' @importFrom stats cov2cor pt
#' @method summary mmlCompositeMeans
#' @export
summary.mmlCompositeMeans <- function(object, gradientHessian=FALSE,
                                      varType=c("Taylor", "Partial Taylor"),
                                      clusterVar=NULL, jkSumMultiplier=1, # cluster
                                      repWeight=NULL, # replicate
                                      strataVar=NULL, PSUVar=NULL, singletonFix=c("drop", "use mean"),# Taylor
                                      ...){
  # EdSurvey 2.7.1 compatibility
  varType0 <- c("consistent", "robust", "cluster", "replicate", "Taylor")
  if(length(varType) == length(varType0) && all(varType == varType0)) {
    varType <- "Taylor"
  }
  varType <- match.arg(varType)
  sumCall <- match.call()
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
  M <- nrow(object$coefficients)
  k <- ncol(object$coefficients)
  rawCoef <- object$coefficients
  # the H_B_prime0 and (later) VCfull are block diagonal matrixes
  # "block" is a list where each element is every index for that block
  block <- list()
  for(i in 1:M) {
    block <- c(block, list(1:k + k*(i-1)))
  }

  # build H_B_prime
  H_B_prime0 <- matrix(0, nrow=k*M, ncol=k*M)
  for(i in 1:nrow(rawCoef)) {
    obji <- list(lnlf = object$lnlfl[[i]],
                 coefficients = rawCoef[i, ],
                 X = object$X[[i]],
                 stuDat = object$stuDat[[i]],
                 rr1 = object$rr1[[i]],
                 nodes = object$nodes,
                 weightVar = object$weightVar)
    ih <- getIHessian.mmlMeans(obji, gradientHessian=gradientHessian)
    H_B_prime0[block[[i]], block[[i]]] <- ih
  }
  # get weighted obs, if relevant
  if(is.null(object$weightVar)) {
    wo <- NA
  } else {
    wo <- object$weightedObs
  }

  if(varType=="consistent") {
    V0 <- -1 * H_B_prime0
  }
  if(varType=="robust") {
    varsL <- list()
    V0 <- 0 * H_B_prime0
    for(i in 1:M) {
      obj <- list(X = object$X[[i]],
                  stuDat = object$stuDat[[i]],
                  coefficients = rawCoef[i, ],
                  rr1 = object$rr1[[i]],
                  nodes = object$nodes)
      Vi <- getVarRobust(obj, H_B_prime0[block[[i]], block[[i]]])
      V0[block[[i]], block[[i]]] <- Vi
    }
  }
  if(varType=="Partial Taylor") {
    V0 <- 0 * H_B_prime0
    dof <- Inf
    for(i in 1:M) {
      obj <- list(X = object$X[[i]],
                  stuDat = object$stuDat[[i]],
                  coefficients = rawCoef[i, ],
                  rr1 = object$rr1[[i]],
                  nodes = object$nodes,
                  weightVar = object$weightVar)
      Vi <- getVarTaylor(object=obj, H_B_prime = H_B_prime0[block[[i]], block[[i]]],
                         strataVar = strataVar, PSUVar = PSUVar,
                         singletonFix = singletonFix,
                         returnVecs=FALSE)
      V0[block[[i]], block[[i]]] <- Vi$VC
      # somewhat conservative
      dof <- min(dof, Vi$simple_dof[1])
    }
    dof <- rep(dof, ncol(Vi$VC))
    td <- object$testScale
    VCfull <- V0
    sVC <- object$SubscaleVC
    for(i in 1:M) {
      for(j in 1:M) {
        if(i > j) {
          # block i, j
          VCfull[block[[i]], block[[j]]] <- td$scale[i] * td$scale[j] * cov2cor(sVC)[i,j] *
                                            sign(VCfull[block[[i]], block[[i]]] * VCfull[block[[j]], block[[j]]]) *
                                            sqrt(abs(VCfull[block[[i]], block[[i]]] * VCfull[block[[j]], block[[j]]]))
          # symetric matrix, so also set block j, i
          VCfull[ block[[j]],block[[i]]] <- t(VCfull[block[[i]], block[[j]]])
        }      
      }
    }
    # now scale within moddle covariances
    for(i in 1:M) {
      VCfull[block[[i]], block[[i]]] <- td$scale[i]^2 * VCfull[block[[i]], block[[i]]]
    }
  }
  if(varType=="Taylor") {
    Vi <- list()
    dof <- Inf
    for(i in 1:M) {
      obj <- list(X = object$X[[i]],
                  stuDat = object$stuDat[[i]],
                  coefficients = rawCoef[i, ],
                  rr1 = object$rr1[[i]],
                  nodes = object$nodes,
                  weightVar = object$weightVar)
      Vii <- getVarTaylor(object=obj, H_B_prime = H_B_prime0[block[[i]], block[[i]]],
                          strataVar = strataVar, PSUVar = PSUVar,
                          singletonFix = singletonFix,
                          returnVecs=TRUE)
      Vi[[i]] <- Vii$vec
      dof <- min(dof, Vii$simple_dof)
    }
    dof <- rep(dof, ncol(object$coef))
    str <- lapply(1:length(Vi[[1]]), function(strati) {
      lapply(1:length(Vi[[1]][[strati]]), function(psui) {
        res <- c()
        for(vi in 1:length(Vi)) {
          res <- c(res, Vi[[vi]][[strati]][[psui]])
        }
        return(res)
      })
    })
    # aggregate V_a
    num2 <- num <- rep(0, ncol(H_B_prime0))
    numStrata <- 0
    Vi <- 0 * H_B_prime0
    for(i in seq_along(str)) {
      st <- str[[i]]
      npsu <- length(st)
      if(npsu > 1) {
        v_a <- lapply(st, FUN=function(s){
          (s) %*% t(s)
        })
        st$V_a <- (npsu/(npsu-1))*Reduce("+", v_a)
        numStrata <- numStrata + 1
        numi <- diag(Vii <- H_B_prime0 %*% st$V_a %*% H_B_prime0)
        num <- num + numi 
        num2 <- num2 + numi^2
        Vi <- Vi + Vii
      } else {
        stop("error in getVarTaylor, please contact Dire developers.")
      }
    }
    #sum variance across over all strata
    td <- object$testScale
    V0 <- Vi
    VCfull <- V0
    sVC <- object$SubscaleVC
    # add covariance here,
    # update VCfull along diag below so we don't muddle the update of the covariances
    for(i in 1:M) {
      for(j in 1:M) {
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
  }
  if(is.null(colnames(VCfull))) {
    # renames
    colnames(H_B_prime0) <- rownames(H_B_prime0) <- colnames(VCfull) <- rownames(VCfull) <- 1:nrow(H_B_prime0)
    for(i in 1:M) {
      rownames(H_B_prime0)[1:k + k*(i-1)] <- 
        colnames(H_B_prime0)[1:k + k*(i-1)] <-
        rownames(VCfull)[1:k + k*(i-1)] <-
        colnames(VCfull)[1:k + k*(i-1)] <- paste0(colnames(object$coefficients), "_", names(object$X)[i])
    }
  }
  # make a containder for scaled weights
  scaledCoef <- rawCoef

  # scale them
  for(i in 1:M) {
    scaledCoef[i, ] <- mml.remap.coef(rawCoef[i, ], td$location[i], td$scale[i])
  }
  # remap to weighted levels
  compositeCoef <- coef(object)
  compositeCoef <- coef.mmlCompositeMeans(object)
  names(compositeCoef) <- colnames(object$coefficients)
  VCsmall <- matrix(0, nrow=k, ncol=k)
  for(i in 1:k) {
    vi <- rep(0, k*M)
    vi[i+((1:M)-1)*k] <- td$subtestWeight
    for(j in 1:k) {
      vj <- rep(0, k*M)
      vj[j+((1:M)-1)*k] <- td$subtestWeight
      VCsmall[i,j] <- t(vj) %*% VCfull %*% vi
    }
  }
  # this is the composite VC
  VC <- VCsmall
  rownames(VC) <- colnames(VC) <- colnames(object$coefficients)
  se <- sqrt(diag(VC))
 
  # this cannot be estimated
  se[names(compositeCoef) == "Population SD"] <- NA
  tval <- as.vector(compositeCoef/se)
  TAB <- data.frame(Estimate = compositeCoef,
               StdErr = se,
               t.value = tval)
  if(!missing("dof")) {
    TAB$dof <- dof
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
  class(res) <- "summary.mmlCompositeMeans"
  return(res)
}

#' @method summary mmlMeans
#' @export
summary.mmlMeans <- function(object, gradientHessian=FALSE,
                             varType=c("consistent", "robust", "cluster", "Taylor"),
                             clusterVar=NULL, jkSumMultiplier=1, # cluster
                             repWeight=NULL, # replicate
                             strataVar=NULL, PSUVar=NULL, singletonFix=c("drop", "use mean"),# Taylor
                             ...){
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
  H_B_prime <- getIHessian.mmlMeans(object, gradientHessian)
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
  # called from summary.mmlCompositeMeans
  if(inherits(object, "mmlCompositeMeans" )) {
    return(structure(list("VC" = VC,
                          "iHessian" = H_B_prime),
                     class="summary.mmlMeans"))
  }
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
                           "weightedObs" = wo))
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

getIHessian.comp <- function(object) {
  k <- length(object$X)
  Xl <- object$X
  co <- object$coefficients
  p <- length(co)/k
  Wsum <- rep(NA, k)
  for(ki in 1:k) {
    coki <- object$coefficients[1:p + p*(ki-1)]
    X <- Xl[[ki]]
    x_ind <- split(X, as.factor(1:nrow(X)), drop=FALSE)
    # dl/dBeta * dl/dBeta ' evaluated for each individual 
    stuDat <- object$stuDat[[ki]]
    stuDat$one <- 1
    weightVar <- object$weightVar
    if(is.null(object$weightVar)) {
      weightVar <- "one"
    }
    vars <- lapply(1:nrow(X), FUN=function(i) {
      stuDat[i,weightVar] * gradInd(location=coki,
                                    ii=i,
                                    X_subset=matrix(x_ind[[i]], nrow=1),
                                    weightVar="one",
                                    rr1=object$rr1[[ki]],
                                    stuDat=stuDat,
                                    nodes=object$nodes)
    })
    varsMat <- data.frame(id=stuDat[,object$idVar], do.call(rbind, vars))
    colnames(varsMat) <- c("id", paste0(ki,":",1:(ncol(varsMat)-1)))
    if(ki == 1) {
      varsMatOverall <- varsMat
    } else {
      varsMatOverall <- merge(varsMatOverall, varsMat, by="id", all=TRUE)
    }
    Wsum[ki] <- sum(stuDat[,weightVar])
  }
  varsMatOverall[is.na(varsMatOverall)] <- 0
  varsMatList <- split(varsMatOverall[,-1], varsMatOverall$id)
  varsMatOuterList <- lapply(varsMatList, FUN=function(gri) {
    gri <- as.numeric(gri)
    gri %*% t(gri)
  })
  Wsum <- mean(Wsum)
  H_B_prime <- -1 * Wsum/(Wsum-1) * solve(Reduce("+", varsMatOuterList))
  return(H_B_prime)
}

getIHessian.mmlMeans <- function(object, gradientHessian=FALSE, returnVars=FALSE) {
  if(gradientHessian) {
    X <- object$X
    x_ind <- split(X, as.factor(1:nrow(X)), drop=FALSE)
    # dl/dBeta * dl/dBeta ' evaluated for each individual 
    stuDat <- object$stuDat
    stuDat$one <- 1
    weightVar <- object$weightVar
    if(is.null(object$weightVar)) {
      weightVar <- "one"
    }
    vars <- lapply(1:nrow(X), FUN=function(i) {
      stuDat[i,weightVar] * gradgradT(location=object$coefficients,
                                      ii=i,
                                      X_subset=matrix(x_ind[[i]], nrow=1),
                                      weightVar="one",
                                      rr1=object$rr1,
                                      stuDat=stuDat,
                                      nodes=object$nodes)
    })
    Wsum <- sum(stuDat[,weightVar])
    H_B_prime <- -1 * Wsum/(Wsum-1) * solve(Reduce("+", vars))
  } else {
    # lnlf is actually the deviance function, so -1/2 maps it back to lnl
    # text is the warning that will happen if nearPD2 finds a condition number over 400 and so warns the user about an unstable Hessian matrix.
    H_B_prime <- solve(-1/2*nearPD2(getHessian(object$lnlf, object$coefficients), warn="Apparently singular Hessian matrix, standard error estimates are likely unstable."))
  }
  return(H_B_prime)
}

# estimate covariance matrix (Standard error)
# -1/2 turns a deviance into a likelihood
getVarConsistent <- function(object, H_B_prime) {
  #consistent  estimate
  return(-1 * H_B_prime )
}

getVarRobust <- function(object, H_B_prime) {
  X <- object$X
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
  X <- object$X
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
  X <- object$X
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
                         returnVecs=FALSE) {
  #find PSU per strata, and warn about dropping strata with one PSU
  singletonFix <- match.arg(singletonFix)
  stuDat <- object$stuDat
  X <- object$X
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
    x_groups <- lapply(split(X_strata, stuDat_strata[[PSUVar]]), matrix, ncol=ncol(X))
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
  if(singletonFix=="drop") {
    # when we reduce these, we will need a zero matrix
    # so grab a valid V_a and multiply it by 0
    zeroV_a <- NULL
    zeros <- NULL
    i <- 1
    while(all(is.null(zeroV_a), is.null(zeros))) {
      if(returnVecs) {
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

#' @method coef mmlCompositeMeans
#' @export
coef.mmlCompositeMeans <- function(object, ...) {
  M <- nrow(object$coefficients)
  k <- ncol(object$coefficients)
  rawCoef <- object$coefficients
  td <- object$testScale
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
