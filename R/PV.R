drawPVs <- function(x) {
  UseMethod("drawPVs")
}

#' @export
formPosteriors <- function(mu, sd, VC, id, priorEst) {
  # inverse of prior, robust to infinite sd
  dinv <- diag(1/sd^2, nrow=length(sd))
  VCi <- solve(VC)
  # posterior VC
  pVC <- solve(VCi + dinv)
  # posterior mean
  pmu <- pVC %*% (VCi %*% priorEst + dinv %*% mu)
  return(list(mu=pmu, VC=pVC))
}

#' @method drawPVs mmlMeans
#' @export
drawPVs.mmlMeans <- function(x, npv=5L, stochasticBeta=FALSE, normalApprox=TRUE, thetaScale=FALSE, newStuDat=NULL, newStuItems=NULL, returnPosterior=FALSE, construct) {
  if(!normalApprox) {
    stop("Only the normal approximation to the posterior distribution is supported.")
  }
  if(!any(inherits(x, c("mmlMeans", "summary.mmlMeans")))) {
    stop("Argument ", dQuote("x"), " must be an object of class ", dQuote("mmlMeans"), " or ", dQuote("summary.mmlMeans"))
  }
  npv <- as.integer(npv)
  if(npv <= 0) {
    stop("Must generate a positive number of plausible values.")
  }
  
  beta <- x$coefficients
  k <- ifelse(inherits(stochasticBeta, "data.frame"), ncol(stochasticBeta), length(stochasticBeta))
  if(!inherits(stochasticBeta, "logical")) {
    if(inherits(stochasticBeta, "data.frame")) {
      # not everything that has class data.frame fully implements data.frame
      stochasticBeta <- as.data.frame(stochasticBeta)
      beta <- stochasticBeta
      stochasticBeta <- TRUE 
      if(nrow(beta) != npv | ncol(beta) != k) {
        stop("The argument ", dQuote("stochasticBeta"), " must be a data frame with row column per nvp (", npv, ") and one column per coefficient (", k, ").")
      }
      if(any(is.na(beta))) {
        # beta is numeric and the correct length, but it contains NAs
        stop("The argument ", dQuote("stochasticBeta"), " must not have any NA values.")
      }
    } else {
      stop("The argument ", dQuote("stochasticBeta"), " must be a data frame or logical (TRUE or FALSE).")
    }
  }
  
  if(stochasticBeta & !any(inherits(x, c("mmlMeans", "summary.mmlMeans")))) {
    stop("The argument ", dQuote("summary"), " must be defined when using stochasticBeta, which uses the variance estimate from that summary for the distribution of beta.")
  }
  # these must both exist or not.
  if( (is.null(newStuDat) + is.null(newStuItems)) == 1) {
    stop("Both ", dQuote("newStuDat"), " and ", dQuote("newStuItems"), " must be defined to use new data.")
  }
  if(!is.null(newStuDat)) {
    # calculate new likelihood function
    stuDat <- newStuDat
    if(length(class(stuDat))>1) {
      # not all data.frames are the same, recast
      stuDat <- as.data.frame(stuDat)
    }
    idVar <- x$idVar
    formula <- x$formula
    stuItems <- newStuItems
    stuItems[[idVar]] <- as.character(stuItems[[idVar]])
    stuItems$key <- as.character(stuItems$key)
    stuDat[[idVar]] <- as.character(stuDat[[idVar]])
    stuItems <- stuItems[order(stuItems[[idVar]]), ]
    stuDat <- stuDat[order(stuDat[[idVar]]), ]
    cc <- complete.cases(stuDat[,c(all.vars(formula))])
    stuDat <- stuDat[cc, ]
    stuItems <- stuItems[stuItems[[idVar]] %in% stuDat[[idVar]], ]
    stuItems <- stuItems[ , c(idVar, "key", "score")]
    # drop rows with no test data on them, they get dropped from the likelihood function anyways
    stuItems <- stuItems[!is.na(stuItems$score),]
    # only keep items in the paramTab
    stuItems <- stuItems[stuItems$key %in% x$paramTab$ItemID,]
    stu <- stuItems[order(stuItems$key), ]
    stu <- split(stu, stu[[idVar]])
    stuIDs <- stuDat[[idVar]]
    missingIDs <- stuIDs[!stuIDs %in% names(stu)]
    nullData <- data.frame(oppID=character(),key=character(), score=integer())
    missingBlanks <- lapply(missingIDs, function(id) {
      return(nullData)
    })
    names(missingBlanks) <- missingIDs
    stu <- c(stu, missingBlanks)
    rr1 <- calcRR1(stu, length(x$nodes), x$polyModel, x$paramTab, x$nodes, x$fast)
    # now form X
    trms <- delete.response(terms(formula))
    m <- model.frame(trms, data=stuDat)
    X <- model.matrix(formula, m, contrasts.arg=x$contrasts, xlab=x$xlevels)
    stuDat$one <- 1
    lnlf <- fn.regression(X_=X, i=NULL, wv="one", rr1=rr1, nodes=x$nodes, stuDat=stuDat)
    
  } else {
    lnlf <- x$lnlf
    X <- x$X
    rr1 <- x$rr1
  }
  # if there is one PV, use just that
  if(stochasticBeta == FALSE) {
    posterior <- lnlf(beta, returnPosterior=TRUE)
    if(returnPosterior) {
      return(list(posterior=posterior, X=X, rr1=rr1))
    }
  }
  if(stochasticBeta & returnPosterior) {
    beta_pvs <- list()
    for(pv in 1:npv){
      # over each pv 
      posterior <- lnlf(unname(as.numeric(beta[pv,])), returnPosterior=TRUE)
      colnames(posterior)[2:3] <- paste0("pv", pv, "_" ,colnames(posterior)[2:3], "_",construct) # naming convention 
      beta_pvs <- c(beta_pvs, list(posterior))
    }
    # reduce list of dfs 
    beta_pvs <- Reduce(function(x, y) merge(x, y, by="id", suffixes=c("", "")), beta_pvs)
    return(list(posterior=beta_pvs, X=X, rr1=rr1))
  }
  # for each pv
  res <- data.frame(id=posterior$id, stringsAsFactors=FALSE)
  for(pvi in 1:npv) {
    if(stochasticBeta) {
      if(inherits(beta, "data.frame")) {
        betai <- unname(as.numeric(beta[pvi,]))
      } else {
        # x is a summary, so grab just column 1
        beta <- x$coefficients
        # need Sigma
        Sigma <- x$VC
#TODO: don't randomize sd 
        betai <- mvrnorm(1, mu=beta, Sigma=Sigma)
      }
      posterior <- lnlf(betai, returnPosterior=TRUE)
      # otherwise posterior already defined above, does not change because beta is fixed
    }
    if(thetaScale) {
      res[ , paste0("pv",pvi)] <- rnorm(nrow(res), mean=posterior$mu, sd=posterior$sd)
    } else {
      res[ , paste0("pv",pvi)] <- x$location + x$scale * rnorm(nrow(res), mean=posterior$mu, sd=posterior$sd)
    }
  }
  return(res)
}


#' @importFrom MASS mvrnorm
#' @method drawPVs mmlCompositeMeans
#' @export
drawPVs.mmlCompositeMeans <- function(x, npv=5L, stochasticBeta=FALSE, normalApprox=TRUE,
                                      thetaScale=FALSE, newStuDat=NULL, newStuItems=NULL, verbose=TRUE) {
  if(!normalApprox) {
    stop("Only the normal approximation to the posterior distribution is supported.")
  }
  npv <- as.integer(npv)
  if(npv <= 0) {
    stop("Must generate a positive number of plausible values.")
  }
  if(!inherits(stochasticBeta, "logical")) {
    stop("The argument ", dQuote("stochasticBeta"), " must be logical.")
  }
  if(stochasticBeta & !inherits(x, "summary.mmlCompositeMeans")) {
    stop("The argument ", dQuote("x"), " must be the result of a call to summary when using stochasticBeta, which uses the variance estimate from that summary for the distribution of beta.")
  }
  if(is.null(newStuDat) + is.null(newStuItems) == 1) {
    stop("Both ", dQuote("newStuDat"), " and ", dQuote("newStuItems"), " must be defined to use new data.")
  }
  if(is.null(newStuDat)) {
    newStuDat <- do.call(rbind, x$stuDat)
    newStuDat <- newStuDat[!duplicated(newStuDat[[x$idVar]]), ]
    rownames(newStuDat) <- newStuDat$sid # fix rownames 
  }
  if(is.null(newStuItems)) {
    newStuItems <- x$stuItems
  }
  beta <- x$coefficients
  # for each pv
  if(stochasticBeta) {
    beta <- x$coefficients[ , 1]
  }
  # estimate mean, VC parameters per respondent
  nConstructs <- length(x$X)
  posteriors <- list()
  rr1 <- list()
  # if stochasticBeta, one beta per PV
  if(stochasticBeta) {
    betaList <- do.call(c, lapply(x$resl, function(x){
      x$coefficients
    }))
    noSD <- !grepl("Population SD", names(betaList))
    # betaList <- betaList[noSD]
    SigmaFull <- x$VCfull
    SigmaFull <- SigmaFull[noSD, noSD]
    cat("NOTE: bad SigmaFull\n")
    for(i in 1:nrow(SigmaFull)) {
      for(j in 1:nrow(SigmaFull)) {
        if(i == j) {
          SigmaFull[i,j] <- 1
        } else {
          SigmaFull[i,j] <- 0
        }
      }
    }
    betaFull <- data.frame(matrix(betaList, byrow = TRUE, nrow = npv, ncol = length(betaList))) 
    colnames(betaFull) <- names(betaList)
    betaFull[noSD] <- mvrnorm(npv, mu=betaList[noSD], Sigma=SigmaFull) # only update regression coef, keep standard deviation 
    k <- sum(!row.names(x$coef) %in% "Population SD")
  }
  # posterior distribution for each construct
  for(ci in 1:nConstructs) {
    construct <- names(x$resl)[ci]
    cat("Calculating posterior distribution for construct ", construct, "(", ci, " of ", nConstructs, ")\n")
    args <- list(x=x$resl[[ci]],
                 npv = 1,
                 stochasticBeta = stochasticBeta,
                 normalApprox = normalApprox,
                 thetaScale = thetaScale,
                 newStuDat = newStuDat, # potentially from full data
                 newStuItems = newStuItems,
                 returnPosterior = TRUE,
                 construct)
    if(stochasticBeta) {
      args$npv <- npv
      args$stochasticBeta <- as.data.frame( betaFull[,grepl(names(x$resl)[ci], colnames(betaFull))] )
    }
    # send each construct args to mmlMenas 
    posteriorCi <-  do.call(drawPVs.mmlMeans, args)
    rr1 <- c(rr1, list(posteriorCi$rr1))
    if(ci == 1) {
      X <- posteriorCi$X
      posterior <- posteriorCi$posterior
#TODO: need to clean up so naming convention is consistent across stochastic cases 
      if(stochasticBeta){
        names(posterior)[names(posterior) != "id"] <- names(posterior)[names(posterior) != "id"]
      } else {
        names(posterior)[names(posterior) != "id"] <- paste0(names(posterior)[names(posterior) != "id"], ci)
      }
    } else {
      if(! length(all.equal(X, posteriorCi$X) == 1)) {
        stop("PV generation Xs broken, please contact package developers.")
      }
      posteriorCi <- posteriorCi$posterior
      if( !all(posterior$id == posteriorCi$id) ) {
        stop("PV generation sorting broken, please contact package developers.")
      }
#TODO: need to clean up so naming convention is consistent across stochastic cases 
      if(stochasticBeta){
        names(posteriorCi)[names(posteriorCi) != "id"] <- names(posteriorCi)[names(posteriorCi) != "id"]
      } else {
        names(posteriorCi)[names(posteriorCi) != "id"] <- paste0(names(posteriorCi)[names(posteriorCi) != "id"], ci)
      }
      
      
      posteriorCi$id <- NULL
      posterior <- cbind(posterior, posteriorCi)
    }
  }
  # fix names to not reference a construct
  posterior$id <- newStuDat[[x$idVar]]
  coefNames <- rownames(x$coefficients)
  cNames <- names(x$resl)
  # calculate correlation between constructs, per student
  if(stochasticBeta){
    # for each construct 
    for(ci in 1:(nConstructs-1)) {
      name1 <- cNames[ci]
      c1_stats <- betaFull[, paste(name1,coefNames, sep=".")]
      # for every other construct 
      for(cj in (ci+1):nConstructs) {
        name2 <- cNames[cj]
        c2_stats <- betaFull[, paste(name2,coefNames, sep=".")]
        cat("Calculating posterior correlation between construct ", name1, " and ", name2, "\n")
        for(pv in 1:nrow(betaFull)){
          c1_stats <- betaFull[pv,]
          c2_stats <- betaFull[pv,]
          s1 <- c1_stats[length(c1_stats)][[1]]
          s2 <- c2_stats[length(c2_stats)][[1]]
          fnC <- fnCor(Xb1 = as.numeric(X %*% unlist(unname(c1_stats[1:(k)])) ),
                       Xb2 = as.numeric(X %*% unlist(unname(c2_stats[1:(k)])) ),
                       # these need to account for stochastib beta
                       # this whole loop has to be 1 per npv index, and s1 and s2 must be updated per npv index
                       s1 = s1,
                       s2 = s2,
                       w = rep(1,nrow(X)),
                       rr1 = rr1[[ci]],
                       rr2 = rr1[[cj]],
                       nodes = x$nodes,
                       fine = TRUE,
                       fast = TRUE)
          FishRho0 <- atanh(x$SubscaleVC[ci,cj]/sqrt(x$SubscaleVC[ci,ci] * x$SubscaleVC[cj,cj]))
#TODO: how should we handle zero division in calculating posterior.
          posterior[paste0("pv", pv, "_rho_", name1, "_", name2)] <- fnC(FishRho0, returnPosterior=TRUE)
        }
      }
    }
  } else {
# TODO: I'm unsure if this non-Stachastic case for correlation between constructs is correct 
    s2 <- x$coefficients[ , ncol(beta)]
    for(ci in 1:(nConstructs-1)) {
      for(cj in (ci+1):nConstructs) {
        cat("Calculating posterior correlation between construct ", names(x$resl)[ci], " and ", names(x$resl)[cj], "\n")
        fnC <- fnCor(Xb1 = as.numeric(X %*% beta[ci, 1:(ncol(beta)-1)]),
                     Xb2 = as.numeric(X %*% beta[cj, 1:(ncol(beta)-1)]),
                     s1 = beta[ci, ncol(beta)],
                     s2 = beta[cj, ncol(beta)],
                     w = rep(1,nrow(X)),
                     rr1 = rr1[[ci]],
                     rr2 = rr1[[cj]],
                     nodes = x$nodes,
                     fine = TRUE,
                     fast = TRUE)
        FishRho0 <- atanh(x$SubscaleVC[ci,cj]/(s2[ci] * s2[cj]))
        posterior[ , paste0("rho_", ci, "_", cj)] <- fnC(FishRho0, returnPosterior=TRUE)
      }
    }
  }
  
  # get posteriors, by student, make PVs
  muStencil <- grepl("mu[\\d]", colnames(posterior), perl=TRUE)
  mat0 <- matrix(0, nrow=nConstructs, ncol=nConstructs)
  ts <- x$testScale
  
  cat("Calculating plausible values.\n")
  if(stochasticBeta){
    #helpers 
    mat0 <- matrix(0, nrow=nConstructs, ncol=nConstructs)
    pvListByID <- list()
    # loop through students 
    pvListByID <- lapply(1:nrow(posterior), function(rowi) {
      row_pv <- posterior[rowi,]
      pvs <- list()
      # loop through pvs 
      for(pv in 1:nrow(betaFull)){
        mu <- as.numeric(row_pv[paste0("pv", pv, "_mu_",cNames)])
        S <- mat0
        diag(S) <- (as.numeric(row_pv[paste0("pv", pv, "_sd_",cNames)])^2)
        for(i in 1:(nrow(S)-1)) {
          for(j in (i+1):ncol(S)) {
            S[i,j] <- S[j,i] <- as.numeric(row_pv[1, paste0("pv", pv, "_rho_", cNames[i], "_", cNames[j])]) * sqrt(S[i,i] * S[j,j])
          }
        }
        pvi <- data.frame(id=row_pv$id, pv=1:npv, mvrnorm(npv, mu=mu, Sigma=S))
        pvs <- c(pvs, list(pvi))
      }
      do.call(rbind, pvs) # combine all pvs for a students 
    })
    
  } else {
# TODO: non-Stochastic and Stochastic could probably be simplified 
    pvListByID <- lapply(1:nrow(posterior), function(rowi) {
      row <- posterior[rowi, ]
      mu <- as.numeric(row[muStencil])
      S <- mat0
      for(i in 1:nrow(S)) {
        S[i,i] <- as.numeric(row[1, paste0("sd",i)])^2
      }
      for(i in 1:(nrow(S)-1)) {
        for(j in (i+1):ncol(S)) {
          S[i,j] <- S[j,i] <- as.numeric(row[1, paste0("rho_", i, "_", j)]) * sqrt(S[i,i] * S[j,j])
        }
      }
      pvi <- data.frame(id=row$id, pv=1:npv, mvrnorm(npv, mu=mu, Sigma=S))
    })
  }

  pv <- do.call(rbind, pvListByID) # combine all student's pvs 
  colnames(pv) <- c(colnames(pv)[1:2], make.names(names(x$resl)))
  if(!thetaScale) {
    comp <- rep(0, nrow(pv))
    for(i in 3:ncol(pv)) {
      pv[,i] <- x$resl[[i-2]]$location + pv[,i] * x$resl[[i-2]]$scale
      w <- ts$subtestWeight[ts$subtest == names(x$resl)[[i-2]]]
      comp <- comp + w * pv[,i]
    }
    pv$comp__tmp_ <- comp
    colnames(pv)[colnames(pv) == "comp__tmp_"] <- ts$test[1]
  }
  pv <- reshape(data=pv, idvar=c("id"), direction="wide", v.names=colnames(pv)[-(1:2)], timevar="pv", sep="_")
  return(pv)
}

