#' Draw plausible values (PVs) from an mml fit
#' @param x a fit from a call to \code{\link{mml}}
#' @param npv integer indicating the number of plausible values to draw
#' @param pvVariableNameSuffix suffix added to new PV variables after construct name and before the plausible value ID. For example, if there is a construct \code{math} and the suffix is the default \code{_dire}, then the fourth plausible value would have a column name, \code{math_dire4}.
#' @param stochasticBeta logical when \code{TRUE} the regression coefficients will be drawn from their posterior distribution. Can also be a data frame of values (see Details).
#' @param normalApprox logical must be \code{TRUE} to use the normal approximation to the posterior distribution rather than drawing from the actual posterior distribution.
#' @param newStuDat new \code{stuDat} object, (see \code{\link{mml}}) for which plausible values will be drawn
#' @param newStuItems new \code{stuItems} object, (see \code{\link{mml}}); unlike in \code{mml} students with no items can be passed to this function
#' @param returnPosterior logical set to \code{TRUE} to change output to include two additional data frames (see Value).
#' @param construct character, changes the name of the columns in the final data frame
#' @param verbose logical set to \code{TRUE} to see the status of the processing
#' @param ... additional parameters
#' 
#' @details
#' 
#' When the argument passed to \code{stocasticBeta} is a data frame then each column is an element that will be used as a
#' regression coefficient for that index of the coefficients vector. The row index used for the nth PV will be the nth row.
#' 
#' @return
#' when \code{returnPosterior} is FALSE returns a object of class \code{DirePV} which is a list of two elements.
#' first, a data frame with a row for every row of newStuDat (or the original stuDat object)
#' \itemize{
#' \item{\code{id} the value of \code{idVar} in the model run}
#' \item{\code{[construct][pvVariableNameSuffix][L]} every other column is a plausible value of this format.
#'                                             The \code{[construct]} is the name of the construct,
#'                                             the \code{[pvVariableNameSuffix]} is the value of the \code{pvVariableNameSuffix} argument, and
#'                                             the \code{[L]} part is the plausible value index, from 1 to \code{npv}.}
#' }
#' The second argument is named \code{newpvvars} and is a list with an element for each set of construct that lists all of the variables in that construct.
#' 
#' When \code{returnPosterior} is TRUE returns \code{list} with three elements. One is named \code{posterior} and has one
#' row per \code{idvar} level in the \code{newStuDat} argument and three columns:
#' \itemize{
#' \item{\code{id} the value of \code{idVar} in the model run}
#' \item{\code{mu} the posterior mean}
#' \item{\code{sd} the posterior standard deviation}
#' }
#' the second list element is named \code{X} that is the design matrix for \code{newStuDat} (see Value for \code{\link{mml}}). The third list element is 
#' the \code{rr1} element returned from \code{\link{mml}} with one column for each individual in \code{newStuDat} (see Value in \code{\link{mml}}).
#' @example man/examples/PV.R
#' @author Paul Bailey, Sun-joo Lee, and Eric Buehler 
#' @export
drawPVs <- function(x, npv, pvVariableNameSuffix="_dire", ...) {
  UseMethod("drawPVs")
}

#' @rdname drawPVs
#' @method drawPVs summary.mmlMeans
#' @export
drawPVs.summary.mmlMeans <- function(x, npv=5L, pvVariableNameSuffix="_dire", ...) {
  drawPVs.mmlMeans(x=x,
                   npv=npv,
                   pvVariableNameSuffix=pvVariableNameSuffix,
                   ...)
}

#' @rdname drawPVs
#' @method drawPVs mmlMeans
#' @importFrom stats model.frame
#' @export
drawPVs.mmlMeans <- function(x, npv=5L,
                             pvVariableNameSuffix="_dire",
                             stochasticBeta=FALSE, normalApprox=TRUE, newStuDat=NULL,
                             newStuItems=NULL, returnPosterior=FALSE,
                             construct=NULL, ...) {
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
  if(inherits(beta, "matrix")) {
    beta <- x$latentCoef
  }
  k <- ifelse(inherits(stochasticBeta, "data.frame"), ncol(stochasticBeta), length(stochasticBeta))
  if(!inherits(stochasticBeta, "logical")) {
    if(inherits(stochasticBeta, "data.frame")) {
      # not everything that has class data.frame fully implements data.frame
      stochasticBeta <- as.data.frame(stochasticBeta)
      beta <- stochasticBeta
      stochasticBeta <- TRUE 
      if(nrow(beta) != npv | ncol(beta) != k) {
        stop("The argument ", dQuote("stochasticBeta"), " must be a data frame with row column per npv (", npv, ") and one column per coefficient (", k, ").")
      }
      if(any(is.na(beta))) {
        # beta is numeric and the correct length, but it contains NAs
        stop("The argument ", dQuote("stochasticBeta"), " must not have any NA values.")
      }
    } else {
      stop("The argument ", dQuote("stochasticBeta"), " must be a data frame or logical (TRUE or FALSE).")
    }
  }

  if(stochasticBeta & !returnPosterior) {
    if(!inherits(x, "summary.mmlMeans")) {
      stop(paste0("Argument ", dQuote("x"), " must be an object of class summary.mmlMeans when stochasticBeta is TRUE and returnPosterior is FALSE."))
    }
  } else {
    if(!inherits(x, c("summary.mmlMeans", "mmlMeans"))) {
      stop(paste0("Argument ", dQuote("x"), " must be an object of class mmlMeans or summary.mmlMeans"))
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
    # keeps rr1 and X in same order
    stu <- stu[ order(names(stu)) ]
    rr1 <- calcRR1(stu, length(x$nodes), x$polyModel, x$paramTab, x$nodes, x$fast)
    # now form X
    trms <- delete.response(terms(formula))
    m <- model.frame(trms, data=stuDat, drop.unused.levels=TRUE)
    X <- model.matrix(formula, m, contrasts.arg=x$contrasts, xlab=x$xlevels)
    # fn.regression uses rownames to name posterior
    rownames(X) <- names(stu)
    lnlf <- fn.regression(X_=X, i=NULL, wv=x$weightVar, rr1=rr1, nodes=x$nodes, stuDat=stuDat, inside=FALSE)
  } else { #   if(!is.null(newStuDat)) {
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
  if(stochasticBeta != FALSE & returnPosterior) {
    beta_pvs <- list()
    for(pv in 1:npv){
      # over each pv 
      if(inherits(beta, "data.frame")) {
        betai <- as.numeric(beta[pv, ])
      } else {
        # directly from x
        betai <- beta 
      }
      posterior <- lnlf(unname(betai), returnPosterior=TRUE)
      if(!is.null(construct)) {
        colnames(posterior)[2:3] <- paste0("pv", pv, "_" ,colnames(posterior)[2:3], "_", construct) # naming convention 
      } else {
        colnames(posterior)[2:3] <- paste0("pv", pv, "_" ,colnames(posterior)[2:3]) # naming convention 
      }
      beta_pvs <- c(beta_pvs, list(posterior))
    }
    # reduce list of dfs 
    beta_pvs <- Reduce(function(x, y) merge(x, y, by="id", suffixes=c("", "")), beta_pvs)
    return(list(posterior=beta_pvs, X=X, rr1=rr1))
  }
  # stochastibBeta is TRUE but returnPosterior is FALSE
  # no posterior needed
  # for each pv
  res <- data.frame(id=rownames(X), stringsAsFactors=FALSE)
  for(pvi in 1:npv) {
    if(stochasticBeta) {
      if(inherits(beta, "data.frame")) {
        betai <- unname(as.numeric(beta[pvi,]))
      } else {
        beta <- x$coefficients
        if(inherits(x, "summary.mmlMeans")) {
          # x is a summary, so grab just column 1
          beta <- x$latentCoef
        }
        # need Sigma
        Sigma <- x$VC
        # do not randomize SD
        Sigma[nrow(Sigma), ] <- 0
        Sigma[ , nrow(Sigma)] <- 0
        Sigma <- Sigma
        betai <- mvrnorm(1, mu=beta, Sigma=Sigma)
      }
      posterior <- lnlf(betai, returnPosterior=TRUE)
      # otherwise posterior already defined above, does not change because beta is fixed
    }
    res[ , paste0("pv", pvi)] <- x$location + x$scale * rnorm(nrow(res), mean=posterior$mu, sd=posterior$sd)
  }
  # update names
  finalNames <- paste0(construct, pvVariableNameSuffix)
  newPV <- lapply(finalNames, function(x) {
    return(list(varnames=paste0(x,1:npv)))
  })
  names(newPV) <- finalNames
  colnames(res) <- c("id", newPV[[1]]$varnames)
  return(structure(list(data=res, newpvvars=newPV), class="DirePV"))
}

#' @method drawPVs summary.mmlCompositeMeans
#' @export
drawPVs.summary.mmlCompositeMeans <- function(x, npv=5L, pvVariableNameSuffix="_dire", ...) {
  drawPVs.mmlCompositeMeans(x=x,
                            npv=npv,
                            pvVariableNameSuffix=pvVariableNameSuffix,
                            ...)
}

getMMLCBeta <- function(x) {
  if(inherits(x, "summary.mmlCompositeMeans")) {
    return(x$rawCoef)
  }
  return(x$coefficients)
}


#' @rdname drawPVs
#' @importFrom MASS mvrnorm
#' @method drawPVs mmlCompositeMeans
#' @export
drawPVs.mmlCompositeMeans <- function(x, npv=5L, pvVariableNameSuffix="_dire", 
                                      stochasticBeta=FALSE, normalApprox=TRUE,
                                      newStuDat=NULL, newStuItems=NULL, verbose=TRUE, ...) {
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
    rownames(newStuDat) <- newStuDat[[x$idVar]] # fix rownames 
  }
  # this is key because the posterior and newStuDat must be in the same order
  newStuDat <- newStuDat[order(as.character(newStuDat[[x$idVar]])),]
  if(is.null(newStuItems)) {
    newStuItems <- x$stuItems
  }
  beta <- getMMLCBeta(x)
  # estimate mean, VC parameters per respondent
  nConstructs <- length(x$X)
  posteriors <- list()
  rr1 <- list()
  # if stochasticBeta, one beta per PV
  if(stochasticBeta) {
    betaList <- do.call(c, lapply(x$resl, function(x){
      getMMLCBeta(x)
    }))
    noSD <- !grepl("Population SD", names(betaList))
    # betaList <- betaList[noSD]
    SigmaFull <- x$VCfull
    SigmaFull <- SigmaFull[noSD, noSD]
    if(is.null(colnames(SigmaFull))) {
      stop("no column names on x$VCfull")
    }
    # SigmaFull is in score space, but we want it in theta space
    ts <- x$testScale
    for(tsi in 1:ncol(ts)) {
      tsiCols <- grepl(paste0("_", ts[tsi, "subtest"],"$"), colnames(SigmaFull), perl=TRUE) 
      SigmaFull[tsiCols, ] <- SigmaFull[tsiCols,] / ts$scale[tsi]
      SigmaFull[ , tsiCols] <- SigmaFull[ , tsiCols] / ts$scale[tsi]
    }
    betaFull <- data.frame(matrix(betaList, byrow = TRUE, nrow = npv, ncol = length(betaList))) 
    colnames(betaFull) <- names(betaList)
    betaFull[noSD] <- mvrnorm(npv, mu=betaList[noSD], Sigma=SigmaFull) # only update regression coef, keep standard deviation 
    k <- sum(!row.names(x$coef) %in% "Population SD")
  }
  # posterior distribution for each construct
  for(ci in 1:nConstructs) {
    construct <- names(x$resl)[ci]
    if(verbose) {
      message("Calculating posterior distribution for construct ", construct, " (", ci, " of ", nConstructs, ")")
    }
    args <- list(x=x$resl[[ci]],
                 npv = 1,
                 pvVariableNameSuffix="",
                 construct=construct,
                 stochasticBeta = stochasticBeta,
                 normalApprox = normalApprox,
                 newStuDat = newStuDat, # potentially from full data
                 newStuItems = newStuItems,
                 returnPosterior = TRUE)
    if(stochasticBeta) {
      args$npv <- npv
      args$stochasticBeta <- as.data.frame( betaFull[,grepl(paste0("^",names(x$resl)[ci],"[.]"), colnames(betaFull))] )
    }
    # send each construct args to mmlMenas 
    posteriorCi <-  do.call(drawPVs.mmlMeans, args)
    rr1 <- c(rr1, list(posteriorCi$rr1))
    if(ci == 1) {
      X <- posteriorCi$X
      posterior <- posteriorCi$posterior
      if(nrow(posterior) != nrow(X) || any(posterior$id != rownames(X))) {
        stop("error ordering variables in X, please contact package developers.")
      }
      if(stochasticBeta){
        names(posterior)[names(posterior) != "id"] <- names(posterior)[names(posterior) != "id"]
      } else {
        names(posterior)[names(posterior) != "id"] <- paste0(names(posterior)[names(posterior) != "id"], ci)
      }
    } else {
      if(! length(all.equal(X, posteriorCi$X)) == 1) {
        stop("PV generation Xs not formated correctly, please contact package developers.")
      }
      posteriorCi <- posteriorCi$posterior
      if( !all(posterior$id == posteriorCi$id) ) {
        stop("PV generation sorting broken, please contact package developers.")
      }
      if(stochasticBeta){
        names(posteriorCi)[names(posteriorCi) != "id"] <- names(posteriorCi)[names(posteriorCi) != "id"]
      } else {
        names(posteriorCi)[names(posteriorCi) != "id"] <- paste0(names(posteriorCi)[names(posteriorCi) != "id"], ci)
      }
      posterior <- merge(posterior, posteriorCi, by="id", all=TRUE)
    }
  }
  # fix names to not reference a construct
  if(nrow(posterior) != nrow(newStuDat) || any(posterior$id != newStuDat[[x$idVar]])) {
    stop("error ordering variables in newStuDat, please contact package developers")
  }

  coefNames <- rownames(x$coefficients)
  cNames <- names(x$resl)
  # calculate correlation between constructs, per student
  if(stochasticBeta){
    # for each construct 
    sbind <- 1
    sbtot <- (nConstructs - 1) * nConstructs/2
    for(ci in 1:(nConstructs-1)) {
      name1 <- cNames[ci]
      c1_stats0 <- betaFull[, paste(name1, coefNames, sep=".")]
      # for every other construct 
      for(cj in (ci+1):nConstructs) {
        name2 <- cNames[cj]
        if(verbose) {
          message(paste0("Calculating posterior correlation between construct ", name1, " and ", name2, " (", sbind, " of ", sbtot, ")"))
          sbind <- sbind + 1
        }
        c2_stats0 <- betaFull[, paste(name2, coefNames, sep=".")]
        for(pv in 1:nrow(betaFull)){
          c1_stats <- c1_stats0[pv, ]
          c2_stats <- c2_stats0[pv, ]
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
                       fast = x$resl[[1]]$fast)
          if( x$SubscaleVC[ci,ci] * x$SubscaleVC[cj,cj] <= 0) {
            stop("subscaleVC without variance, cannot process")
          }
          FishRho0 <- atanh(x$SubscaleVC[ci,cj]/sqrt(x$SubscaleVC[ci,ci] * x$SubscaleVC[cj,cj]))
          posterior[paste0("pv", pv, "_rho_", name1, "_", name2)] <- fnC(FishRho0, returnPosterior=TRUE)
        }
      }
    }
  } else {
    s2 <- getMMLCBeta(x)
    # got population SD, final column
    s2 <- s2[ , colnames(s2) == "Population SD"]
    sbind <- 1
    sbtot <- (nConstructs - 1) * nConstructs/2
    for(ci in 1:(nConstructs-1)) {
      name1 <- cNames[ci]
      for(cj in (ci+1):nConstructs) {
        name2 <- cNames[cj]
        if(verbose) {
          message(paste0("Calculating posterior correlation between construct ", name1, " and ", name2, " (", sbind, " of ", sbtot, ")"))
          sbind <- sbind + 1
        }
        fnC <- fnCor(Xb1 = as.numeric(X %*% beta[ci, 1:(ncol(beta)-1)]),
                     Xb2 = as.numeric(X %*% beta[cj, 1:(ncol(beta)-1)]),
                     s1 = beta[ci, ncol(beta)],
                     s2 = beta[cj, ncol(beta)],
                     w = rep(1, nrow(X)),
                     rr1 = rr1[[ci]],
                     rr2 = rr1[[cj]],
                     nodes = x$nodes,
                     fine = TRUE,
                     fast = x$resl[[1]]$fast)
        FishRho0 <- atanh(x$SubscaleVC[ci,cj]/(s2[ci] * s2[cj]))
        posterior[ , paste0("rho_", ci, "_", cj)] <- fnC(FishRho0, returnPosterior=TRUE)
      }
    }
  }
  
  # get posteriors, by student, make PVs
  muStencil <- grepl("mu[\\d]", colnames(posterior), perl=TRUE)
  mat0 <- matrix(0, nrow=nConstructs, ncol=nConstructs)
  ts <- x$testScale
  
  if(verbose) {
    message("Generating plausible values.")
  }
  
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
        S <- nearPD2(S)
        pvi <- data.frame(id=row_pv$id, pv=pv, t(mvrnorm(1, mu=mu, Sigma=S)))
        pvs <- c(pvs, list(pvi))
      }
      do.call(rbind, pvs) # combine all pvs for a students 
    })
    
  } else {
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
      S <- nearPD2(S)
      pvi <- data.frame(id=row$id, pv=1:npv, mvrnorm(npv, mu=mu, Sigma=S))
    })
  }
  pv <- do.call(rbind, pvListByID) # combine all student's pvs 

  colnames(pv) <- c(colnames(pv)[1:2], make.names(names(x$resl))) # columns ID, PV will be present regardless of number of pvs 
  comp <- rep(0, nrow(pv))
  for(i in 3:ncol(pv)) {
    pv[,i] <- x$resl[[i-2]]$location + pv[,i] * x$resl[[i-2]]$scale
    w <- ts$subtestWeight[ts$subtest == names(x$resl)[[i-2]]]
    comp <- comp + w * pv[,i]
  }
  pv$comp__tmp_ <- comp
  colnames(pv)[colnames(pv) == "comp__tmp_"] <- ts$test[1]
  constructNames <- colnames(pv)[-(1:2)]
  finalNames <- paste0(constructNames, pvVariableNameSuffix)
  newPV <- lapply(finalNames, function(x) {
    return(list(varnames=paste0(x,1:npv)))
  })
  names(newPV) <- finalNames
  pv <- reshape(data=pv, idvar=c("id"), direction="wide", v.names=colnames(pv)[-(1:2)], timevar="pv", sep=pvVariableNameSuffix)
  return(structure(list(data=pv, newpvvars=newPV), class="DirePV"))
}

#' @importFrom Matrix nearPD
nearPD2 <- function(X, tol=400, warn="") {
  eig <- eigen(X)
  if(min(eig$values) <= 0 || max(eig$values)/min(eig$values) >= 1/((.Machine$double.eps)^0.25)) {
    if(nchar(warn) > 0) {
      warning(warn)
    }
    X <- nearPD(X,  posd.tol=tol*sqrt(.Machine$double.eps))$mat
  }
  return(X)
}