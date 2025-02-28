mmlu <- function(formula,
                 scored_data,
                 verbose=0,
                 optimizer=c("EM", "QN"),
                 fixedSD=NULL) {
  formula <- reEnvFormula(formula)
  call <- match.call()
  optimizer <- match.arg(optimizer)
  # check for parrallel if multiCore True 
  if(length(formula) <= 2) {
    stop("must specify an outcome")
  }
  outcome <- as.character(formula[[2]])
  is_ri <- inherits(scored_data,"reducedInformationScoredTest")
  if(is_ri) {
    scored_data0 <- scored_data
    scored_data <- scored_data0$scored_test
  }
  if(outcome %in% scored_data$paramTab$subtest) {
    paramTab <- scored_data$paramTab[scored_data$paramTab$subtest %in% outcome,]
  } else {
    paramTab <- scored_data$paramTab[scored_data$paramTab$test %in% outcome,]
  }
  
  # this also checks for valid testScale
  ls <- getLocationScale(scored_data$testScale, outcome)
  score_lik <- getLikelihood(scored_data, outcome)
  formula[[2]] <- NULL
  if(is_ri) {
    formula <- scored_data0$formula0
    formula[[2]] <- NULL
  }
  stuDat <- scored_data$stuDat
  # done with the left hand side of the formula, if it exists, drop it
  # subset stuDat to valid data
  idVar <- scored_data$idVar
  cc <- complete.cases(stuDat[ , c(all.vars(formula), scored_data$weightVar, scored_data$strataVar, scored_data$PSUVar, idVar)])
  stuDat <- stuDat[cc, ]
  # subset the scores to 
  score_lik <- score_lik[rownames(score_lik) %in% stuDat[ , idVar], ]
  # go the other way too
  stuDat <- stuDat[stuDat[ , idVar] %in% rownames(score_lik), ]
  if(nrow(stuDat) == 0) {
    stop(paste0("no complete cases in ", dQuote("stuDat"), "."))
  }
  # order
  score_lik <- score_lik[order(rownames(score_lik), method="radix"), ]
  stuDat <- stuDat[order(stuDat[,idVar], method="radix"), ]
  rr1 <- t(score_lik)
  
  # now form X
  trms <- delete.response(terms(formula))
  m <- model.frame(trms, data=stuDat, drop.unused.levels = TRUE)
  X <- model.matrix(formula, m)
  Xnames <- colnames(X)
  #for prediction
  contrasts <- attributes(X)$contrasts
  xlevels <- .getXlevels(trms, m)
  V <- diag(ncol(X))
  skipCols <- c()
  if(is_ri) {
    # notice this happens after we captured how X used contrasts and levels
    V <- scored_data0$V
    skipCols <- scored_data0$skip_columns
    if(length(skipCols)>0) {
      X <- cbind(X[,skipCols], X[ , -skipCols, drop=FALSE] %*% scored_data0$V)
    }
    nPC <- 1:ncol(scored_data0$V)
    Xnames <- c(Xnames[skipCols], paste0("PC",nPC))
  }
  
  tryCatch(eig <- eigen(crossprod(X)), error=function(e) {
    print(crossprod(X))
    cat("eigen(crossprod(X)) in mml.R about line 74\n")
    print(e)
  })
  ev <- eig$values/eig$values[1]
  if(verbose >= 1) {
    message(paste0("design matrix condition number = ", round(eig$values[1] / eig$values[length(eig$values)],4) ))
  }
  Xf <- rep(TRUE, ncol(X)) # default, keep all X columns
  if(any(ev <= 100 * .Machine$double.eps)) {
    # remove columns that qr would remove
    w <- stuDat[ , scored_data$weightVar]
    cf <- qr.coef(qr(sqrt(w)*X), rep(0, nrow(X)))
    Xf <- !is.na(cf)
  }
  
  X <- X[ , Xf, drop=FALSE]
  Xnames <- Xnames[Xf]
  rownames(X) <- stuDat[[idVar]]
  colnames(X) <- Xnames
  makeGetX <- function(formula, trms, contrasts, xlevels, Xf, idVar, V, skipCols, Xnames) {
    function(stuDat) {
      m <- model.frame(trms, data=stuDat, drop.unused.levels = TRUE)
      X <- model.matrix(formula, stuDat, xlev=xlevels, contrasts.arg=contrasts)
      X <- X[, Xf, drop=FALSE]
      if(length(skipCols)>0) {
        X <- cbind(X[ , skipCols, drop=FALSE], X[ , -skipCols, drop=FALSE] %*% V)
      }
      colnames(X) <- Xnames
      rownames(X) <- stuDat[[idVar]]
      return(X)
    }
  }
  getX <- makeGetX(formula=formula, trms=trms, contrasts=contrasts, xlevels=xlevels, Xf=Xf, idVar=idVar, V=V, skipCols=skipCols, Xnames=Xnames)
  
  K <- ncol(X) # number of fixed parameters to estimate   
  nms <- c(Xnames, 'Population SD')  
  if(!is.null(fixedSD)) {
    # no SD needed
    startVal <- rep(0, K)
    # but may need to may fixedSD back
    if(fixedSD < 1) {
      fixedSD <- 2*log(fixedSD) + 1
      # undoes this, which is done after fitting
      #opt$par[length(opt$par)] <- sqrt(ifelse(opt$par[length(opt$par)] < 1, exp(opt$par[length(opt$par)] - 1), opt$par[length(opt$par)]^2))
    }
  } else {
    # default SD is 1
    startVal <- c(rep(0, K), 1)
  }
  
  if(nrow(X) != nrow(stuDat)) {
    stop("Missing values not allowed in independent variables.")
  }
  ### These are the exported functions passed to the dopar function
  nodes <- scored_data$nodes
  if(!all.equal(colnames(rr1), stuDat[[idVar]])) {
    stop("Sorting error in mml.")
  }
  fn2 <- fn.regression(X_=X, i=NULL, wv=scored_data$weightVar, rr1=rr1, nodes=nodes, stuDat=stuDat)
  func_lnl <- fn.regression_lnl(X_=X, w=stuDat[[scored_data$weightVar]], rr1=rr1, nodes=nodes)
  func_gr_lnl <- fn.regression_gr_lnl(X_=X, w=stuDat[[scored_data$weightVar]], rr1=rr1, nodes=nodes)
  func_hess_lnl <- fn.regression_hess_lnl(X_=X, w=stuDat[[scored_data$weightVar]], rr1=rr1, nodes=nodes)
  
  if("EM" == optimizer) {
    if(verbose > 0) {
      cat("Finding MLE using EM\n")
    }
    opt <- updateGamma(rr1 = rr1, nodes = nodes, X = X, w=stuDat[[scored_data$weightVar]], verbose=verbose-1)
    # func_lnl must be used here (not fn2) because s2 is not mapped like it is in fn2
    opt$value <- func_lnl(opt$par)
    assign("insd", FALSE, envir=environment(fun=fn2))
  } else {
    if(verbose > 0) {
      cat("Finding MLE using QN\n")
    }
    opt <- robustOptim(fn2, startVal, verbose=verbose, X=X, fixedSD=fixedSD) 
  }
  posteriorEsts <- fn2(opt$par, returnPosterior=TRUE)
  names(opt$par) <- c(Xnames, "s")
  
  weightVar <- scored_data$weightVar
  if(is.null(weightVar)) {
    obs <- sum(apply(rr1, 2, sum)>0)
    weightedObs <- obs
  } else {
    obs <- sum(stuDat[ , weightVar] > 0 & apply(rr1, 2, sum) > 0) # number with positive weight, at least one response
    weightedObs <- sum(stuDat[apply(rr1, 2, sum) > 0, weightVar]) # sum of weights for those with positive weight, at least one response
  }
  coefficients <- opt$par
  names(coefficients) <- nms
  assign("insd", FALSE, envir=environment(fun=fn2))
  pred <- as.vector(X %*% coefficients[-length(coefficients)])
  names(pred) <- rownames(X)
  oneD_E_sd <- getOneD_E_sd(pred, rr1, coefficients[length(coefficients)], nodes)
  n <- nrow(stuDat)
  values <- list(n = n, 
                 k = 1,
                 Q = length(nodes),
                 idVar=idVar,
                 id_vec = stuDat[,idVar],
                 nodes=nodes,
                 w=stuDat[,weightVar],
                 Sigma = matrix(1),
                 E = matrix(oneD_E_sd$E, nrow=n, ncol=1),
                 Xb = matrix(pred, ncol=1),
                 X = X,
                 cov = matrix(oneD_E_sd$V, nrow=n, ncol=1),
                 rr1_list = rr1,
                 rr1_splines= NULL,
                 sigma_list=lapply(oneD_E_sd$V, function(x) {sqrt(x)})) # not used
  res <- structure(list(call = call,
                        coefficients = coefficients,
                        LogLik = -1/2*func_lnl(opt$par),
                        getX = getX,
                        Convergence = opt$convergence,
                        location = ls$location,
                        scale = ls$scale,
                        lnlf= fn2,
                        funcs = list(lnl=func_lnl, gr=func_gr_lnl, hess=func_hess_lnl),
                        rr1= rr1,
                        stuDat = stuDat,
                        likelihood_stus = score_lik,
                        weightVar = weightVar,
                        nodes = scored_data$nodes,
                        iterations = opt$iter,
                        obs = obs,
                        weightedObs = weightedObs,
                        strataVar = scored_data$strataVar,
                        contrasts = contrasts,
                        xlevels = xlevels,
                        PSUVar = scored_data$PSUVar,
                        formula = formula,
                        polyModel = scored_data$polyModel,
                        paramTab = paramTab,
                        idVar = idVar,
                        posteriorEsts = posteriorEsts,
                        itemScorePoints = scored_data$responseRanges,
                        outcome=outcome,
                        fast = scored_data$fast,
                        V = V,
                        testScale=scored_data$testScale,
                        pred=pred,
                        values=values),
                   class = "mmlMeans")
  return(res)
}


mmlcomp <- function(mmlList, multiCore=FALSE, optimizer=c("MLE"), verbose=0, vc=NULL) {
  
  is_mmlMeans <- vapply(X=mmlList, FUN=function(x) { inherits(x, "mmlMeans")}, FUN.VALUE=logical(1))
  if(!all(is_mmlMeans)) {
    stop(paste0("mmlList must be a list of ", sQuote("mmlMeans") , "objects"))
  }
  call <- match.call()
  optimizer <- match.arg(optimizer)
  resL <- mmlList
  # nodes are all the same
  nodes <- resL[[1]]$nodes
  k <- length(resL[[1]]$coefficients)
  subtests <- unlist(lapply(resL, function(x) {
    x$outcome
  }))

  idVar <- resL[[1]]$idVar
  ids <- resL[[1]]$stuDat[[idVar]]
  stuDatl <- list()
  for(sti in 1:length(resL)) {
    # fix posteriors
    nidv <- colnames(resL[[sti]]$posteriorEsts) != "id" # non id variable names
    colnames(resL[[sti]]$posteriorEsts)[nidv] <- make.names(paste0(colnames(resL[[sti]]$posteriorEsts)[nidv], "_", subtests[sti]))
    if(sti == 1) {
      posteriorEsts <- resL[[sti]]$posteriorEsts
    } else {
      posteriorEsts <- merge(posteriorEsts, resL[[sti]]$posteriorEsts, by="id", all=TRUE)
    }
    # remove, condensation complete
    resL[[sti]]$posteriorEsts <- NULL
    # get Xb
    co <- resL[[sti]]$coefficients
    # remove standard deviation
    Xbi <- as.vector(resL[[sti]]$pred)
    Xbdf <- data.frame(id=resL[[sti]]$stuDat[,idVar], xb=Xbi, stringsAsFactors=FALSE)
    colnames(Xbdf)[2] <- paste0("Xb", sti)
    if(sti > 1) {
      Xb <- merge(Xb, Xbdf, by="id", all=TRUE)
    } else {
      Xb <- Xbdf
    }
    resL[[sti]]$stdev <- resL[[sti]]$coefficients[length(resL[[sti]]$coefficients)]
    stuDatl <- c(stuDatl, list(resL[[sti]]$stuDat))
    resL[[sti]]$iter <- resL[[sti]]$iterations
    resL[[sti]]$iterations <- NULL
    resL[[sti]]$coef <- resL[[sti]]$coefficients
    resL[[sti]]$coefficients <- NULL
    resL[[sti]]$call <- NULL
    resL[[sti]]$itemScorePoints <- NULL
    resL[[sti]]$wobs <- resL[[sti]]$weightedObs
    resL[[sti]]$weightedObs <- NULL
    if(sti == 1) {
      weightVar <- resL[[sti]]$weightVar
      idVar <- resL[[sti]]$idVar
      testScale <- resL[[sti]]$testScale
      testScale <- testScale[testScale$subtest %in% subtests, ]
      reorder <- function(var, ref) {
        res <- 1:length(ref)
        for(i in 1:length(res)) {
          res[i] <- which(var == ref[i])
        }
        return(res)
      }
      testScale <- testScale[reorder(testScale$subtest, subtests), ]
      if(verbose > 2) {
        cat("testScale after filtering to relevant tests:\n")
        print(testScale)
      }
      strataVar <- resL[[sti]]$strataVar
      PSUVar <- resL[[sti]]$PSUVar
    } else {
      if(resL[[sti]]$weightVar != weightVar) {
        stop("You must use same weight variable to make a composite. Every observation needs a unique weight for calculating correlation.")
      }
      if(resL[[sti]]$idVar != idVar) {
        stop("The idVar argument must be the same across univariate calls. Try refitting with idVar set to the same variable name in every fit.")
      }
      if(!is.null(resL[[sti]]$strataVar) && resL[[sti]]$strataVar != strataVar) {
        stop("The strataVar argument must be the same across univariate calls. Try refitting with strataVar set to the same variable name in every fit.")
      }
      if(!is.null(resL[[sti]]$PSUVar) && resL[[sti]]$PSUVar != PSUVar) {
        stop("The PSUVar argument must be the same across univariate calls. Try refitting with PSUVar set to the same variable name in every fit.")
      }
    }
    resL[[sti]]$strataVar <- NULL
    resL[[sti]]$PSUVar <- NULL
    resL[[sti]]$scored_data <- NULL
  }
  # store full res list

  for(sti in 1:length(resL)) {
    # drop nodes
    resL[[sti]]$nodes <- NULL
  }
  # merge stuDat
  stuDat <- do.call(rbind, stuDatl)
  stuDat <- stuDat[!duplicated(stuDat[,idVar]),]
  stuDat <- stuDat[order(stuDat[,idVar]),]

  subsets <- list()
  for(i in seq_along(resL)) {
    subsets <- c(subsets, list(stuDat[[idVar]] %in% resL[[i]]$stuDat[[idVar]]))
  }
  stuDatSubsets <- do.call(cbind, subsets)
  colnames(stuDatSubsets) <- subtests
  for(i in seq_along(resL)) {
    resL[[i]]$stuDat <- NULL
  }
  # exclude is for variables we don't want to set globally,
  # but need resL to retain for post hoc

  # assign variables from the many calls back onto the current frame (the env for this call to mml)
  setRes(resL, exclude=c("formula","polyModel","paramTab","fast","scale",
                         "location", "weightVar", "idVar"))
  # merge together posterior estimates
  # nodes are the same for every run
  if(is.null(weightVar)) {
    stuDat$one <- 1
    weightVar <- "one"
  }
  testScale <- filterTestScale(testScale, subtests)
  if(is.null(vc)) {
    vc <- matrix(0, nrow=length(subtests), ncol=length(subtests))
    diag(vc) <- unlist(lapply(resL, function(x) { x$stdev }))^2
    colnames(vc) <- rownames(vc) <- subtests
  }
  # this object holds all the results. For now much of it is blank.
  # the rest of this function fills in components of it
  k <- length(resL)
  n <- nrow(stuDat)
  values <- list(n = n, 
                 k = k,
                 Q = length(nodes),
                 idVar=idVar,
                 id_vec = stuDatl[[1]][,idVar],
                 nodes=nodes,
                 w=stuDatl[[1]][,weightVar],
                 Sigma = vc,
                 E = matrix(NA, nrow=n, ncol=k),
                 Xb = Xb[,paste0("Xb",1:k)],
                 X = resL[[1]]$getX(stuDat),
                 cov = matrix(NA, nrow=n, ncol=sum(1:k)),
                 rr1_list = lapply(resL, function(x) { x$rr1 }) ,
                 rr1_splines= NULL) # not used

  if(multiCore) {
    stop("no multicore yet")
  } else { # end: if(multiCore)
    for(i in 1:(length(subtests)-1)) {
      for(j in (i+1):length(subtests)) {
        if(verbose > 0) {
          cat("optimizing element ",i,",",j,"\n")
        }
        values <- mml_MLE_2D(values, i, j, verbose=verbose-2)
      }
    }
  } # end else for if(multicore)

  # updates E_i, E_j, vari_vec, varj_vec, cov_vec
  values <- e_step_nD_by2D(values) 

  # add names to everything
  names(iter) <- names(obs) <- names(wobs) <- subtests
  coefM <- matrix(coef, nrow=length(subtests), byrow=TRUE)
  rownames(coefM) <- subtests
  colnames(coefM) <- names(coef)[1:ncol(coefM)]
  res <- structure(list(call = call,
                        coefficients = coefM,
                        ids = ids,
                        Convergence = Convergence,
                        stuDat = stuDat,
                        weightVar = weightVar,
                        nodes = nodes,
                        iterations = iter,
                        obs = obs,
                        testScale = testScale,
                        weightedObs = wobs,
                        idVar = idVar,
                        resl = resL,
                        strataVar = strataVar,
                        PSUVar = PSUVar,
                        posteriorEsts = posteriorEsts,
                        SubscaleVC = values$Sigma,
                        stuDatSubsets=stuDatSubsets,
                        values=values
                        ),
                   class = "mmlCompositeMeans")
  return(res)
}


filterTestScale <- function(ts, subtests) {
  ts <- ts[ts$subtest %in% subtests, ]
  ts$test_subtest <- paste(ts$test, ts$subtest, ":")
  ts <- ts[!duplicated(ts$test_subtest),]
  ts$test_subtest <- NULL
  if(nrow(ts) != length(subtests)) {
    stop("testScale should have exactly one row for each subscale in this composite.")
  }
  return(ts)
}

getLikelihood <- function(scored_data, outcome) {
  idVar <- scored_data$idVar
  paramTab <- scored_data$paramTab
  if(!outcome %in% c(paramTab$test, paramTab$subtest)) {
    stop("the outcome variable must be a test or subtest on the test data.")
  }
  if(outcome %in% paramTab$test) {
    subtests <- unique(paramTab$subtest[paramTab$test == outcome])
    for(i in seq_along(subtests)) {
      st <- subtests[i]
      liki <- scored_data$likelihood_list[[st]]
      if(i == 1) {
        lik <- liki
      } else {
        lik <- rbind(lik, liki)
      }
    }
    # condense this by person ID, summing them
    lik <- split(lik, lik[[idVar]])
    lik <- lapply(lik, function(x) {
      exp(apply(log(x[,-1]), 2, sum))
    })
    lik <- do.call(rbind, lik)
  } else {
    lik <- scored_data$likelihood_list[[outcome]]
    rownames(lik) <- lik[,1]
    lik <- lik[,-1]
  }
  return(lik)
}

getLocationScale <- function (testScale, outcome) {
  if(outcome %in% testScale$test) {
    ts <- testScale[testScale$test %in% outcome, ]
    if(nrow(ts) > 1) {
      ts <- ts[is.na(ts$subtest) | ts$subtest %in% outcome, ]
      if(nrow(ts) > 1) {
        print(ts)
        stop("test scale unclear for outcome variable ", dQuote(outcome), " the testScale (passed as an argument to scoreTest) must have a single row that has missing subtest. The locaiton and scale from that line will then be used to scale this outcome.")
      }
    }
  } else {
    if(outcome %in% testScale$subtest) {
      ts <- testScale[testScale$subtest %in% outcome, ]
      if(nrow(ts) > 1) {
        print(ts)
        stop("test scale unclear for outcome variable ", dQuote(outcome), " the testScale (passed as an argument to scoreTest) must have a single row for this subscale.")
      }
    } else {
      ts <- data.frame(scale=1, locaiton=0)
    }
  }
  scale <- ts$scale
  location <- ts$location

  if(is.na(scale) | is.na(location)) {
    if(is.na(scale)) {
      warning(paste0("Could not find a valid scale. Resetting to theta scale. Check ", dQuote("testDat"), " argument."))
      scale <- 1
    } else {
      warning(paste0("Could not find a valid location. Resetting to theta scale. Check ", dQuote("testDat"), " argument."))
      location <- 0
    }
  }
  return(list(location=location, scale=scale))
}


