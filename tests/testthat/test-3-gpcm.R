context("GPCM") 
# setup GPCM test
require(testthat)

set.seed(142857)
n <- 2000
theta <- rnorm(n)
x1 <- runif(n)
theta <- theta + x1 * 0.2

csv2paramList <- function(pModel = c('GPCM', 'GRM', 'PCM'), dat, subsetItems = NULL){
  # If converting modelParameters tab from config, only grab
  # relevant variables and rename ItemID to ITEM_ID
  if (all(c('MODEL', 'ItemID', 'ScorePoints', paste0('P', 0:10)) %in% names(dat))) {
    dat <- dat[, c('ItemID', 'MODEL', 'ScorePoints', paste0('P', 0:10))]
    names(dat)[1] <- 'ITEM_ID'
  }
  
  params2use <- with(dat, dat[order(dat$ITEM_ID), ])
  
  if(!is.null(subsetItems)){
    dat <- dat[which(dat$ITEM_ID %in% subsetItems), ]
    params2use <- dat[match(subsetItems, dat$ITEM_ID), ] #sort according to item order input
  }
  
  pModel <- toupper(pModel)
  
  names(params2use) <- tolower(names(params2use))
  finalItems <- params2use$item_id
  
  # any polytymous items?-- this will take care of IRTGPC sp = 2, IRTGPC sp = 3, IRTPCL sp = 2, IRTPCL sp = 3
  any.poly <- any(params2use$scorepoints > 1)
  # identifly all other items-- IRTGPC sp = 1, IRT3pln sp = 1, IRTPCL sp = 1, IRT3pl sp = 1
  # in config file, IRTPCL scorepoints = 1 items have b parameter stored in p0, unlike other 3pl items
  # for example, parameters for IRT3pl sp = 1 are stored as p0 = a(1), p1 = b, p2 = c(0)
  pos <- which(params2use$scorepoints == 1 & grepl('3pl', tolower(params2use$model))) # 3pl
  pos1 <- which(params2use$scorepoints == 1 & tolower(params2use$model) %in% 'irtgpc') # gpc
  pos <- c(pos, pos1)
  pos2 <- which(params2use$scorepoints == 1 & tolower(params2use$model) %in% 'irtpcl') # pcl
  
  mcParams <- params2use[pos, paste0('p', 0:2)]
  mcParams$p2 <- with(mcParams, ifelse(is.na(p2), 0, p2)) # make c-parameter if it is missing, this is needed for IRTGPC sp = 1
  
  mcParams2 <- params2use[pos2, paste0('p', 0:2)]
  if(nrow(mcParams2) > 0){
    mcParams2$p1 <- mcParams2$p0 # b is stored in p0 for IRTPCL sp = 1 items
    mcParams2$p0 <- 1 # set a to 1
    mcParams2$p2 <- 0 # set c to 0
  }
  
  mcParams <- rbind(mcParams, mcParams2)
  binary <- c(pos, pos2)
  
  ## inserted 4 lines for modifying
  binaryID <- params2use$item_id[binary]  
  polyID <- params2use$item_id[-binary]
  
  if(any.poly){
    ### This is the code to get the CR items organized for the params list
    numSteps <- length( grep("[p][0-9]", colnames(params2use)) )
    if(pModel == 'PCM') steps <- paste('p', 0:(numSteps-1), sep='')
    if(pModel %in% c('GRM', 'GPCM')) steps <- paste('p', 1:(numSteps-1), sep='')
    if(length(binary) > 0){
      crParams <- params2use[-binary, steps]
    } else {
      crParams <- params2use[, steps]
    } 
    if(pModel == 'PCM') aVar <- rep(1, nrow(crParams))
    if(length(binary) > 0){
      if(pModel %in% c('GRM', 'GPCM')) aVar <- params2use[-binary, 'p0']
    } else{
      if(pModel %in% c('GRM', 'GPCM')) aVar <- params2use[, 'p0']
    }
    crParams$base <- 0
    if(pModel == 'PCM') crParams <- crParams[, c('base', steps)]
    if(pModel == 'GPCM') crParams <- crParams[, c('base', steps)]
    if(pModel == 'GRM') crParams <- crParams[, steps]
    crParams <- as.data.frame(t(crParams))
    crList <- as.list(crParams)
    crList <- lapply(crList, function(x) x[!is.na(x)])
    if(pModel == 'PCM') pModel <- 'GPCM' # so that model can be automatically called by irt.ability
    if(length(binary) > 0){ # test has both CM and polytomous items
      params <- list('3pl' = list(a = mcParams$p0, b = mcParams$p1, c = mcParams$p2), 
                     gpcm = list(a = aVar , d = crList, polyID=polyID), model = pModel, items = finalItems, mcID=binaryID, pyID=polyID, mcPos = binary, 
                     score.pts = params2use$scorepoints)
    } else { # test has only polytomous items
      params <- list('3pl' = NULL, 
                     gpcm = list(a = aVar , d = crList), model = pModel, items = finalItems, mcID=NULL, pyID=polyID, mcPos = NULL, score.pts = params2use$scorepoints)
    }
  } else { # test has only MC items
    if(pModel == 'PCM') pModel <- 'GPCM' # so that model can be automatically called by irt.ability
    params <- list('3pl' = list(a = mcParams$p0, b = mcParams$p1, c = mcParams$p2), 
                   gpcm = NULL, items = finalItems, mcID=binaryID, pyID=NULL, mcPos = binary, score.pts = params2use$scorepoints, model = pModel)
  }
  class(params) <- 'param.list'
  params     
}

gen <- function(paramTab, theta, key) {
  gpcm <- function (theta, d, score, a) {
    if (score > length(d)) {
      stop ("Score higher than maximum")
    }
    Da <- con$D * a
    exp(sum(Da * (theta - d[1:score]))) / sum(exp(cumsum(Da * (theta - d))))
  }
  params <- csv2paramList(pModel='GPCM', paramTab, key)
  # get parameter values from items
  a <- params$'3pl'$a
  b <- params$'3pl'$b
  c <- params$'3pl'$c
  resb <- t(sapply(theta, function(thetai) {
    rbinom(length(a),1, c + (1 - c) / (1 + exp(-1.7 * a * (thetai - b))))
  }))
  aa <- params$gpcm$a
  d <- params$gpcm$d
  ind.dichot <- params$mcPos
  #x2 <- x[!1:length(x) %in% ind.dichot] +1 
  con <- list(D=1.7)
  coli <- ncol(resb) +1
  for(ii in 1:(length(aa))) {
    resb <- cbind(resb, rep(NA, nrow(resb)))
    for(thetai in 1:length(theta)) {
      pr <- rep(0,length(d[[ii]]))
      for(ri in 1:length(d[[ii]])) {
        pr[ri] <- gpcm(theta[thetai], d[[ii]], ri, aa)
      }
      resb[thetai, coli] <- sample(-1+1:length(pr), size=1, prob=pr)
    }
    coli <- coli + 1
  }
  colnames(resb) <- key
  return(resb)
}
# paramTab
paramTab <- structure(list(ItemID = structure(1:14,
                                              .Label = c("m017401", "m017701", "m017901", "m018201", "m018401", "m018501", "m018601", "m020001", "m020501", "m046301", "m046501", "m051501", "n202831", "m073601"), class = "factor"),
                      P0 = c(0.25 , 1    , 1.15 , 0.52 , 1.11 , 1.64 , 0.78 , 0.72 , 0.72 , 0.89 , 0.92 , 1.2  , 0.75 , 1.13),
                      P1 = c(-5.16, -1.01, -0.93, -1.21, -1.03, 0.34 , 0.9  , -0.49, -0.62, -1.07, -0.23, 1.22 , -2.58, 1.46),
                      P2 = c(0.19 , 0.16 , 0.15 , 0.03 , 0.24 , 0.26 , 0.12 , 0    , 0    , 0.28 , 0.33 , 0.2  , 0.25 , 1.50),
                      P3 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , 0.60),
                      P4 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA),
                      P5 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA),
                      P6 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA),
                      P7 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA),
                      P8 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA),
                      P9 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA), 
                      P10 = c(NA  , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA),
                      ScorePoints = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 4L),
                      MODEL = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L), .Label = c("3pl", "GPCM"), class = "factor")), class = "data.frame", row.names = c(NA, -14L))
mat <- gen(paramTab, theta, key=paramTab$ItemID)
colnames(mat) <- c("m017401", "m017701", "m017901", "m018201", "m018401", "m018501", 
"m018601", "m020001", "m020501", "m046301", "m046501", "m051501", 
"n202831", "m073601")

dichotParamTab <- structure(list(ItemID = structure(1:13,
                                              .Label = c("m017401", "m017701", "m017901", "m018201", "m018401", "m018501", "m018601", "m020001", "m020501", "m046301", "m046501", "m051501", "n202831"), class = "factor"),
                      slope = c(0.25 , 1    , 1.15 , 0.52 , 1.11 , 1.64 , 0.78 , 0.72 , 0.72 , 0.89 , 0.92 , 1.2  , 0.75 ),
                      difficulty = c(-5.16, -1.01, -0.93, -1.21, -1.03, 0.34 , 0.9  , -0.49, -0.62, -1.07, -0.23, 1.22 , -2.58),
                      guessing = c(0.19 , 0.16 , 0.15 , 0.03 , 0.24 , 0.26 , 0.12 , 0    , 0    , 0.28 , 0.33 , 0.2  , 0.25 ),
                      scorePoints = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
                      D=rep(1.7,13),
                      MODEL = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("3pl", "GPCM"), class = "factor")), class = "data.frame", row.names = c(NA, -13L))
dichotParamTab$test <- "overall"

polyParamTab <- structure(list(ItemID = structure(1,
                                              .Label = c("m073601"), class = "factor"),
                      slope = c( 0.58 ),
                      itemLocation = c(1.13),
                      d0 = c(0),
                      d1 = c(-0.33),
                      d2 = c(-0.39),
                      d3 = c( 0.72),
                      D = 1.7,
                      scorePoints = c(3L),
                      MODEL = structure(c(2L), .Label = c("3pl", "GPCM"), class = "factor")), class = "data.frame", row.names = c(NA, -1L))

polyParamTab$test <- "overall"

testDat <- data.frame(test="overall",
                      location=277.1563,
                      scale=37.7297)

mat <- data.frame(mat)
rownames(mat) <- paste0("pseudo-student",1:nrow(mat))
mat$origwt <- 1
nperstratum <- 10
nstrata <- length(theta)/nperstratum
mat$repgrp1 <- rep(1:nstrata, each=nperstratum)
mat$jkunit <- rep(rep(1:2, each=nperstratum/2), nstrata)
stuItems <- mat[,1:14]
stuItems$oppID <- factor(rownames(mat), levels=rownames(mat))
stuItems <- reshape(data=stuItems, varying=c(paramTab$ItemID), idvar=c("oppID"), direction="long", v.names="score", times=paramTab$ItemID, timevar="key")
stuDat <- mat[, c('origwt', 'repgrp1', 'jkunit', as.character(paramTab$ItemID))]
stuDat$oppID <- rownames(stuDat)
mat$x1 <- stuDat$x1 <- x1
stuDat$origwt <- mat$origwt <- runif(nrow(stuDat)) * 4 * abs(stuDat$x1 + 3)
############### test functions ###############
context("expect to error")
test_that("expect to error", {
  stuDat_bad <- stuDat
  stuDat_bad$m073601[1] <- 5
  expect_error(capture.output(suppressWarnings(mml1 <- mml(formula = overall ~ 1, stuDat=stuDat_bad, 
                                            dichotParamTab=dichotParamTab, polyParamTab=polyParamTab, 
                                            Q=34, idVar="oppID", multiCore=FALSE, testScale=testDat))), "inconsistent with expectations")
  stuDat_bad <- stuDat
  stuDat_bad$m020501[1] <- 2
  expect_error(capture.output(suppressWarnings(mml1 <- mml(stuDat=stuDat_bad, 
                                                           dichotParamTab=dichotParamTab, polyParamTab=polyParamTab, 
                                                           Q=34, idVar="oppID", multiCore=FALSE, testScale=testDat, verbose=0))), "inconsistent with expectations")
})


context("mml agrees with scoreTest + mmlu")
test_that("mml agrees with scoreTest + mmlu", {
  suppressWarnings(mml1 <- mml(formula = overall ~ 1, stuDat=stuDat, 
                               dichotParamTab=dichotParamTab, 
                               polyParamTab=polyParamTab, Q=34, idVar="oppID", testScale=testDat))
  stuDat <- cbind(stuDat, mat[,1:14])
  scoreTest <<- getFromNamespace("scoreTest", "Dire")
  mmlu <<- getFromNamespace("mmlu", "Dire")
  scored_test <- scoreTest(stuDat=stuDat,
                           dichotParamTab=dichotParamTab,
                           polyParamTab=polyParamTab,
                           testScale=testDat,
                           idVar="oppID",
                           strataVar="repgrp1",
                           PSUVar="jkunit",
                           minNode = -5, maxNode = 5)
  
  mmlu1 <- mmlu(overall ~ 1, scored_test)
  expect_equal(coef(mml1), coef(mmlu1))
})


