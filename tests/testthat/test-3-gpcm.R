skip_on_cran()
context("GPCM") 
# setup GPCM test
require(testthat)
library(doParallel)
# known bug on mac-os for parallel
# This is a known workaround: https://github.com/rstudio/rstudio/issues/6692
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}

cl <- makeCluster(2, outfile="")
registerDoParallel(cl, cores=2)
set.seed(142857)
n <- 2000
theta <- rnorm(n)
x1 <- runif(n)
theta <- theta + x1 * 0.2
mean(theta)
gen <- function(paramTab, theta, key) {
  gpcm <- function (theta, d, score, a) {
    if (score > length(d)) {
      stop ("Score higher than maximum")
    }
    Da <- con$D * a
    exp(sum(Da * (theta - d[1:score]))) / sum(exp(cumsum(Da * (theta - d))))
  }
  params <- Dire:::csv2paramList(pModel='GPCM', paramTab, key)
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
                      P0 = c(0.25 , 1    , 1.15 , 0.52 , 1.11 , 1.64 , 0.78 , 0.72 , 0.72 , 0.89 , 0.92 , 1.2  , 0.75 , 0.58     ),
                      P1 = c(-5.16, -1.01, -0.93, -1.21, -1.03, 0.34 , 0.9  , -0.49, -0.62, -1.07, -0.23, 1.22 , -2.58, 1.14-0.10),
                      P2 = c(0.19 , 0.16 , 0.15 , 0.03 , 0.24 , 0.26 , 0.12 , 0    , 0    , 0.28 , 0.33 , 0.2  , 0.25 , 1.14-0.16),
                      P3 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , 1.14-0.06),
                      P4 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , 1.14+0.32),
                      P5 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA),
                      P6 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA),
                      P7 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA),
                      P8 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA),
                      P9 = c(NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA), 
                      P10 = c(NA  , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA   , NA),
                      ScorePoints = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 5L),
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
                      d1 = c(0.13),
                      d2 = c( 0.24),
                      d3 = c( -0.38),
                      d4 = c( 0.01),
                      D = 1.7,
                      scorePoints = c(5L),
                      MODEL = structure(c(2L), .Label = c("3pl", "GPCM"), class = "factor")), class = "data.frame", row.names = c(NA, -1L))

polyParamTab$test <- "overall"

testDat <- data.frame(location=c(277.1563),
                      scale=c(37.7297))

# for comparing to NAEP
gpcm2 <- function (theta, b, d, score, a) {
  if (score > length(d)) {
    stop ("Score higher than maximum")
  }
  Da <- con$D * a
  exp(sum(Da * (theta - b + d[1:score]))) / sum(exp(cumsum(Da * (theta - b + d))))
}

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
stuDat <- mat[, c('origwt', 'repgrp1', 'jkunit')]
stuDat$oppID <- rownames(stuDat)
mat$x1 <- stuDat$x1 <- x1
stuDat$origwt <- mat$origwt <- runif(nrow(stuDat)) * 4 * abs(stuDat$x1 + 3)
############### test functions ###############
# write here to compare to AM
# library(haven)
# write_sav(mat,'Q:/Direct Estimation/QC/math 2003 g8 number subtest/math2003g8PBRandom2k_x1_poly.sav')
# tests:
mml1 <- mml(stuItems=split(stuItems, stuItems$oppID), stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab, Q=34, idVar="oppID", multiCore=FALSE, testScale=testDat)
mml1Taylor <- summary(mml1, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit", gradientHessian=TRUE)

expect_equal(coef(mml1Taylor)[,1], c(`(Intercept)` =  281.432554746317, `Population SD` = 36.387298620653), tolerance=20*sqrt(.Machine$double.eps))
     
expect_equal(coef(mml1Taylor)[,2], c(`(Intercept)` = 0.962322764662336, `Population SD` = 0.85340608987368), tolerance=5*(.Machine$double.eps)^0.25)
mml1mc <- mml(stuItems=split(stuItems, stuItems$oppID), stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab, Q=34, idVar="oppID", multiCore=TRUE, testScale=testDat)
mml1mcTaylor <- summary(mml1mc, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit", gradientHessian=TRUE)
expect_equal(coef(mml1), coef(mml1mc), tolerance=20*sqrt(.Machine$double.eps))
expect_equal(coef(mml1Taylor)[,2], coef(mml1mcTaylor)[,2], tolerance=20*sqrt(.Machine$double.eps)) 

mml2 <- mml(~x1,stuItems=split(stuItems, stuItems$oppID), stuDat=stuDat,  dichotParamTab=dichotParamTab, polyParamTab=polyParamTab, Q=34, idVar="oppID", weightVar="origwt", multiCore=TRUE, testScale=testDat)
mml2Taylor <- summary(mml2, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit", gradientHessian=TRUE)

expect_equal(coef(mml2Taylor)[,1], c(`(Intercept)` = 276.190834717205, x1 = 10.4808589128737, `Population SD` = 35.9063749547461), tolerance=20*sqrt(.Machine$double.eps))

expect_equal(coef(mml2Taylor)[,2], c(`(Intercept)` = 2.0454799919857, x1 = 3.45532352630566, `Population SD` = 0.95252258578783), tolerance=2*(.Machine$double.eps)^0.25)

stopCluster(cl)