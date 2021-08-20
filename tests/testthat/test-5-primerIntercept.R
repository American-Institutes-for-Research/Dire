# AM Results updated on 10/15/2020
skip_on_cran()
context("NAEPPrimer composite")
require(EdSurvey)

## Dichotomous Parameter Tab 
# paramTab$test
# paramTab$subtest
# paramTab$a/slope <-
# paramTab$g/asymptote/guessing <- 
# paramTab$d/difficulty <- 
# param$D <- 1.7
dichotParamTab <- structure(list(ItemID = structure(1:11,
                                                    .Label = c("m109801", "m020001", "m111001", "m046301", "m046501", "m051501", "m111601", "m111301", "m111201", "m110801", "m110101"), class = "factor"),
                                 test = rep("composite",11),
                                 subtest = c(rep("num",6),rep("alg",5)),
                                 slope = c(0.96, 0.69, 0.83, 0.99, 1.03, 0.97, 1.45, 0.59, 0.34, 0.18, 1.20),
                                 difficulty = c(-0.37, -0.55, 0.85 , -0.97, -0.14, 1.21, 0.53, -1.84, -0.46, 2.43, 0.70),
                                 guessing = c(0.16,0.00,0.17,0.31,0.37,0.18, 0.28, 0.15, 0.09, 0.05, 0.18),
                                 D = rep(1.7, 11),
                                 MODEL = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "3pl", class = "factor")), class = "data.frame", row.names = c(NA, -11L))


## GPCM Parameter Tab 
# a / slope
# d1,2,...,(param-1)
# b/itemLocation, then map to d1, d2
# D / D
polyParamTab <- structure(list(ItemID = structure(c(1:2), .Label = c("m0757cl", "m066501"), class = "factor"),
                               test = rep("composite",2),
                               subtest = c(rep("alg",2)),
                               slope = c( 0.43, 0.52 ), # could also be called "a"
                               itemLocation = c( -1.21, -0.96), # could also not exist, added to d1 ... dn
                               d1 = c(2.38, -0.56), 
                               d2 = c(-0.57, 0.56),
                               d3 = c(-1.18, NA),
                               D = c(1.7, 1.7),
                               scorePoints = c(4L, 3L) # number of score points, read d1 to d(n-1)
), class = "data.frame", row.names=c(1L,2L))

## Test Data information
testDat <- data.frame(
  subtest=c("num", "alg"),
  location=c(277.1563, 280.2948),
  scale=c(37.7297, 36.3887),
  subtestWeight=c(0.3,0.7),
  test=c("composite","composite"))

## Read-in NAEP Primer and Subset to items in Parameter Tabs 
sdf <- readNAEP(system.file("extdata/data", "M36NT2PM.dat", package = "NAEPprimer"))
items <- c(as.character(dichotParamTab$ItemID), as.character(polyParamTab$ItemID))
edf <- getData(data=sdf, varnames=c(items, "dsex", "origwt", "repgrp1", "jkunit"), omittedLevels = FALSE)
edf$sid <- paste0("S",1:nrow(edf))

# apply super basic scoring to Primer 
for(i in 1:length(items)) {
  coli <- items[i]
  # save the original
  edf[,paste0(coli,"raw")] <- edf[,coli]
  if( coli %in% dichotParamTab$ItemID) {
    edf[,coli] <- ifelse(grepl("[ABCDE]", edf[,paste0(coli,"raw")]), 0, NA)
    edf[,coli] <- ifelse(grepl("*", edf[,paste0(coli,"raw")], fixed=TRUE), 1, edf[,coli])
  } else {
    edf[,coli] <- ifelse(grepl("Incorrect", edf[,paste0(coli,"raw")], fixed=TRUE), 0, NA)
    edf[,coli] <- ifelse(grepl("None correct", edf[,paste0(coli,"raw")], fixed=TRUE), 0, edf[,coli])
    edf[,coli] <- ifelse(grepl("Incorrect", edf[,paste0(coli,"raw")], fixed=TRUE), 0, edf[,coli])
    edf[,coli] <- ifelse(grepl("One correct", edf[,paste0(coli,"raw")], fixed=TRUE), 1, edf[,coli])
    edf[,coli] <- ifelse(grepl("Partial", edf[,paste0(coli,"raw")], fixed=TRUE), 1, edf[,coli])
    edf[,coli] <- ifelse(grepl("Two correct", edf[,paste0(coli,"raw")], fixed=TRUE), 2, edf[,coli])
    edf[,coli] <- ifelse(grepl("Correct", edf[,paste0(coli,"raw")], fixed=TRUE), 2, edf[,coli])
    edf[,coli] <- ifelse(grepl("Three correct", edf[,paste0(coli,"raw")], fixed=TRUE), 3, edf[,coli])
  }
  # print results
  #print(table(edf[,paste0(coli,"raw")], edf[,coli], useNA="always"))
  edf[,paste0(coli,"raw")] <- NULL # delete original
}

# Separate items and student data from eachother 
stuItems <- edf[,c("sid", items)]
stuDat <- edf[,c("sid", "dsex", "origwt", "repgrp1", "jkunit", getWeightJkReplicates(data=sdf, "origwt"))]

testDat <- data.frame(
  subtest=c("num", "alg"),
  location=c(277.1563, 280.2948),
  scale=c(37.7297, 36.3887),
  subtestWeight=c(0.3,0.7),
  test=c("composite","composite"))

# make stuItems long
stuItems <- reshape(data=stuItems, 
                    varying=c(items), idvar=c("sid"), direction="long", v.names="score", times=items, timevar="key")

if(FALSE) {
  ## Run to output direct estimation ready data to AM ready data (dct and sav file)
  saveSav(stuItems, edf, "IRTparams", getwd(), idVar="sid")
  paramTab2AMdct(dichotParamTab,
                 polyParamTab,
                 testDat,
                 weightVar="origwt",
                 stratumVar="repgrp1",
                 PSUVar="jkunit",
                 testName="IRTparams",
                 saveFilePath=getwd())
}

# algebral score
mmlA <- mml(alg ~ 1, stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
            Q=34, idVar="sid", testScale=testDat, composite=FALSE, weightVar="origwt",
            minNode = -5, maxNode = 5, bobyqaControl=list(maxfun=1e5),
            strataVar="repgrp1", PSUVar="jkunit")
# Eric Buehler updated results on 10/15/2020
mmlATaylor <- summary(mmlA, varType="Taylor")

# AM output on 10/21/2020
#  Observations:        8826
#  Strata Variable:     REPGRP1
#  Cluster Variable:    JKUNIT
#  Weight Variable:     ORIGWT
# Iterations: 7 
# Log Likelihood: -16394.5  
# Dependent Variable: composite, alg  
# --  Weighted N  Estimate  Standard Error  
# All 8782.713764235  278.334734151 1.035753180 
# Root MSE  ... 32.173725366  1.132270387 
# AM Statistical Software Beta Version 0.06.04. (c) The American Institutes for Research and Jon Cohen  

expect_equal(unname(mmlATaylor$coef[,1]), c(278.334734151, 32.173725366), tolerance=25*sqrt(.Machine$double.eps))
expect_equal(unname(mmlATaylor$coef[1,2]), c(1.035753180), tolerance=5*(.Machine$double.eps)^0.25)
expect_equal(mmlA$LogLik, -16394.5, tol=1E-5)
# check that replicate weights approximately agree with Taylor series estimates
repw <- colnames(stuDat)[grep("srwt", colnames(stuDat))]
mmlARep <- summary(mmlA, varType = "replicate", repWeight = repw)
expect_equal(mmlATaylor$coef[,1], mmlARep$coef[,1], tolerance=sqrt(.Machine$double.eps))
expect_equal(mmlATaylor$coef[,2], mmlARep$coef[,2], tolerance=0.002)

# numeracy score
mmlN <- mml(num ~ 1, stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
            Q=34, idVar="sid", testScale=testDat, composite=FALSE, weightVar="origwt",
            minNode = -5, maxNode = 5, bobyqaControl=list(maxfun=1e5),
            strataVar="repgrp1", PSUVar="jkunit")

# AM output, 10/21/2020
#  Observations:        6378
#  Strata Variable:     REPGRP1
#  Cluster Variable:    JKUNIT
#  Weight Variable:     ORIGWT
# Iterations: 8 
# Log Likelihood: -10299.1  
# Dependent Variable: composite, num  
# --  Weighted N  Estimate  Standard Error  
# All 6254.724474587  285.351404123 1.244369007 
# Root MSE  ... 32.354923506  1.366182380 
# AM Statistical Software Beta Version 0.06.04. (c) The American Institutes for Research and Jon Cohen  
# Eric Buehler updated results on 10/15/2020
mmlNTaylor <- summary(mmlN, varType="Taylor")
repw <- colnames(stuDat)[grep("srwt", colnames(stuDat))]
#mmlNRep <- summary(mmlN, varType = "replicate", repWeight = repw)
expect_equal(unname(mmlNTaylor$coef[,1]), c(285.351403043, 32.354933129), tolerance=25*sqrt(.Machine$double.eps))
expect_equal(unname(mmlNTaylor$coef[1,2]), c(1.244369241), tolerance=0.004)
expect_equal(mmlN$LogLik, -10299.1, tol=1e-5)

#expect_equal(mmlNTaylor$coef[,1], mmlNRep$coef[,1], tolerance=25*sqrt(.Machine$double.eps))
#expect_equal(mmlNTaylor$coef[,2], mmlNRep$coef[,2], tolerance=0.002)


if(FALSE){
  library(doParallel)
  if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
      Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
    parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
  }
    
    cl <- makeCluster(4, outfile="")
    registerDoParallel(cl, cores=4)
    clusterExport(cl, c("fnCor"))
    set.seed(142857)
    
    mmlC <- mml(composite ~ 1, stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
                Q=34, idVar="sid", testScale=testDat, composite=TRUE, weightVar="origwt",
                minNode = -5, maxNode = 5, bobyqaControl=list(maxfun=1e5),
                multiCore=TRUE,
                strataVar="repgrp1", PSUVar="jkunit")
} else {
  # composite
  mmlC <- mml(composite ~ 1, stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
              Q=34, idVar="sid", testScale=testDat, composite=TRUE, weightVar="origwt",
              minNode = -5, maxNode = 5, bobyqaControl=list(maxfun=1e5),
              strataVar="repgrp1", PSUVar="jkunit", fast=TRUE)
}
#mmlCRep <- summary(mmlC, varType = "replicate", repWeight = repw)
# check that composite coefficients really are composite
# only [1] because the SD is not formed in this way
expect_equal(0.3* coef(mmlN)[1] + 0.7*coef(mmlA)[1], coef(mmlC)[1])

# AM run 10/21/2020
#  Observations:        16915
#  Strata Variable:     REPGRP1
#  Cluster Variable:    JKUNIT
#  Weight Variable:     ORIGWT  
# Composite Results 
# --  Weighted N  Mean  Standard Error  Standard Deviation  
# All 16932.463329904 280.439735226 1.098169602 32.223451497  
# use AM means
# means correctly estimates correlation, AM's regression finds a false maximum for the correlation
expect_equal(unname(coef(mmlC)[1]), c(280.4397365300, 32.223451497)[1], tolerance=1e-6)
expect_equal(unname(coef(mmlC)[2]), c(280.4397365300, 32.223451497)[2], tolerance=0.001)

context("NAEPPrimer composite variance")

# Taylor variance
mmlCTaylor <- summary(mmlC, varType="Taylor", gradientHessian=TRUE)
# updated on 10/16/2020 to mml means (which gets the correct correlation)
expect_equal(mmlCTaylor$coef[1,2], c(1.098169521), 
             tolerance=20*(.Machine$double.eps)^0.25)

# Replicate variance, regression only 
#expect_equal(mmlCRep$coef[1, 2], c(0.978889737365143), 
#             tolerance=200*sqrt(.Machine$double.eps))

#expect_equal(mmlCTaylor$coef[,1], mmlCRep$coef[,1], tolerance=25*sqrt(.Machine$double.eps))
# this is much less precise
#expect_equal(mmlCTaylor$coef[1,2], mmlCRep$coef[1,2], tolerance=0.15)
